import gzip
import os
import sys
import shutil
import traceback
import time
import numpy as np
from numba import njit
from typing import Literal

from structman.base_utils.config_class import Config
from structman.base_utils.base_utils import is_alphafold_model, Residue_Map, Slotted_obj
from structman.lib.rin_utils import calcRIN, get_entry_path, call_in_out_flow
from structman.lib.sdsc.consts import residues, ligands
from structman.lib.database.afdb_retrieval import get_afdb_rin_and_cent


class Interaction_type(Slotted_obj):
    def addEdge(self, interaction_type: str, score: float) -> None:
        # The leaf of the edge dispatch, called once per interaction of every RIN
        # file, so it neither rebuilds the pair nor walks the exception path when
        # the slot already holds one. setScore mutates the same list in place, so
        # updating rather than replacing it is the behaviour that was already there.
        try:
            interaction_data = self.__getattribute__(interaction_type)
        except AttributeError:
            interaction_data = None
        if type(interaction_data) is list:
            interaction_data[0] += 1
            interaction_data[1] += score
            return
        try:
            self.__setattr__(interaction_type, [1, score])
        except AttributeError:
            raise KeyError(f'Unknown interactiontype: {interaction_type}')

    def __init__(self):
        #...
        for slot in self.__slots__:
            self.__setattr__(slot, None)

    def to_list(self) -> list[int | float]:
        data_list: list[int | float] = []
        for slot in self.__slots__:
            try:
                dt = list(self.__getattribute__(slot))
            except (AttributeError, TypeError):
                dt = [float('nan'), float('nan')]
            try:
                data_list += dt
            except TypeError:
                data_list += [float('nan'), float('nan')]
        return data_list

    def set_data_from_list(self, value_list):
        left = 0
        #right = 2
        for slot in self.__slots__:
            try:
                deg = int(value_list[left])
            except ValueError:
                deg = 0
            dt = [deg, value_list[left+1]]
            self.__setattr__(slot, dt)
            left += 2
            #right += 2

    def getScore(self, interaction_type: str, debug_print=False) -> float:
        try:
            dt = self.__getattribute__(interaction_type)
            if debug_print:
                print(f'base level {interaction_type=} {dt=}')
            score = dt[1]
            if score is None:
                return np.nan
            return score
        except (AttributeError, TypeError, IndexError) as err:
            #print(f'{err=} {interaction_type=} {self.__getattribute__(interaction_type)=}')
            self.__setattr__(interaction_type, [0, np.nan])
            return self.__getattribute__(interaction_type)[1]
        
    def getDegree(self, interaction_type: str) -> int:
        try:
            degree = self.__getattribute__(interaction_type)[0]
            if np.isnan(degree) or degree is None:
                return 0
            return degree
        except (AttributeError, TypeError, IndexError):
            self.__setattr__(interaction_type, [0, np.nan])
            return self.__getattribute__(interaction_type)[0]

    def setScore(self, interaction_type: str, score: float) -> None:
        try:
            self.__getattribute__(interaction_type)[1] = score
        except (AttributeError, TypeError):
            self.__setattr__(interaction_type, [0, score])
            
    def setDegree(self, interaction_type: str, degree: int) -> None:
        try:
            self.__getattribute__(interaction_type)[0] = degree
        except (AttributeError, TypeError):
            self.__setattr__(interaction_type, [degree, np.nan])

class Ligand_types(Interaction_type):
    __slots__ = ['ligand', 'ion', 'metal']

ligand_types = set(['ligand', 'ion', 'metal'])

class Intrachain_types(Interaction_type):
    __slots__ = ['neighbor', 'short', 'long']

intrachain_types = set(['neighbor', 'short', 'long'])

class Interchain_types(Interaction_type):
    __slots__ = ['Protein', 'DNA', 'RNA', 'Peptide', 'Homomer']

interchain_types = set(['Protein', 'DNA', 'RNA', 'Peptide', 'Homomer'])

interaction_types = set(['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide', 'Homomer'])

non_intra_interaction_types = set(['ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide', 'Homomer'])

chain_type_interactions = set(['Protein', 'DNA', 'RNA', 'Peptide', 'Homomer'])

small_molecule_interactions = set(['ligand', 'ion', 'metal'])


# Which of Interaction_type_profile's three holders owns each interaction type.
_INTERACTION_HOLDER: dict[str, str] = {}
for _type in interchain_types:
    _INTERACTION_HOLDER[_type] = 'interchain'
for _type in intrachain_types:
    _INTERACTION_HOLDER[_type] = 'intrachain'
for _type in ligand_types:
    _INTERACTION_HOLDER[_type] = 'sw_ligand'


class Interaction_type_profile(Slotted_obj):
    __slots__ = ['sw_ligand', 'interchain', 'intrachain']

    def __init__(self):
        self.sw_ligand = Ligand_types()
        self.interchain = Interchain_types()
        self.intrachain = Intrachain_types()

    def to_list(self) -> list[int | float]:
        #3 +3 +5 = 11
        data_list = self.sw_ligand.to_list() + self.interchain.to_list() + self.intrachain.to_list()
        return data_list
    
    def set_data_from_list(self, value_list):
        self.sw_ligand.set_data_from_list(value_list[:6])
        self.interchain.set_data_from_list(value_list[6:16])
        self.intrachain.set_data_from_list(value_list[16:])

    def addEdge(self, interaction_type, score):
        # One dict lookup instead of walking up to three set memberships.
        holder = _INTERACTION_HOLDER.get(interaction_type)
        if holder is None:
            raise TypeError(f'Unknown interactiontype: {interaction_type}')
        return getattr(self, holder).addEdge(interaction_type, score)

    def getScore(self, interaction_type, debug_print=False):
        if debug_print:
            print(f'type profile level: {interaction_type=}')
        if interaction_type in interchain_types:
            return self.interchain.getScore(interaction_type, debug_print=debug_print)
        elif interaction_type in intrachain_types:
            return self.intrachain.getScore(interaction_type, debug_print=debug_print)
        elif interaction_type in ligand_types:
            return self.sw_ligand.getScore(interaction_type, debug_print=debug_print)
        else:
            return np.nan

    def getDegree(self, interaction_type):
        if interaction_type in interchain_types:
            return self.interchain.getDegree(interaction_type)
        elif interaction_type in intrachain_types:
            return self.intrachain.getDegree(interaction_type)
        elif interaction_type in ligand_types:
            return self.sw_ligand.getDegree(interaction_type)
        else:
            return 0

    def setScore(self, interaction_type, score):
        if interaction_type in interchain_types:
            return self.interchain.setScore(interaction_type, score)
        elif interaction_type in intrachain_types:
            return self.intrachain.setScore(interaction_type, score)
        elif interaction_type in ligand_types:
            return self.sw_ligand.setScore(interaction_type, score)
        else:
            return f'Unknown interactiontype: {interaction_type}'

    def setDegree(self, interaction_type, degree):
        if interaction_type in interchain_types:
            return self.interchain.setDegree(interaction_type, degree)
        elif interaction_type in intrachain_types:
            return self.intrachain.setDegree(interaction_type, degree)
        elif interaction_type in ligand_types:
            return self.sw_ligand.setDegree(interaction_type, degree)
        else:
            return f'Unknown interactiontype: {interaction_type}'


bond_types = {'cnt', 'hbond', 'ovl'}


class Chain_types(Slotted_obj):
    __slots__ = ['contact', 'h_bond', 'overlap']

    def __init__(self):
        self.contact = Interaction_type_profile()
        self.h_bond = Interaction_type_profile()
        self.overlap = Interaction_type_profile()

    def to_list(self) -> list[int | float]:
        #11 + 11 + 11 = 33
        data_list = self.contact.to_list() + self.h_bond.to_list() + self.overlap.to_list()
        return data_list

    def set_data_from_list(self, value_list):
        self.contact.set_data_from_list(value_list[:22])
        self.h_bond.set_data_from_list(value_list[22:44])
        self.overlap.set_data_from_list(value_list[44:])

    def addEdge(self, bondtype, interaction_type, score):
        if bondtype == 'cnt':
            return self.contact.addEdge(interaction_type, score)
        elif bondtype == 'hbond':
            return self.h_bond.addEdge(interaction_type, score)
        elif bondtype == 'ovl':
            return self.overlap.addEdge(interaction_type, score)
        else:
            return f'Unknown bondtype: {bondtype}'

    def getScore(self, bondtype, interaction_type, debug_print = False):
        if debug_print:
            print(f'chain level: {bondtype=} {interaction_type=}')
        if bondtype == 'cnt':
            return self.contact.getScore(interaction_type, debug_print=debug_print)
        elif bondtype == 'hbond':
            return self.h_bond.getScore(interaction_type, debug_print=debug_print)
        elif bondtype == 'ovl':
            return self.overlap.getScore(interaction_type, debug_print=debug_print)
        else:
            return f'Unknown bondtype: {bondtype}'

    def getDegree(self, bondtype, interaction_type):

        if bondtype == 'cnt':
            return self.contact.getDegree(interaction_type)
        elif bondtype == 'hbond':
            return self.h_bond.getDegree(interaction_type)
        elif bondtype == 'ovl':
            return self.overlap.getDegree(interaction_type)
        else:
            return f'Unknown bondtype: {bondtype}'

    def setScore(self, bondtype, interaction_type, score):

        if bondtype == 'cnt':
            return self.contact.setScore(interaction_type, score)
        elif bondtype == 'hbond':
            return self.h_bond.setScore(interaction_type, score)
        elif bondtype == 'ovl':
            return self.overlap.setScore(interaction_type, score)
        else:
            return f'Unknown bondtype: {bondtype}'

    def setDegree(self, bondtype, interaction_type, degree):

        if bondtype == 'cnt':
            return self.contact.setDegree(interaction_type, degree)
        elif bondtype == 'hbond':
            return self.h_bond.setDegree(interaction_type, degree)
        elif bondtype == 'ovl':
            return self.overlap.setDegree(interaction_type, degree)
        else:
            return f'Unknown bondtype: {bondtype}'

    def getCombiScore(self, interaction_type):
        ind_scores = np.array([
            self.contact.getScore(interaction_type),
            self.h_bond.getScore(interaction_type),
            self.overlap.getScore(interaction_type)])
        return np.nansum(ind_scores)

    def getCombiDegree(self, interaction_type):
        
        return self.contact.getDegree(interaction_type) + self.h_bond.getDegree(interaction_type) + self.overlap.getDegree(interaction_type)


chain_types = {'mc', 'sc'}


class Interaction_profile(Slotted_obj):
    __slots__ = [
        'mainchain', 'sidechain', '_class',
        '_simple_class', 'interacting_chains', 'interacting_ligands',
        'total_score'
        ]

    slot_mask = [True, True, False, False, False, False, False]

    def __init__(self):
        #self.mainchain = Chain_types()
        #self.sidechain = Chain_types()
        self.interacting_chains = set()
        self.interacting_ligands = set()
        self.total_score = None

    def to_list(self):
        try:
            m_list = self.mainchain.to_list()
        except AttributeError:
            m_list = [float('nan')] * 66
        try:
            s_list = self.sidechain.to_list()
        except AttributeError:
            s_list = [float('nan')] * 66

        return m_list + s_list

    def set_data_from_list(self, value_list):
        self.mainchain = Chain_types()
        self.mainchain.set_data_from_list(value_list[0:66])
        self.sidechain = Chain_types()
        self.sidechain.set_data_from_list(value_list[66:])

    def addEdge(self, chaintype, bondtype, interaction_type, score, chain_id, res_id):
        #print(f'Call of addEdge: {chaintype=} {bondtype=} {interaction_type=} {score=} {chain_id=} {res_id=}')
        if bondtype != 'ovl':
            if interaction_type in small_molecule_interactions:
                self.interacting_ligands.add(f'{chain_id}:{res_id}')
            elif interaction_type in chain_type_interactions:
                self.interacting_chains.add(chain_id)

        if chaintype == 'mc':
            try:
                return self.mainchain.addEdge(bondtype, interaction_type, score)
            except AttributeError:
                self.mainchain = Chain_types()
                return self.mainchain.addEdge(bondtype, interaction_type, score)
        elif chaintype == 'sc' or chaintype == 'ligand':
            try:
                return self.sidechain.addEdge(bondtype, interaction_type, score)
            except AttributeError:
                self.sidechain = Chain_types()
                return self.sidechain.addEdge(bondtype, interaction_type, score)
        else:
            return f'Unknown chaintype in addEdge: {chaintype}'

    def getScore(self, chain_type, bond_type, interaction_type, debug_print = False):
        if debug_print:
            print(f'profile level: {chain_type=} {bond_type=} {interaction_type=}')
        if chain_type == 'mc':
            try:
                return self.mainchain.getScore(bond_type, interaction_type, debug_print=debug_print)
            except AttributeError:
                self.mainchain = Chain_types()
                return self.mainchain.getScore(bond_type, interaction_type, debug_print=debug_print)
            
        elif chain_type == 'sc':
            try:
                return self.sidechain.getScore(bond_type, interaction_type, debug_print=debug_print)
            except AttributeError:
                self.sidechain = Chain_types()
                return self.sidechain.getScore(bond_type, interaction_type, debug_print=debug_print)
        else:
            return f'Unknown chaintype in getScore: {chain_type}'

    def getDegree(
            self,
            chain_type: Literal['sc', 'mc'],
            bond_type: Literal['cnt', 'hbond', 'ovl'], 
            interaction_type: Literal['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide', 'Homomer']
            ):
        if chain_type == 'mc':
            try:
                return self.mainchain.getDegree(bond_type, interaction_type)
            except AttributeError:
                self.mainchain = Chain_types()
                return self.mainchain.getDegree(bond_type, interaction_type)

        elif chain_type == 'sc':
            try:
                return self.sidechain.getDegree(bond_type, interaction_type)
            except AttributeError:
                self.sidechain = Chain_types()
                return self.sidechain.getDegree(bond_type, interaction_type)
        else:
            return f'Unknown chaintype in getDegree: {chain_type}'

    def setScore(self, chain_type, bond_type, interaction_type, score):
        if chain_type == 'mc':
            try:
                return self.mainchain.setScore(bond_type, interaction_type, score)
            except AttributeError:
                self.mainchain = Chain_types()
                return self.mainchain.setScore(bond_type, interaction_type, score)
        elif chain_type == 'sc':
            try:
                return self.sidechain.setScore(bond_type, interaction_type, score)
            except AttributeError:
                self.sidechain = Chain_types()
                return self.sidechain.setScore(bond_type, interaction_type, score)
        else:
            return f'Unknown chaintype in setScore: {chain_type}'

    def setDegree(self, chain_type, bond_type, interaction_type, degree):
        if chain_type == 'mc':
            try:
                return self.mainchain.setDegree(bond_type, interaction_type, degree)
            except AttributeError:
                self.mainchain = Chain_types()
                return self.mainchain.setDegree(bond_type, interaction_type, degree)
        elif chain_type == 'sc':
            try:
                return self.sidechain.setDegree(bond_type, interaction_type, degree)
            except AttributeError:
                self.sidechain = Chain_types()
                return self.sidechain.setDegree(bond_type, interaction_type, degree)
        else:
            return f'Unknown chaintype in setDegree: {chain_type}'

    def getChainSpecificCombiScore(self, chain_type, interaction_type):
        if chain_type == 'mc':
            try:
                return self.mainchain.getCombiScore(interaction_type)
            except AttributeError:
                self.mainchain = Chain_types()
                return self.mainchain.getCombiScore(interaction_type)
        elif chain_type == 'sc':
            try:
                return self.sidechain.getCombiScore(interaction_type)
            except AttributeError:
                self.sidechain = Chain_types()
                return self.sidechain.getCombiScore(interaction_type)
        else:
            return None

    def getTotalInterfaceInteractionScore(self, config=None):
        if self.total_score is None:
            score = 0.
            for chain_type in chain_types:
                for i_type in non_intra_interaction_types:
                    
                    s = self.getChainSpecificCombiScore(chain_type, i_type)
                    if config is not None:
                        config.logger.info(f'In getTotalInterfaceInteractionScore: {chain_type=} {i_type=} {s=}')

                    if s is not None and not isinstance(s, str) and not np.isnan(s):
                        score += s
            self.total_score = score
        return self.total_score

    def getChainSpecificCombiDegree(self, chain_type, interaction_type):
        if chain_type == 'mc':
            try:
                return self.mainchain.getCombiDegree(interaction_type)
            except AttributeError:
                self.mainchain = Chain_types()
                return self.mainchain.getCombiDegree(interaction_type)
        elif chain_type == 'sc':
            try:
                return self.sidechain.getCombiDegree(interaction_type)
            except AttributeError:
                self.sidechain = Chain_types()
                return self.sidechain.getCombiDegree(interaction_type)
        else:
            return None

    def computeClass(self, debug_print = False):
        if debug_print:
            print(f'Call of computeClass: {self.to_list()}')
        interactions = []
        b_types = ['cnt', 'hbond']#, 'ovl']  # ignore overlaps
        i_types = ['ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide', 'Homomer']
        for interaction_type in i_types:
            max_score_interaction = None
            max_score = 0.
            for chain_type in chain_types:
                for bond_type in b_types:
                    score = self.getScore(chain_type, bond_type, interaction_type)
                    if np.isnan(score):
                        continue

                    if debug_print and score > 0.:
                        score = self.getScore(chain_type, bond_type, interaction_type, debug_print=debug_print)
                        print(f'in compute class {interaction_type=} {chain_type=} {bond_type=} {score=}')

                    if score > max_score:
                        max_score_interaction = (chain_type, bond_type)
                        max_score = score
            if max_score_interaction is not None:
                (chain_type, bond_type) = max_score_interaction
                interactions.append((max_score, chain_type, bond_type, interaction_type))

        s = len(interactions)
        if s == 0:
            self._class = 'No interaction'
            self._simple_class = 'No interaction'
        elif s == 1:
            score, chain_type, bond_type, interaction_type = interactions[0]
            self._class = f'{chain_type} {bond_type} {interaction_type}'
            self._simple_class = f'{interaction_type} interaction'
        else:
            i_str_parts = []
            max_score = 0.
            for score, chain_type, bond_type, interaction_type in interactions:
                i_str_parts.append(f'{chain_type} {bond_type} {interaction_type}')
                if score > max_score:
                    max_score = score
                    self._simple_class = f'{interaction_type} interaction'
            self._class = 'Multi interaction: {}'.format(' and '.join(i_str_parts))

    def getClass(self, debug_print=False):
        try:
            self.computeClass(debug_print=debug_print)
            return self._class, self._simple_class
        except AttributeError as err:
            self._class = None
            self._simple_class = None
            if debug_print:
                print(f'{err=}')
            self.computeClass(debug_print=debug_print)
            return self._class, self._simple_class


# --- canonical encoded-profile layout ----------------------------------------
#
# An Interaction_profile round-trips losslessly through its 132-float encoding
# (set_data_from_list -> to_list is the identity), so the flat array is not a
# serialisation of the profile - it *is* the profile, and the 15-object tree is
# a transient on either side of it. Everything that only reads a profile can
# therefore work off the array and skip the tree entirely.
#
# PROFILE_INDEX maps each (chain_type, bond_type, interaction_type) to the flat
# index of its degree; the matching score lives at that index + 1. The table is
# read back out of the real set_data_from_list/getDegree/getScore path with a
# sentinel array, so it cannot drift from the class definitions above, and the
# import-time checks below fail loudly if the layout ever changes shape.

PROFILE_ENCODED_LENGTH: int = 132


def _build_profile_index_table() -> dict[tuple[str, str, str], int]:
    sentinel = list(range(PROFILE_ENCODED_LENGTH))
    probe = Interaction_profile()
    probe.set_data_from_list(sentinel)

    table: dict[tuple[str, str, str], int] = {}
    covered: set[int] = set()

    for chain_type in sorted(chain_types):
        for bond_type in sorted(bond_types):
            for interaction_type in sorted(interaction_types):
                degree = probe.getDegree(chain_type, bond_type, interaction_type)
                score = probe.getScore(chain_type, bond_type, interaction_type)
                try:
                    degree_index = int(degree)
                    score_index = int(score)
                except (TypeError, ValueError):
                    raise RuntimeError(
                        f'encoded profile layout probe failed for {chain_type=} {bond_type=} {interaction_type=}'
                    )

                if score_index != degree_index + 1:
                    raise RuntimeError(
                        f'encoded profile layout is not [degree, score] adjacent for '
                        f'{chain_type=} {bond_type=} {interaction_type=}: {degree_index=} {score_index=}'
                    )

                table[(chain_type, bond_type, interaction_type)] = degree_index
                covered.add(degree_index)
                covered.add(score_index)

    if len(covered) != PROFILE_ENCODED_LENGTH:
        raise RuntimeError(
            f'encoded profile layout does not tile its {PROFILE_ENCODED_LENGTH} values '
            f'exactly once: covered {len(covered)}'
        )

    return table


PROFILE_INDEX: dict[tuple[str, str, str], int] = _build_profile_index_table()


def profile_degree_index(chain_type: str, bond_type: str, interaction_type: str) -> int:
    return PROFILE_INDEX[(chain_type, bond_type, interaction_type)]


def profile_score_index(chain_type: str, bond_type: str, interaction_type: str) -> int:
    return PROFILE_INDEX[(chain_type, bond_type, interaction_type)] + 1


# --- the stored profile ------------------------------------------------------
#
# RIN_based_features.profile used to hold the whole object tree, built by
# set_data_from_list off a flat 132 value array and read back through four levels
# of string dispatch. That tree cost 69us per position to pack and unpack, and the
# 66 values add_to_output_object reads out of it cost 218us per SNV, none of which
# buys anything the flat form does not already carry. The slot holds the flat list
# now and the tree is only built where something still needs the object.
#
# set_data_from_list is not a plain copy, so the flat form is not the encoded array
# either: it truncates every degree with int(), and a NaN degree raises ValueError
# there and lands on 0. Normalising once on the way in keeps every reader below a
# bare index, and makes the stored list exactly what to_list() used to return.

# Each degree sits at an even index with its score directly after it. That follows
# from PROFILE_INDEX tiling all 132 values with 66 adjacent [degree, score] pairs,
# which _build_profile_index_table checks, but the pairing is what the loop below
# relies on, so it is spelled out.
assert all(index % 2 == 0 for index in PROFILE_INDEX.values()), 'profile degrees are not at even indices'


def stored_profile_from_encoded(encoded_profile) -> list:
    """The flat form Interaction_profile.to_list() produced, without the tree."""
    try:
        values = encoded_profile.tolist()
    except AttributeError:
        values = list(encoded_profile)

    for degree_position in range(0, PROFILE_ENCODED_LENGTH, 2):
        try:
            # int() and not trunc(): None still has to raise TypeError here, the way
            # set_data_from_list did, rather than being quietly read as a degree.
            values[degree_position] = int(values[degree_position])
        except ValueError:
            values[degree_position] = 0

    return values


# The three values add_to_output_object reads per (chain type, interaction type):
# the combi score, which getCombiScore nansums over cnt/hbond/ovl; the combi degree,
# which getCombiDegree adds over the same three; and the hbond score on its own.
# Names and index triples are both fixed, so they are laid out once here in the
# order the output object is filled.
PROFILE_OUTPUT_CHAIN_TYPES: tuple[str, ...] = ('mc', 'sc')
PROFILE_OUTPUT_INTERACTION_TYPES: tuple[str, ...] = (
    'neighbor', 'short', 'long', 'ligand', 'ion', 'metal',
    'Protein', 'DNA', 'RNA', 'Peptide', 'Homomer')
PROFILE_OUTPUT_BOND_TYPES: tuple[str, ...] = ('cnt', 'hbond', 'ovl')


def _build_profile_output_layout() -> tuple[tuple, ...]:
    layout = []
    for chain_type in PROFILE_OUTPUT_CHAIN_TYPES:
        for interaction_type in PROFILE_OUTPUT_INTERACTION_TYPES:
            degree_indices = tuple(
                profile_degree_index(chain_type, bond_type, interaction_type)
                for bond_type in PROFILE_OUTPUT_BOND_TYPES)
            layout.append((
                f'{chain_type} {interaction_type} score',
                f'{chain_type} {interaction_type} degree',
                f'{chain_type} {interaction_type} H-bond score',
                tuple(index + 1 for index in degree_indices),
                degree_indices,
                profile_score_index(chain_type, 'hbond', interaction_type),
            ))
    return tuple(layout)


PROFILE_OUTPUT_LAYOUT: tuple[tuple, ...] = _build_profile_output_layout()

# Same names in the same order, for the branch that has no profile to read.
PROFILE_OUTPUT_FEATURE_NAMES: tuple[str, ...] = tuple(
    name for entry in PROFILE_OUTPUT_LAYOUT for name in entry[:3])


def as_stored_profile(profile):
    """The stored flat form, accepting a profile written by an older version.

    Position_Data is a persisted column, so packed Mappings rows outlive the code
    that wrote them. Until the profile became a flat list it was the object tree,
    and rows already in a database still decode into one. to_list() is exactly what
    the stored form is defined to be, so converting here keeps those rows readable
    for one isinstance per row.
    """
    if isinstance(profile, Interaction_profile):
        return profile.to_list()
    return profile


@njit(cache=True)
def profile_degree_at(profile_row: np.ndarray, degree_index: int) -> float:
    # Mirrors the object path exactly: set_data_from_list truncates the degree with
    # int() when decoding, and getDegree maps a missing (NaN) slot to 0. The raw
    # encoded value is therefore NOT the degree - it has to be truncated here too.
    value = profile_row[degree_index]
    if np.isnan(value):
        return 0.0
    return np.trunc(value)


@njit(cache=True)
def profile_score_at(profile_row: np.ndarray, degree_index: int) -> float:
    return profile_row[degree_index + 1]


def profile_matrix_from_profiles(profiles) -> tuple[np.ndarray, np.ndarray]:
    # Materialises a sequence of Interaction_profile|None into the flat form that
    # njit kernels can consume: an (n, 132) matrix plus a validity mask marking
    # which rows came from a real profile.
    n_profiles = len(profiles)
    matrix = np.full((n_profiles, PROFILE_ENCODED_LENGTH), np.nan, dtype=np.float64)
    valid = np.zeros(n_profiles, dtype=bool)

    for row_index, profile in enumerate(profiles):
        if profile is None:
            continue
        matrix[row_index] = profile.to_list()
        valid[row_index] = True

    return matrix, valid


# getTotalInterfaceInteractionScore sums getCombiScore over both chain types and the
# eight non-intrachain interaction types, and getCombiScore is itself a nansum over
# cnt/hbond/ovl. That flattens to a plain nansum over these 48 score positions.
_INTERFACE_SCORE_INDICES: np.ndarray = np.array(
    sorted(
        profile_score_index(chain_type, bond_type, interaction_type)
        for chain_type in chain_types
        for bond_type in bond_types
        for interaction_type in non_intra_interaction_types
    ),
    dtype=np.int64,
)


@njit(cache=True)
def _total_interface_scores(profile_matrix: np.ndarray, valid: np.ndarray, score_indices: np.ndarray) -> np.ndarray:
    n_rows = profile_matrix.shape[0]
    n_indices = score_indices.shape[0]
    totals = np.zeros(n_rows, dtype=np.float64)

    for row in range(n_rows):
        if not valid[row]:
            continue
        total: float = 0.0
        for k in range(n_indices):
            value = profile_matrix[row, score_indices[k]]
            if not np.isnan(value):
                total += value
        totals[row] = total

    return totals


def total_interface_scores_from_matrix(profile_matrix: np.ndarray, valid: np.ndarray) -> np.ndarray:
    # Batch equivalent of Interaction_profile.getTotalInterfaceInteractionScore over a
    # whole chain. Summation order is fixed here, whereas the object path iterates the
    # chain_types / non_intra_interaction_types *sets*, so its order - and therefore its
    # last-bit rounding - already varies between processes.
    return _total_interface_scores(profile_matrix, valid, _INTERFACE_SCORE_INDICES)


# computeClass only ever reads the 'cnt' and 'hbond' scores of the eight
# interaction types below - 32 of the 132 encoded values, with the whole 'ovl'
# branch deliberately ignored. Decoding an encoded profile into the full
# 15-object tree just to serve those 32 lookups is the dominant cost of the
# aggregation hot loop, so the layout below lets the class be computed straight
# off the encoded array instead.
_CLASS_B_TYPES: list[str] = ['cnt', 'hbond']
_CLASS_I_TYPES: list[str] = ['ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide', 'Homomer']


def _build_class_score_layout() -> list[tuple[str, tuple[tuple[int, str, str], ...]]]:
    layout: list[tuple[str, tuple[tuple[int, str, str], ...]]] = []
    for interaction_type in _CLASS_I_TYPES:
        entries: list[tuple[int, str, str]] = []
        for chain_type in sorted(chain_types):
            for bond_type in _CLASS_B_TYPES:
                entries.append((profile_score_index(chain_type, bond_type, interaction_type), chain_type, bond_type))
        layout.append((interaction_type, tuple(entries)))
    return layout


_CLASS_SCORE_LAYOUT: list[tuple[str, tuple[tuple[int, str, str], ...]]] = _build_class_score_layout()


# --- simple class as an integer code -----------------------------------------
#
# The aggregation hot loop needs only the *simple* class of each mapping, and
# only to compare it against the position's aggregated simple class. Building
# the class strings per mapping meant several billion string constructions per
# chunk, so the loop now works with codes and the strings are materialised once
# per position. Codes >= 0 index _CLASS_I_TYPES; negative codes carry the
# location that 'No interaction' falls back to.
SIMPLE_CLASS_NO_LOCATION: int = -1
SIMPLE_CLASS_CORE: int = -2
SIMPLE_CLASS_BURIED: int = -3
SIMPLE_CLASS_SURFACE: int = -4

LOCATION_CODES: dict[str | None, int] = {
    None: SIMPLE_CLASS_NO_LOCATION,
    'Core': SIMPLE_CLASS_CORE,
    'Buried': SIMPLE_CLASS_BURIED,
    'Surface': SIMPLE_CLASS_SURFACE,
}

# Maps whatever computeClass()/classify() produced as a simple class string back
# onto the same code space, so the aggregated value can be compared with the
# per-mapping codes.
SIMPLE_CLASS_CODE_OF_STRING: dict[str | None, int] = dict(LOCATION_CODES)
for _type_index, _interaction_type in enumerate(_CLASS_I_TYPES):
    SIMPLE_CLASS_CODE_OF_STRING[f'{_interaction_type} interaction'] = _type_index

# Inverse of the above, for the few places that still need the string itself.
SIMPLE_CLASS_STRING_OF_CODE: dict[int, str | None] = {
    code: text for text, code in SIMPLE_CLASS_CODE_OF_STRING.items()
}


@njit(cache=True)
def locate_code(rsa: float, surface_threshold: float, buried_threshold: float) -> int:
    # locate() in code space. An encoded RSA is never None, and NaN fails both
    # comparisons and so lands on Core - exactly as the original chain of
    # comparisons did.
    if rsa > surface_threshold:
        return SIMPLE_CLASS_SURFACE
    if rsa > buried_threshold:
        return SIMPLE_CLASS_BURIED
    return SIMPLE_CLASS_CORE


def _build_class_score_index_matrix() -> np.ndarray:
    # (n_interaction_types x max_entries) of flat profile indices, -1 padded, so
    # the kernel below can walk the layout without python objects.
    widths = [len(entries) for _, entries in _CLASS_SCORE_LAYOUT]
    matrix = np.full((len(_CLASS_SCORE_LAYOUT), max(widths)), -1, dtype=np.int64)
    for type_index, (_, entries) in enumerate(_CLASS_SCORE_LAYOUT):
        for entry_index, (flat_index, _, _) in enumerate(entries):
            matrix[type_index, entry_index] = flat_index
    return matrix


_CLASS_SCORE_INDICES: np.ndarray = _build_class_score_index_matrix()


@njit(cache=True)
def _simple_class_code(profile: np.ndarray, score_indices: np.ndarray, location_code: int) -> int:
    # Mirrors compute_class_from_encoded()'s choice of simple class exactly:
    # per interaction type the strongest score above 0 wins, and across types
    # the first type reaching the highest score wins (strict >, so ties keep the
    # earlier type). No interaction at all falls back to the location.
    best_type: int = -1
    best_score: float = 0.0

    n_types, n_entries = score_indices.shape
    for type_index in range(n_types):
        type_max: float = 0.0
        for entry_index in range(n_entries):
            flat_index = score_indices[type_index, entry_index]
            if flat_index < 0:
                break
            score = profile[flat_index]
            # NaN fails both comparisons, matching the score != score guard.
            type_max = max(type_max, score)
        if type_max > best_score:
            best_score = type_max
            best_type = type_index

    if best_type < 0:
        return location_code
    return best_type


def simple_class_code_from_encoded(encoded_profile: np.ndarray, location_code: int) -> int:
    return _simple_class_code(encoded_profile, _CLASS_SCORE_INDICES, location_code)


@njit(cache=True)
def _simple_class_codes_for_columns(
        profiles: np.ndarray, score_indices: np.ndarray, sc_rsa: np.ndarray,
        surface_threshold: float, buried_threshold: float, out: np.ndarray) -> np.ndarray:
    # The same two calls the aggregation kernel makes per mapping, over a whole
    # chain at once. Both are called rather than reimplemented, so a code stored
    # by this path cannot drift from the one recomputed by that one.
    for row in range(profiles.shape[0]):
        out[row] = _simple_class_code(
            profiles[row], score_indices,
            locate_code(sc_rsa[row], surface_threshold, buried_threshold))
    return out


def simple_class_codes_from_columns(
        profiles: np.ndarray, sc_rsa: np.ndarray,
        surface_threshold: float, buried_threshold: float) -> np.ndarray[np.int8]:
    """One simple class code per row, from the profile block and the side chain RSA column.

    The codes span -4..7 (see SIMPLE_CLASS_CODE_OF_STRING), so int8 holds them with
    room to spare - the point of storing them at all is that they are small next to
    the 195 float32 encoded row they are derived from.
    """
    out: np.ndarray[np.int8] = np.empty(profiles.shape[0], dtype=np.int8)
    if profiles.shape[0] == 0:
        return out
    return _simple_class_codes_for_columns(
        profiles, _CLASS_SCORE_INDICES, sc_rsa, surface_threshold, buried_threshold, out)


# Public handles for callers that run the classification inside a kernel of their
# own. The aggregation loop does exactly that: calling locate_code and
# _simple_class_code from python cost it two dispatches per mapping, billions of
# times per run, so it inlines both into its own njit function and needs the score
# index matrix as an argument.
CLASS_SCORE_INDICES: np.ndarray = _CLASS_SCORE_INDICES
simple_class_code_kernel = _simple_class_code


def compute_class_from_encoded(encoded_profile) -> tuple[str | None, str | None]:
    if encoded_profile is None:
        return None, None

    try:
        values = encoded_profile.tolist()
    except AttributeError:
        values = encoded_profile

    interactions = []
    for interaction_type, entries in _CLASS_SCORE_LAYOUT:
        max_score = 0.
        max_score_interaction = None
        for index, chain_type, bond_type in entries:
            score = values[index]
            if score is None or score != score:  # NaN
                continue
            if score > max_score:
                max_score_interaction = (chain_type, bond_type)
                max_score = score
        if max_score_interaction is not None:
            interactions.append((max_score, max_score_interaction[0], max_score_interaction[1], interaction_type))

    s = len(interactions)
    if s == 0:
        return 'No interaction', 'No interaction'

    if s == 1:
        score, chain_type, bond_type, interaction_type = interactions[0]
        return f'{chain_type} {bond_type} {interaction_type}', f'{interaction_type} interaction'

    i_str_parts = []
    max_score = 0.
    simple_class = None
    for score, chain_type, bond_type, interaction_type in interactions:
        i_str_parts.append(f'{chain_type} {bond_type} {interaction_type}')
        if score > max_score:
            max_score = score
            simple_class = f'{interaction_type} interaction'

    return 'Multi interaction: {}'.format(' and '.join(i_str_parts)), simple_class


def classify_encoded(encoded_profile, location) -> tuple[str | None, str | None]:
    if encoded_profile is None:
        return None, None

    raw_rin_class, raw_rin_simple_class = compute_class_from_encoded(encoded_profile)

    if raw_rin_class == 'No interaction':
        rin_class = location
    else:
        rin_class = raw_rin_class
    if raw_rin_simple_class == 'No interaction':
        rin_simple_class = location
    else:
        rin_simple_class = raw_rin_simple_class
    return rin_class, rin_simple_class


class Centrality_scores(Slotted_obj):
    __slots__ = ['AbsoluteCentrality', 'LengthNormalizedCentrality', 'MinMaxNormalizedCentrality', 'AbsoluteCentralityWithNegative',
                 'LengthNormalizedCentralityWithNegative', 'MinMaxNormalizedCentralityWithNegative', 'AbsoluteComplexCentrality',
                 'LengthNormalizedComplexCentrality', 'MinMaxNormalizedComplexCentrality', 'AbsoluteComplexCentralityWithNegative',
                 'LengthNormalizedComplexCentralityWithNegative', 'MinMaxNormalizedComplexCentralityWithNegative']
    
    feature_names = ['AbsoluteCentrality', 'LengthNormalizedCentrality', 'MinMaxNormalizedCentrality', 'AbsoluteCentralityWithNegative',
                 'LengthNormalizedCentralityWithNegative', 'MinMaxNormalizedCentralityWithNegative', 'AbsoluteComplexCentrality',
                 'LengthNormalizedComplexCentrality', 'MinMaxNormalizedComplexCentrality', 'AbsoluteComplexCentralityWithNegative',
                 'LengthNormalizedComplexCentralityWithNegative', 'MinMaxNormalizedComplexCentralityWithNegative']

    def __init__(self):
        pass

    def to_list(self):
        data_list = []
        for slot in self.__slots__:
            try:
                d = self.__getattribute__(slot)
            except AttributeError:
                d = float('nan')
            data_list.append(d)
        return data_list
    
    def set_data_from_list(self, value_list):
        for index, slot in enumerate(self.__slots__):
            self.__setattr__(slot, value_list[index])

class Inflow_values(Centrality_scores):
    __slots__ = [
        'InflowValue', 'OutflowValue',
        'WeightedInflowValue', 'WeightedOutflowValue',
        'InflowValueFirstPass', 'OutflowValueFirstPass',
        'WeightedInflowValueFirstPass', 'WeightedOutflowValueFirstPass'
        ]
    
    feature_names = [x for x in __slots__]
    

# called by database
# called by sdsc
def calculateAverageProfile(profiles):
    #print(f'call of calculateAverageProfiels: {len(profiles)=}')
    if len(profiles) == 0:
        return None
    average_profile = Interaction_profile()
    for chain_type in chain_types:
        for bond_type in bond_types:
            for interaction_type in interaction_types:
                weight_sum = 0.
                degree_sum = 0.
                score_sum = 0.
                for weight, profile in profiles:
                    if profile is None:
                        #print('Profile is None')
                        continue
                    degree = profile.getDegree(chain_type, bond_type, interaction_type)
                    score = profile.getScore(chain_type, bond_type, interaction_type)
                    if np.isnan(degree) or np.isnan(score):
                        continue
                    weight_sum += weight
                    try:
                        degree_sum += float(degree) * weight
                    except TypeError:
                        ...
                    score_sum += score * weight
                if weight_sum == 0.:
                    continue
                average_degree = degree_sum / weight_sum
                average_score = score_sum / weight_sum
                average_profile.setDegree(chain_type, bond_type, interaction_type, average_degree)
                average_profile.setScore(chain_type, bond_type, interaction_type, average_score)
    return average_profile

def disect_rinerator_node_id(node_id: str):
    chain, res_nr, insertioncode, res_name = node_id.split(':')
    return chain, res_nr, insertioncode, res_name

def _contig_index(res_contig_map: dict) -> dict[str, dict]:
    """Per chain, a plain {residue number: contig index} dict.

    The parse loop asked res_contig_map[chain].contains(x) and then
    .get_item(x)[0], which is two full Residue_Map lookups - each one an isinstance
    test, an int() attempt that raises for insertion-code keys, a SparseArray probe
    and a dict fallback - and it did that four or more times per line. Millions of
    lines later that was the single biggest cost in the parse. The map is built once
    here and every lookup afterwards is a dict.get.

    The keys are exactly what Residue_Map holds: ints for plain residue numbers,
    strings for the ones carrying an insertion code, which is the same pair of forms
    the parsed res_nr can take.
    """
    index = {}
    for chain, residue_map in res_contig_map.items():
        chain_index = {}
        for key in residue_map.get_keys():
            value = residue_map.get_item(key)
            if value is not None:
                chain_index[key] = value[0]
        index[chain] = chain_index
    return index


def parse_rinerator_interaction_score_file(
        interaction_score_file: str,
        struct_ligands: set[tuple[str, str | int]],
        metals: set[tuple[str, str | int]],
        ions: set[tuple[str, str | int]],
        chain_type_map: dict[str, str],
        homomer_set: set[str],
        res_contig_map: dict[str, Residue_Map[tuple[int, str]]],
        page: str | None = None,
        ) -> tuple[list[tuple[str, int, str, int, str, str , float]], dict[str, Residue_Map[Interaction_profile]]]:

    if page is None:
        f = gzip.open(interaction_score_file, 'rt')
        lines = f.readlines()
        f.close()
    else:
        lines = page.splitlines(keepends=True)

    interactions: list[tuple[str, int | str, str, int | str, str, str , float]] = []

    profiles: dict[str, list[Interaction_profile]] = {}
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}

    contig_index = _contig_index(res_contig_map)
    profile_lengths = {chain: len(index) for chain, index in contig_index.items()}
    boring_ligands = ligands.BORING_LIGANDS
    three_to_one = residues.THREE_TO_ONE

    for line in lines[1:]:
        # A:39:_:ASP (cnt:mc_sc) A:41:_:LEU    2.220600
        parts = line.split('\t')
        if len(parts) == 1:
            print(f'{interaction_score_file=}')
            break
        score = float(parts[1])

        #IGNORE ALL OVERLAPS
        if score < 0.:
            continue

        # One split of the edge instead of a split plus two node splits: the node
        # ids are colon separated and the edge is space separated, so the eight
        # fields come out in a single pass.
        res_a, interaction_type, res_b = parts[0].split()
        chain_a, res_nr_a, insertioncode_a, res_name_a = res_a.split(':')
        chain_b, res_nr_b, insertioncode_b, res_name_b = res_b.split(':')
        interaction_base_type: Literal['cnt', 'ovl', 'hbond', 'combi']
        interaction_base_type, interaction_sub_type = interaction_type[1:-1].split(':')
        interaction_type_a: Literal['sc', 'mc', 'ligand', 'all']
        interaction_type_b: Literal['sc', 'mc', 'ligand', 'all']
        interaction_type_a, interaction_type_b = interaction_sub_type.split('_')

        # An insertion code of '_' means there is none, which is almost every
        # residue, so that case skips the format and the failed int() entirely.
        if insertioncode_a == '_':
            try:
                res_nr_a = int(res_nr_a)
            except ValueError:
                pass
        else:
            res_nr_a = f'{res_nr_a}{insertioncode_a}'
            try:
                res_nr_a = int(res_nr_a)
            except ValueError:
                pass
        if insertioncode_b == '_':
            try:
                res_nr_b = int(res_nr_b)
            except ValueError:
                pass
        else:
            res_nr_b = f'{res_nr_b}{insertioncode_b}'
            try:
                res_nr_b = int(res_nr_b)
            except ValueError:
                pass

        if res_name_a in boring_ligands and res_name_a not in three_to_one:
            continue
        if res_name_b in boring_ligands and res_name_b not in three_to_one:
            continue

        index_a = contig_index.get(chain_a)
        index_b = contig_index.get(chain_b)

        if index_a is not None:
            res_conti_a = index_a.get(res_nr_a)
            if res_conti_a is not None:

                if index_b is not None:
                    res_conti_b = index_b.get(res_nr_b)
                    if res_conti_b is not None:

                        if chain_a != chain_b: #Only PPI
                            interactions.append((chain_a, res_conti_a, chain_b, res_conti_b, interaction_base_type, interaction_type_a, score))
                            interactions.append((chain_b, res_conti_b, chain_a, res_conti_a, interaction_base_type, interaction_type_b, score))

                chain_profiles = profiles.get(chain_a)
                if chain_profiles is None:
                    chain_profiles = [None] * profile_lengths[chain_a]
                    profiles[chain_a] = chain_profiles
                profile = chain_profiles[res_conti_a-1]
                if profile is None:
                    profile = Interaction_profile()
                    chain_profiles[res_conti_a-1] = profile

                add_interaction_to_profile(chain_a, res_nr_a, chain_b, res_nr_b, struct_ligands, metals, ions,
                                                    chain_type_map, homomer_set,
                                                    index_a, profile,
                                                    interaction_type_a, interaction_base_type, score, ligand_profiles,
                                                    metal_profiles,
                                                    ion_profiles,
                                                    chain_chain_profiles)

        if index_b is not None:
            res_conti_b = index_b.get(res_nr_b)
            if res_conti_b is not None:

                chain_profiles = profiles.get(chain_b)
                if chain_profiles is None:
                    chain_profiles = [None] * profile_lengths[chain_b]
                    profiles[chain_b] = chain_profiles
                profile = chain_profiles[res_conti_b-1]
                if profile is None:
                    profile = Interaction_profile()
                    chain_profiles[res_conti_b-1] = profile

                add_interaction_to_profile(chain_b, res_nr_b, chain_a, res_nr_a, struct_ligands, metals, ions,
                                                    chain_type_map, homomer_set,
                                                    index_b, profile,
                                                    interaction_type_b, interaction_base_type, score,
                                                    ligand_profiles,
                                                    metal_profiles,
                                                    ion_profiles,
                                                    chain_chain_profiles)
    return interactions, profiles, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles


def float_or_none(string : str) -> float | None:
    try:
        return float(string)
    except:
        return None
    
def parse_inflow_file(
        inflow_file: str, 
        res_contig_map: dict[str, Residue_Map[tuple[int, str]]],
        config: Config,
        page=None,
        ):
    if page is None:
        f = open(inflow_file, 'r')
        page = f.read()
        f.close()
    else:
        page: str
    try:
        lines = page.split('\n')
    except TypeError:
        page = page.decode()
        lines = page.split('\n')

    #Node,InflowValue,OutflowValue,Degree,WeightSum,WeightedInflowValue,WeightedOutflowValue,InflowValueFirstPass,OutflowValueFirstPass,WeightedInflowValueFirstPass,WeightedOutflowValueFirstPass

    # -1 , 0         , 1          , 2    , 3       , 4                 , 5                  , 6                  , 7                   , 8                          , 9                                       

    inflow_map: dict[str, np.ndarray] = {}

    for line in lines[1:]:
        if line == '':
            continue
        words = line.split(',')
        node_id = words[0]
        chain, res_nr, insertioncode, res_name = disect_rinerator_node_id(node_id)

        if chain not in res_contig_map:
            continue

        insertioncode = insertioncode.replace('_', '')

        res_id = f'{res_nr}{insertioncode}'

        if chain not in inflow_map:
            inflow_map[chain] = np.zeros((len(res_contig_map[chain].get_keys()), 8), dtype=float) 

        if not res_contig_map[chain].contains(res_id):
            #config.logger.info(f'Residue not in res_contig_map: {inflow_file=} {chain=} {res_id=}')
            continue

        values = np.array(words[1:], dtype=float)
        values = values[[0, 1, 4, 5, 6, 7, 8, 9]]

        res_contig: int = res_contig_map[chain].get_item(res_id)[0]

        inflow_map[chain][res_contig-1] = values
    return inflow_map

def getCentMap(centrality_file: str, res_contig_map: dict[str, Residue_Map[tuple[int, str]]], page: str | None = None):

    if page is None:
        f = gzip.open(centrality_file, 'rt')
        lines = f.read().split('\n')
        f.close()
    else:
        lines = page.split('\n')

    centrality_map = {}

    for line in lines[1:]:
        # Residue        AbsoluteCentrality      LengthNormalizedCentrality      MinMaxNormalizedCentrality      AbsoluteCentralityWithNegative  LengthNormalizedCentralityWithNegative  MinMaxNormalizedCentralityWithNegative  AbsoluteComplexCentrality       LengthNormalizedComplexCentrality       MinMaxNormalizedComplexCentrality       AbsoluteComplexCentralityWithNegative   LengthNormalizedComplexCentralityWithNegative   MinMaxNormalizedComplexCentralityWithNegative
        # A:20:_:TYR      1927.0  0.24084489438820147     1.0     1947.0  0.24334458192725908     1.0     2412.0  0.07389705882352941     1.2516865594187856      2260.0  0.06924019607843138     1.1607601438109914

        # A:9:_:TYR    1036.0    0.150079675503    0.829463570857
        #print(f'{centrality_file} {line}')

        if line == '':
            continue
        if line[0] == '#':
            continue
        if line.count('\t') == 0:
            print(centrality_file, line)
            continue
        try:
            words = line.split('\t')
            res = words[0]

            cent_scores = Centrality_scores()

            for pos, attr_name in enumerate(cent_scores.__slots__, start = 1):
                cent_scores.__setattr__(attr_name, float_or_none(words[pos]))

            [chain, res_nr, insertioncode, res_name] = res.split(':')

            if chain not in res_contig_map:
                continue

            res_nr = "{}{}".format(res_nr, insertioncode.replace('_', ''))
            try:
                res_nr = int(res_nr)
            except:
                pass

            if chain not in centrality_map:
                centrality_map[chain] = [None] * len(res_contig_map[chain].get_keys())

            if not res_contig_map[chain].contains(res_nr):
                continue

            centrality_map[chain][res_contig_map[chain].get_item(res_nr)[0]-1] = cent_scores
        except:
            # some of the rins are not updated, this leads to an error here,-> simply update all rins
            print('Error in getCentMap: ', centrality_file)
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            print('\n'.join([str(e), str(f), str(g)]))
            continue

    return centrality_map

def add_to_mol_type_profile(chain, res, score, mol_type_profile):
    if chain not in mol_type_profile:
        mol_type_profile[chain] = Residue_Map()
    if not mol_type_profile[chain].contains(res):
        mol_type_profile[chain].add_item(res, (1, score))
    else:
        deg, total_score = mol_type_profile[chain].get_item(res)
        deg += 1
        total_score += score
        mol_type_profile[chain].add_item(res, (deg, total_score))

def add_interaction_to_profile(
        chain: str,
        res: int | str,
        chain_b: str,
        res_b: int | str,
        ligands: set[tuple[str, int | str]],
        metals: set[tuple[str, int | str]],
        ions: set[tuple[str, int | str]],
        chain_type_map: dict[str, str],
        homomer_set: set[str],
        chain_contig_index: dict,
        profile: Interaction_profile,
        chaintype: str,
        bondtype: str,
        score: float,
        ligand_profiles,
        metal_profiles,
        ion_profiles,
        chain_chain_profiles
        ) -> None:
    
    if bondtype == 'combi':
        if chaintype != 'all':
            return
        if (chain_b, res_b) in ligands:
            add_to_mol_type_profile(chain, res, score, ligand_profiles)
        elif (chain_b, res_b) in metals:
            add_to_mol_type_profile(chain, res, score, metal_profiles)
        elif (chain_b, res_b) in ions:
            add_to_mol_type_profile(chain, res, score, ion_profiles)
        elif chain_b != chain:
            add_to_mol_type_profile(chain, chain_b, score, chain_chain_profiles)
        return
    
    if (chain_b, res_b) in ligands:
        interaction_type = 'ligand'
        
    elif (chain_b, res_b) in metals:
        interaction_type = 'metal'
    elif (chain_b, res_b) in ions:
        interaction_type = 'ion'
    elif chain != chain_b:
        if chain_b not in chain_type_map:
            return
        if chain_b in homomer_set:
            interaction_type = 'Homomer'
        else:
            interaction_type = chain_type_map[chain_b]
    else:
        # Intra-chain, so both residues are in this chain's own contig index. Four
        # Residue_Map lookups became two dict gets.
        res_conti = chain_contig_index.get(res)
        if res_conti is None:
            return
        res_conti_b = chain_contig_index.get(res_b)
        if res_conti_b is None:
            return
        res_dist = abs(res_conti - res_conti_b)
        if res_dist < 2:
            interaction_type = 'neighbor'
        elif res_dist < 6:
            interaction_type = 'short'
        else:
            interaction_type = 'long'

    profile.addEdge(chaintype, bondtype, interaction_type, score, chain_b, res_b)


def calculateIAPProfiles(
        interaction_map: list,
        chains: list[str],
        ligands: set[tuple[str, str]],
        metals: set[tuple[str, str]],
        ions: set[tuple[str, str]]):
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}

    for chain, res in ligands:
        deg = 0
        total_score = 0.0
        if chain not in interaction_map:
            # This can happen if the chain consists only of ligands with at least one non-boring ligand
            # and this non-boring ligand only has bonds with boring ligands
            continue
        if interaction_map[chain].contains(res):
            for chain_b in interaction_map[chain].get_item(res):
                for res_b in interaction_map[chain].get_item(res)[chain_b].get_keys():
                    interactions = interaction_map[chain].get_item(res)[chain_b].get_item(res_b)
                    for bondtype, chaintype, score in interactions:
                        if bondtype != 'combi':
                            continue
                        if chaintype != 'all':
                            continue
                        total_score += score
                        deg += 1
        ligand_profiles[(chain, res)] = [deg, total_score]

    #print(len(ligand_profiles))

    for chain, res in metals:
        deg = 0
        total_score = 0.0
        if interaction_map[chain].contains(res):
            for chain_b in interaction_map[chain].get_item(res):
                for res_b in interaction_map[chain].get_item(res)[chain_b].get_keys():
                    interactions = interaction_map[chain].get_item(res)[chain_b].get_item(res_b)
                    for bondtype, chaintype, score in interactions:
                        if bondtype != 'combi':
                            continue
                        if chaintype != 'all':
                            continue
                        total_score += score
                        deg += 1
        metal_profiles[(chain, res)] = [deg, total_score]

    #print(len(metal_profiles))

    for chain, res in ions:
        deg = 0
        total_score = 0.0
        if interaction_map[chain].contains(res):
            for chain_b in interaction_map[chain].get_item(res):
                for res_b in interaction_map[chain].get_item(res)[chain_b].get_keys():
                    interactions = interaction_map[chain].get_item(res)[chain_b].get_item(res_b)
                    for bondtype, chaintype, score in interactions:
                        if bondtype != 'combi':
                            continue
                        if chaintype != 'all':
                            continue
                        total_score += score
                        deg += 1
        ion_profiles[(chain, res)] = [deg, total_score]

    #print(len(ion_profiles))

    for chain in chains:
        #print(chain)
        if chain  not in interaction_map:
            continue
        for res in interaction_map[chain].get_keys():
            if (chain, res) in ligands:
                continue
            if (chain, res) in metals:
                continue
            if (chain, res) in ions:
                continue
            for chain_b in interaction_map[chain].get_item(res):
                if chain == chain_b:
                    continue
                for res_b in interaction_map[chain].get_item(res)[chain_b].get_keys():
                    if (chain_b, res_b) in ligands:
                        continue
                    if (chain_b, res_b) in metals:
                        continue
                    if (chain_b, res_b) in ions:
                        continue
                    if (chain, chain_b) not in chain_chain_profiles:
                        chain_chain_profiles[(chain, chain_b)] = [0, 0.0]

                    interactions = interaction_map[chain].get_item(res)[chain_b].get_item(res_b)
                    for bondtype, chaintype, score in interactions:
                        if bondtype != 'combi':
                            continue
                        if chaintype != 'all':
                            continue

                        chain_chain_profiles[(chain, chain_b)][0] += 1
                        chain_chain_profiles[(chain, chain_b)][1] += score

    return ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles


# called by templateFiltering

def lookup(
        structure_id: str,
        page: str,
        config: Config,
        ligands: set[tuple[str, int | str]], 
        metals: set[tuple[str, int | str]],
        ions:set[tuple[str, int | str]],
        res_contig_map: dict[str, Residue_Map[tuple[int, str]]],
        base_path: str,
        chain_type_map: dict[str, str],
        homomer_set: set[str],
        model_path : str | None = None,
        keep_tmp_files: bool = True
        ) -> tuple[
            dict[str, Residue_Map[Interaction_profile]],
            dict,
            dict,
            dict,
            dict,
            dict[str, Residue_Map[dict[str, Residue_Map[list[tuple[str, str, float]]]]]]
        ]:

    runtimes = []
    t0 = time.time()
    if config.verbosity >= 4:
        config.logger.info(f'Call of rin.lookup: {structure_id=} {base_path=} {model_path=} {keep_tmp_files=} {config.custom_db_path=}')

    folder_path, path_stem, structure_id = get_entry_path(structure_id, config.rin_db_path, custom_db_path= config.custom_db_path, model_path=model_path, path_to_model_db=config.path_to_model_db, temp_folder=config.tmp_folder)

    if config.verbosity >= 4:
        config.logger.info(f'Paths in rin.lookup: {structure_id=} {folder_path=} {path_stem=}')

    afdb_rin_page = None
    afdb_cent_page = None
    afdb_iof_page = None
    if is_alphafold_model(structure_id):
        if config.afdb is not None:
            rin_and_cent = get_afdb_rin_and_cent(structure_id, config)
            if rin_and_cent is not None:
                afdb_rin_page, afdb_cent_page, afdb_iof_page = rin_and_cent
                if afdb_rin_page is not None:
                    afdb_rin_page = afdb_rin_page.decode('ascii')
                    afdb_cent_page = afdb_cent_page.decode('ascii')
                    try:
                        afdb_iof_page = afdb_iof_page.decode('ascii')
                    except AttributeError as err:
                        raise AttributeError(f'{err=} {structure_id=}')
            else:
                config.errorlog.add_warning(f'RIN and centralities missing from AFDB for {structure_id=}')

    interaction_score_file = f"{path_stem}_intsc.ea.gz"

    if not os.path.isfile(interaction_score_file) and afdb_rin_page is None:
        if config.verbosity >= 3:
            config.logger.info(f'Did not find RIN: {interaction_score_file}, folder path: {folder_path=} {config.rin_db_path=} {config.custom_db_path=} {model_path=} {config.tmp_folder=}' )

        if not os.path.isdir(folder_path):
            try:
                os.makedirs(folder_path)
            except OSError:
                return f'RIN not found and calculation of RIN failed for: {structure_id=}, {model_path=} {interaction_score_file=}'

        calcRIN(config, page.encode(), folder_path, path_stem, structure_id, config.rinerator_path, True, structure_path=model_path)

        interaction_score_file = f"{path_stem}_intsc.ea.gz"
        
        if not os.path.isfile(interaction_score_file):
            message = f'RIN not found and calculation of RIN failed for: {structure_id=}, {model_path=} {interaction_score_file=}'
            if config.verbosity >= 3:
                config.logger.info(message)
            return message

    network_file = f"{path_stem}.sif.gz"
    interaction_count_file = f"{path_stem}_nrint.ea.gz"
    residue_file = f"{path_stem}_res.txt.gz"
    centrality_file = f"{path_stem}_btw_cent.txt.gz"
    inflow_file = f'{path_stem}_in_out_flow.csv'

    t1 = time.time()
    runtimes.append(t1-t0)
    if config.verbosity >= 4:
        config.logger.info(f'Call of parse_rinerator_interaction_score_file: {structure_id} {interaction_score_file} {t1-t0} {res_contig_map=} {chain_type_map=} {ligands=}')

    profiles: dict[str, Residue_Map[Interaction_profile]]
    (
        interaction_list,
        profiles,
        ligand_profiles,
        metal_profiles,
        ion_profiles,
        chain_chain_profiles
    ) = parse_rinerator_interaction_score_file(
        interaction_score_file, ligands, metals, ions,
        chain_type_map, homomer_set, res_contig_map, page = afdb_rin_page)

    t2 = time.time()
    runtimes.append(t2-t1)
    if config.verbosity >= 4:
        config.logger.info(f'Call of getCentMap: {structure_id}  {centrality_file} {t2-t1}')

    centrality_map: dict[str, Residue_Map[Centrality_scores]] = {}

    if os.path.isfile(centrality_file):
        centrality_map = getCentMap(centrality_file, res_contig_map, page=afdb_cent_page)
    else:
        centrality_map = {}

    chains = list(res_contig_map.keys())

    if afdb_iof_page is None:
        if not os.path.isfile(inflow_file):
            call_in_out_flow(f"{path_stem}_intsc.ea.gz", inflow_file)

        inflow_map = parse_inflow_file(inflow_file, res_contig_map, config)
    else:
        inflow_map = parse_inflow_file(inflow_file, res_contig_map, config, page=afdb_iof_page)

    t3 = time.time()
    runtimes.append(t3-t2)
    if config.verbosity >= 4:
        config.logger.info(f'After getCentMap: {structure_id} {centrality_file} {t3-t2}')

    if len(chains) == 0:
        return f'Illegal input, no chains chosen in RIN lookup of {structure_id}'

    t4 = time.time()
    runtimes.append(t4-t3)
    if config.verbosity >= 4:
        config.logger.info(f'Call of calculateIAPProfiles: {structure_id} {t4-t3}')

    if len(profiles) == 0:
        return f'Profiles_map is empty, in RIN lookup of {structure_id}'
    
    #ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles = calculateIAPProfiles(interaction_list, chains, ligands, metals, ions)

    if not keep_tmp_files:
        os.remove(interaction_score_file)
        os.remove(network_file)
        os.remove(interaction_count_file)
        os.remove(residue_file)
        os.remove(centrality_file)
    else:
        if model_path is not None:
            model_directory = os.path.dirname(model_path)
            network_filename = os.path.basename(network_file)
            target_path = f'{model_directory}/{network_filename}'
            if network_file != target_path:
                shutil.copyfile(network_file, target_path)
            target_score_file = f'{model_directory}/{os.path.basename(interaction_score_file)}'
            if interaction_score_file != target_score_file:
                shutil.copyfile(interaction_score_file, target_score_file)


    t5 = time.time()
    runtimes.append(t5-t4)
    if config.verbosity >= 4:
        config.logger.info(f'At end of rin.lookup: {structure_id} {t5-t4} {t5-t0}')

    return profiles, centrality_map, inflow_map, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, interaction_list, runtimes
