import gzip
import os
import sys
import traceback
import time

from structman.base_utils.config_class import Config
from structman.lib.createRINdb import calcRIN, get_entry_path
from structman.lib.sdsc.consts import residues, ligands
from structman.lib.sdsc.sdsc_utils import Slotted_obj
from structman.lib.sdsc import residue as residue_package


class Interaction_type(Slotted_obj):
    def addEdge(self, interaction_type: str, score: float) -> None:
        try:
            try:
                interaction_data = self.__getattribute__(interaction_type)
            except AttributeError:
                interaction_data = [0, 0.0]
            try:
                deg, old_score = interaction_data
            except TypeError:
                deg = 0
                old_score = 0.0
            self.__setattr__(interaction_type, [deg + 1, old_score + score])
        except:
            raise KeyError(f'Unknown interactiontype: {interaction_type}')

    def __init__(self):
        ...
        #for slot in self.__slots__:
        #    self.__setattr__(slot, None)

    def getScore(self, interaction_type: str) -> float:
        try:
            return self.__getattribute__(interaction_type)[1]
        except (AttributeError, TypeError):
            self.__setattr__(interaction_type, [0, 0.0])
            return self.__getattribute__(interaction_type)[1]

    def getDegree(self, interaction_type: str) -> int:
        try:
            return self.__getattribute__(interaction_type)[0]
        except (AttributeError, TypeError):
            self.__setattr__(interaction_type, [0, 0.0])
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
            self.__setattr__(interaction_type, [degree, 0.0])

class Ligand_types(Interaction_type):
    __slots__ = ['ligand', 'ion', 'metal']

ligand_types = set(['ligand', 'ion', 'metal'])

class Intrachain_types(Interaction_type):
    __slots__ = ['neighbor', 'short', 'long']

intrachain_types = set(['neighbor', 'short', 'long'])

class Interchain_types(Interaction_type):
    __slots__ = ['Protein', 'DNA', 'RNA', 'Peptide']

interchain_types = set(['Protein', 'DNA', 'RNA', 'Peptide'])

interaction_types = set(['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide'])

non_intra_interaction_types = set(['ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide'])

chain_type_interactions = set(['Protein', 'DNA', 'RNA', 'Peptide'])

small_molecule_interactions = set(['ligand', 'ion', 'metal'])


class Interaction_type_profile(Slotted_obj):
    __slots__ = ['sw_ligand', 'interchain', 'intrachain']

    def __init__(self):
        self.sw_ligand = Ligand_types()
        self.interchain = Interchain_types()
        self.intrachain = Intrachain_types()


    def addEdge(self, interaction_type, score):
        if interaction_type in interchain_types:
            return self.interchain.addEdge(interaction_type, score)
        elif interaction_type in intrachain_types:
            return self.intrachain.addEdge(interaction_type, score)
        elif interaction_type in ligand_types:
            return self.sw_ligand.addEdge(interaction_type, score)
        else:
            raise TypeError(f'Unknown interactiontype: {interaction_type}')

    def getScore(self, interaction_type):
        if interaction_type in interchain_types:
            return self.interchain.getScore(interaction_type)
        elif interaction_type in intrachain_types:
            return self.intrachain.getScore(interaction_type)
        elif interaction_type in ligand_types:
            return self.sw_ligand.getScore(interaction_type)
        else:
            return None

    def getDegree(self, interaction_type):
        if interaction_type in interchain_types:
            return self.interchain.getDegree(interaction_type)
        elif interaction_type in intrachain_types:
            return self.intrachain.getDegree(interaction_type)
        elif interaction_type in ligand_types:
            return self.sw_ligand.getDegree(interaction_type)
        else:
            return None

    def setScore(self, interaction_type, score):
        if interaction_type in interchain_types:
            return self.interchain.setScore(interaction_type, score)
        elif interaction_type in intrachain_types:
            return self.intrachain.setScore(interaction_type, score)
        elif interaction_type in ligand_types:
            return self.sw_ligand.setScore(interaction_type, score)
        else:
            return 'Unknown interactiontype: %s' % interaction_type

    def setDegree(self, interaction_type, degree):
        if interaction_type in interchain_types:
            return self.interchain.setDegree(interaction_type, degree)
        elif interaction_type in intrachain_types:
            return self.intrachain.setDegree(interaction_type, degree)
        elif interaction_type in ligand_types:
            return self.sw_ligand.setDegree(interaction_type, degree)
        else:
            return 'Unknown interactiontype: %s' % interaction_type


bond_types = set(['cnt', 'hbond', 'ovl'])


class Chain_types(Slotted_obj):
    __slots__ = ['contact', 'h_bond', 'overlap']

    def __init__(self):
        self.contact = Interaction_type_profile()
        self.h_bond = Interaction_type_profile()
        self.overlap = Interaction_type_profile()

    def addEdge(self, bondtype, interaction_type, score):
        if bondtype == 'cnt':
            return self.contact.addEdge(interaction_type, score)
        elif bondtype == 'hbond':
            return self.h_bond.addEdge(interaction_type, score)
        elif bondtype == 'ovl':
            return self.overlap.addEdge(interaction_type, score)
        else:
            return 'Unknown bondtype: %s' % bondtype

    def getScore(self, bondtype, interaction_type):
        if bondtype == 'cnt':
            return self.contact.getScore(interaction_type)
        elif bondtype == 'hbond':
            return self.h_bond.getScore(interaction_type)
        elif bondtype == 'ovl':
            return self.overlap.getScore(interaction_type)
        else:
            return 'Unknown bondtype: %s' % bondtype

    def getDegree(self, bondtype, interaction_type):

        if bondtype == 'cnt':
            return self.contact.getDegree(interaction_type)
        elif bondtype == 'hbond':
            return self.h_bond.getDegree(interaction_type)
        elif bondtype == 'ovl':
            return self.overlap.getDegree(interaction_type)
        else:
            return 'Unknown bondtype: %s' % bondtype

    def setScore(self, bondtype, interaction_type, score):

        if bondtype == 'cnt':
            return self.contact.setScore(interaction_type, score)
        elif bondtype == 'hbond':
            return self.h_bond.setScore(interaction_type, score)
        elif bondtype == 'ovl':
            return self.overlap.setScore(interaction_type, score)
        else:
            return 'Unknown bondtype: %s' % bondtype

    def setDegree(self, bondtype, interaction_type, degree):

        if bondtype == 'cnt':
            return self.contact.setDegree(interaction_type, degree)
        elif bondtype == 'hbond':
            return self.h_bond.setDegree(interaction_type, degree)
        elif bondtype == 'ovl':
            return self.overlap.setDegree(interaction_type, degree)
        else:
            return 'Unknown bondtype: %s' % bondtype

    def getCombiScore(self, interaction_type):
        return self.contact.getScore(interaction_type) + self.h_bond.getScore(interaction_type) + self.overlap.getScore(interaction_type)

    def getCombiDegree(self, interaction_type):
        return self.contact.getDegree(interaction_type) + self.h_bond.getDegree(interaction_type) + self.overlap.getDegree(interaction_type)


chain_types = set(['mc', 'sc'])


class Interaction_profile(Slotted_obj):
    __slots__ = [
        'mainchain', 'sidechain', '_class',
        '_simple_class', 'interacting_chains', 'interacting_ligands'
        ]

    slot_mask = [True, True, False, False, False, False]

    def __init__(self):
        #self.mainchain = Chain_types()
        #self.sidechain = Chain_types()
        self.interacting_chains = set()
        self.interacting_ligands = set()

    def addEdge(self, chaintype, bondtype, interaction_type, score, chain_id, res_id):
        if bondtype != 'ovl':
            if interaction_type in small_molecule_interactions:
                self.interacting_ligands.add('%s:%s' % (chain_id, res_id))
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
            return 'Unknown chaintype in addEdge: %s' % chaintype

    def getScore(self, chain_type, bond_type, interaction_type):
        if chain_type == 'mc':
            try:
                return self.mainchain.getScore(bond_type, interaction_type)
            except AttributeError:
                self.mainchain = Chain_types()
                return self.mainchain.getScore(bond_type, interaction_type)
            
        elif chain_type == 'sc':
            try:
                return self.sidechain.getScore(bond_type, interaction_type)
            except AttributeError:
                self.sidechain = Chain_types()
                return self.sidechain.getScore(bond_type, interaction_type)
        else:
            return 'Unknown chaintype in getScore: %s' % chain_type

    def getDegree(self, chain_type, bond_type, interaction_type):
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
            return 'Unknown chaintype in getDegree: %s' % chain_type

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
            return 'Unknown chaintype in setScore: %s' % chain_type

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
            return 'Unknown chaintype in setDegree: %s' % chain_type

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

    def getTotalInterfaceInteractionScore(self):
        score = 0.
        for chain_type in chain_types:
            for i_type in non_intra_interaction_types:
                s = self.getChainSpecificCombiScore(chain_type, i_type)
                if s is not None and not isinstance(s, str):
                    score += s
        return score

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

    def computeClass(self):
        #print('Call of computeClass')
        interactions = []
        b_types = ['cnt', 'hbond']  # ignore overlaps
        i_types = ['ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']
        for interaction_type in i_types:
            max_score_interaction = None
            max_score = 0.
            for chain_type in chain_types:
                for bond_type in b_types:
                    score = self.getScore(chain_type, bond_type, interaction_type)
                    if score > max_score:
                        max_score_interaction = (chain_type, bond_type)
                        max_score = score
            if max_score > 0.:
                (chain_type, bond_type) = max_score_interaction
                interactions.append((max_score, chain_type, bond_type, interaction_type))
        s = len(interactions)
        if s == 0:
            self._class = 'No interaction'
            self._simple_class = 'No interaction'
        elif s == 1:
            score, chain_type, bond_type, interaction_type = interactions[0]
            self._class = '%s %s %s' % (chain_type, bond_type, interaction_type)
            self._simple_class = '%s interaction' % interaction_type
        else:
            i_str_parts = []
            max_score = 0.
            for score, chain_type, bond_type, interaction_type in interactions:
                i_str_parts.append('%s %s %s' % (chain_type, bond_type, interaction_type))
                if score > max_score:
                    max_score = score
                    self._simple_class = '%s interaction' % interaction_type
            self._class = 'Multi interaction: %s' % (' and '.join(i_str_parts))

    def getClass(self):
        try:
            if self._class is not None:
                return self._class, self._simple_class
            else:
                self.computeClass()
                return self._class, self._simple_class
        except AttributeError:
            self.computeClass()
            return self._class, self._simple_class


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
                    weight_sum += weight
                    degree_sum += float(degree) * weight
                    score_sum += score * weight
                if weight_sum == 0.:
                    continue
                average_degree = degree_sum / weight_sum
                average_score = score_sum / weight_sum
                average_profile.setDegree(chain_type, bond_type, interaction_type, average_degree)
                average_profile.setScore(chain_type, bond_type, interaction_type, average_score)
    return average_profile


def parse_rinerator_interaction_score_file(interaction_score_file: str, struct_ligands, metals, ions,
                                            chain_type_map, res_contig_map) -> list[tuple[str, int | str, str, int | str, str, str , float]]:
    f = gzip.open(interaction_score_file, 'rt')
    lines = f.readlines()
    f.close()

    interactions: list[tuple[str, int | str, str, int | str, str, str , float]] = []

    profiles: dict[str, residue_package.Residue_Map[Interaction_profile]] = {}
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}


    for line in lines[1:]:
        # A:39:_:ASP (cnt:mc_sc) A:41:_:LEU    2.220600
        if line.count('\t') == 0:
            print(f'{interaction_score_file=}')
            break
        words = line.split('\t')
        score = float(words[1])
        edge: str = words[0]
        res_a, interaction_type, res_b = edge.split()
        chain_a, res_nr_a, insertioncode_a, res_name_a = res_a.split(':')
        chain_b, res_nr_b, insertioncode_b, res_name_b = res_b.split(':')
        interaction_base_type, interaction_sub_type = interaction_type[1:-1].split(':')
        interaction_type_a, interaction_type_b = interaction_sub_type.split('_')

        res_nr_a = "%s%s" % (res_nr_a, insertioncode_a.replace('_', ''))
        try:
            res_nr_a = int(res_nr_a)
        except:
            pass
        res_nr_b = "%s%s" % (res_nr_b, insertioncode_b.replace('_', ''))
        try:
            res_nr_b = int(res_nr_b)
        except:
            pass

        if res_name_a in ligands.BORING_LIGANDS and res_name_a not in residues.THREE_TO_ONE:
            continue
        if res_name_b in ligands.BORING_LIGANDS and res_name_b not in residues.THREE_TO_ONE:
            continue

        interactions.append((chain_a, res_nr_a, chain_b, res_nr_b, interaction_base_type, interaction_type_a, score))
        interactions.append((chain_b, res_nr_b, chain_a, res_nr_a, interaction_base_type, interaction_type_b, score))
        
        if chain_a not in profiles:
            profiles[chain_a] = residue_package.Residue_Map()
        if not profiles[chain_a].contains(res_a):
            profiles[chain_a].add_item(res_a, Interaction_profile())
        profile = profiles[chain_a].get_item(res_a)
        add_interaction_to_profile(chain_a, res_a, chain_b, res_b, struct_ligands, metals, ions,
                                            chain_type_map, res_contig_map, profile,
                                            interaction_type_a, interaction_base_type, score, ligand_profiles,
                                            metal_profiles,
                                            ion_profiles,
                                            chain_chain_profiles)
        
        if chain_b not in profiles:
            profiles[chain_b] = residue_package.Residue_Map()
        if not profiles[chain_b].contains(res_b):
            profiles[chain_b].add_item(res_b, Interaction_profile())
        profile = profiles[chain_b].get_item(res_b)
        add_interaction_to_profile(chain_b, res_b, chain_a, res_a, struct_ligands, metals, ions,
                                            chain_type_map, res_contig_map, profile,
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

def getCentMap(centrality_file):

    f = gzip.open(centrality_file, 'rt')
    lines = f.read().split('\n')
    f.close()

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
            res_nr = "%s%s" % (res_nr, insertioncode.replace('_', ''))
            try:
                res_nr = int(res_nr)
            except:
                pass

            if chain not in centrality_map:
                centrality_map[chain] = residue_package.Residue_Map()

            centrality_map[chain].add_item(res_nr, cent_scores)
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
        mol_type_profile[chain] = residue_package.Residue_Map()
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
        ligands: set[str, int | str],
        metals: set[str, int | str],
        ions: set[str, int | str],
        chain_type_map: dict[str, str],
        res_contig_map: dict[str, residue_package.Residue_Map],
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
        if (chain, res) in ligands:
            add_to_mol_type_profile(chain, res, score, ligand_profiles)
        elif (chain, res) in metals:
            add_to_mol_type_profile(chain, res, score, metal_profiles)
        elif (chain, res) in ions:
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
        interaction_type = chain_type_map[chain_b]
    else:
        if chain not in res_contig_map:
            return
        if not res_contig_map[chain].contains(res_b):
            return
        res_dist = abs(res_contig_map[chain].get_item(res)[0] - res_contig_map[chain].get_item(res_b)[0])
        if res_dist < 2:
            interaction_type = 'neighbor'
        elif res_dist < 6:
            interaction_type = 'short'
        else:
            interaction_type = 'long'

    error = profile.addEdge(chaintype, bondtype, interaction_type, score, chain_b, res_b)
    if error is not None:
        return None


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
                    if not (chain, chain_b) in chain_chain_profiles:
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
        res_contig_map: dict[str, residue_package.Residue_Map[tuple[int, str]]],
        base_path: str,
        chain_type_map: dict[str, str],
        model_path : str | None = None,
        keep_tmp_files: bool = True
        ) -> tuple[
            dict[str, residue_package.Residue_Map[tuple[Interaction_profile, Centrality_scores]]],
            dict,
            dict,
            dict,
            dict,
            dict[str, residue_package.Residue_Map[dict[str, residue_package.Residue_Map[list[tuple[str, str, float]]]]]]
        ]:

    runtimes = []
    t0 = time.time()
    if config.verbosity >= 4:
        print(f'Call of rin.lookup: {structure_id=} {base_path=} {model_path=} {keep_tmp_files=}')

    folder_path, path_stem = get_entry_path(structure_id, config.rin_db_path, custom_db_path= config.custom_db_path, model_path=model_path, path_to_model_db=config.path_to_model_db, temp_folder=config.tmp_folder)

    if config.verbosity >= 4:
        print(f'Paths in rin.lookup: {structure_id=} {folder_path=} {path_stem=}')

    interaction_score_file = f"{path_stem}_intsc.ea.gz"

    if not os.path.isfile(interaction_score_file):
        if config.verbosity >= 3:
            print(f'Did not find RIN: {interaction_score_file}, folder path: {folder_path}' )

        if not os.path.isdir(folder_path):
            os.makedirs(folder_path)

        calcRIN(config, page.encode(), folder_path, path_stem, structure_id, config.rinerator_path, True, structure_path=model_path)

        interaction_score_file = "%s_intsc.ea.gz" % (path_stem)
        
        if not os.path.isfile(interaction_score_file):
            return f'RIN not found and calculation of RIN failed for: {structure_id=}, {model_path=} {interaction_score_file=}'

    network_file = "%s.sif.gz" % (path_stem)
    interaction_count_file = "%s_nrint.ea.gz" % (path_stem)
    residue_file = "%s_res.txt.gz" % (path_stem)
    centrality_file = "%s_btw_cent.txt.gz" % (path_stem)

    t1 = time.time()
    runtimes.append(t1-t0)
    if config.verbosity >= 4:
        print(f'Call of parse_rinerator_interaction_score_file: {structure_id} {interaction_score_file} {t1-t0}')

    interaction_list, profiles, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles = parse_rinerator_interaction_score_file(interaction_score_file, ligands, metals, ions, chain_type_map, res_contig_map)

    t2 = time.time()
    runtimes.append(t2-t1)
    if config.verbosity >= 4:
        print(f'Call of getCentMap: {structure_id}  {centrality_file} {t2-t1}')

    centrality_map: dict[str, residue_package.Residue_Map[Centrality_scores]] = {}

    if os.path.isfile(centrality_file):
        centrality_map = getCentMap(centrality_file)
    else:
        centrality_map = {}

    chains = list(res_contig_map.keys())

    t3 = time.time()
    runtimes.append(t3-t2)
    if config.verbosity >= 4:
        print(f'After getCentMap: {structure_id} {centrality_file} {t3-t2}')

    if len(chains) == 0:
        return f'Illegal input, no chains chosen in RIN lookup of {structure_id}'

    t4 = time.time()
    runtimes.append(t4-t3)
    if config.verbosity >= 4:
        print(f'Call of calculateIAPProfiles: {structure_id} {t4-t3}')

    if len(profiles) == 0:
        return f'Profiles_map is empty, in RIN lookup of {structure_id}'
    
    #ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles = calculateIAPProfiles(interaction_list, chains, ligands, metals, ions)

    if not keep_tmp_files:
        os.remove(interaction_score_file)
        os.remove(network_file)
        os.remove(interaction_count_file)
        os.remove(residue_file)
        os.remove(centrality_file)

    t5 = time.time()
    runtimes.append(t5-t4)
    if config.verbosity >= 4:
        print(f'At end of rin.lookup: {structure_id} {t5-t4} {t5-t0}')

    return profiles, centrality_map, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, interaction_list, runtimes
