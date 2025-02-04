import gzip
import os
import sys
import traceback

from structman.lib.createRINdb import calcRIN
from structman.lib.sdsc.consts import residues, ligands

def encode_float(fl_num, precision = 5):
    return f'{fl_num:.{precision}f}'


class Ligand_types:
    __slots__ = ['ligand', 'ion', 'metal']

    def __init__(self, lig_str=None):
        if lig_str is None:
            self.ligand = [0, 0.0]
            self.ion = [0, 0.0]
            self.metal = [0, 0.0]
        else:
            nd, ns, sd, ss, ld, ls = lig_str.split('(')
            self.ligand = [int(float(nd)), float(ns)]
            self.ion = [int(float(sd)), float(ss)]
            self.metal = [int(float(ld)), float(ls)]

    def addEdge(self, interaction_type, score):
        if interaction_type == 'ligand':
            self.ligand[0] += 1
            self.ligand[1] += score
            return None
        elif interaction_type == 'ion':
            self.ion[0] += 1
            self.ion[1] += score
            return None
        elif interaction_type == 'metal':
            self.metal[0] += 1
            self.metal[1] += score
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type

    def encode(self):
        return str(self.ligand[0]) + '(' + encode_float(self.ligand[1]) + '(' + str(self.ion[0]) + '(' + encode_float(self.ion[1]) + '(' + str(self.metal[0]) + '(' + encode_float(self.metal[1])

    def getScore(self, interaction_type):
        if interaction_type == 'ligand':
            return self.ligand[1]
        elif interaction_type == 'ion':
            return self.ion[1]
        elif interaction_type == 'metal':
            return self.metal[1]
        else:
            return 0.0

    def getDegree(self, interaction_type):
        if interaction_type == 'ligand':
            return self.ligand[0]
        elif interaction_type == 'ion':
            return self.ion[0]
        elif interaction_type == 'metal':
            return self.metal[0]
        else:
            return 0

    def setScore(self, interaction_type, score):
        if interaction_type == 'ligand':
            self.ligand[1] = score
            return None
        elif interaction_type == 'ion':
            self.ion[1] = score
            return None
        elif interaction_type == 'metal':
            self.metal[1] = score
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type

    def setDegree(self, interaction_type, degree):
        if interaction_type == 'ligand':
            self.ligand[0] = degree
            return None
        elif interaction_type == 'ion':
            self.ion[0] = degree
            return None
        elif interaction_type == 'metal':
            self.metal[0] = degree
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type


ligand_types = set(['ligand', 'ion', 'metal'])


class Intrachain_types:
    __slots__ = ['neighbor', 'short', 'long']

    def __init__(self, intra_str=None):
        if intra_str is None:
            self.neighbor = [0, 0.0]
            self.short = [0, 0.0]
            self.long = [0, 0.0]
        else:
            nd, ns, sd, ss, ld, ls = intra_str.split('(')
            self.neighbor = [int(float(nd)), float(ns)]
            self.short = [int(float(sd)), float(ss)]
            self.long = [int(float(ld)), float(ls)]

    def addEdge(self, interaction_type, score):
        if interaction_type == 'neighbor':
            self.neighbor[0] += 1
            self.neighbor[1] += score
            return None
        elif interaction_type == 'short':
            self.short[0] += 1
            self.short[1] += score
            return None
        elif interaction_type == 'long':
            self.long[0] += 1
            self.long[1] += score
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type

    def encode(self):
        return str(self.neighbor[0]) + '(' + encode_float(self.neighbor[1]) + '(' + str(self.short[0]) + '(' + encode_float(self.short[1]) + '(' + str(self.long[0]) + '(' + encode_float(self.long[1])

    def getScore(self, interaction_type):
        if interaction_type == 'neighbor':
            return self.neighbor[1]
        elif interaction_type == 'short':
            return self.short[1]
        elif interaction_type == 'long':
            return self.long[1]
        else:
            return 0.0

    def getDegree(self, interaction_type):
        if interaction_type == 'neighbor':
            return self.neighbor[0]
        elif interaction_type == 'short':
            return self.short[0]
        elif interaction_type == 'long':
            return self.long[0]
        else:
            return 0

    def setScore(self, interaction_type, score):
        if interaction_type == 'neighbor':
            self.neighbor[1] = score
            return None
        elif interaction_type == 'short':
            self.short[1] = score
            return None
        elif interaction_type == 'long':
            self.long[1] = score
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type

    def setDegree(self, interaction_type, degree):
        if interaction_type == 'neighbor':
            self.neighbor[0] = degree
            return None
        elif interaction_type == 'short':
            self.short[0] = degree
            return None
        elif interaction_type == 'long':
            self.long[0] = degree
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type


intrachain_types = set(['neighbor', 'short', 'long'])


class Interchain_types:
    __slots__ = ['Protein', 'DNA', 'RNA', 'Peptide']

    def __init__(self, interchain_str=None):
        if interchain_str is None:
            self.Protein = [0, 0.0]
            self.DNA = [0, 0.0]
            self.RNA = [0, 0.0]
            self.Peptide = [0, 0.0]
        else:
            nd, ns, sd, ss, ld, ls, pd, ps = interchain_str.split('(')
            self.Protein = [int(float(nd)), float(ns)]
            self.DNA = [int(float(sd)), float(ss)]
            self.RNA = [int(float(ld)), float(ls)]
            self.Peptide = [int(float(pd)), float(ps)]

    def addEdge(self, interaction_type, score):
        if interaction_type == 'Protein':
            self.Protein[0] += 1
            self.Protein[1] += score
            return None
        elif interaction_type == 'DNA':
            self.DNA[0] += 1
            self.DNA[1] += score
            return None
        elif interaction_type == 'RNA':
            self.RNA[0] += 1
            self.RNA[1] += score
            return None
        elif interaction_type == 'Peptide':
            self.Peptide[0] += 1
            self.Peptide[1] += score
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type

    def encode(self):
        return str(self.Protein[0]) + '(' + encode_float(self.Protein[1]) + '(' + str(self.DNA[0]) + '(' + encode_float(self.DNA[1]) + '(' + str(self.RNA[0]) + '(' + encode_float(self.RNA[1]) + '(' + str(self.Peptide[0]) + '(' + encode_float(self.Peptide[1])

    def getScore(self, interaction_type):
        if interaction_type == 'Protein':
            return self.Protein[1]
        elif interaction_type == 'DNA':
            return self.DNA[1]
        elif interaction_type == 'RNA':
            return self.RNA[1]
        elif interaction_type == 'Peptide':
            return self.Peptide[1]
        else:
            return 0.0

    def getDegree(self, interaction_type):
        if interaction_type == 'Protein':
            return self.Protein[0]
        elif interaction_type == 'DNA':
            return self.DNA[0]
        elif interaction_type == 'RNA':
            return self.RNA[0]
        elif interaction_type == 'Peptide':
            return self.Peptide[0]
        else:
            return 0

    def setScore(self, interaction_type, score):
        if interaction_type == 'Protein':
            self.Protein[1] = score
            return None
        elif interaction_type == 'DNA':
            self.DNA[1] = score
            return None
        elif interaction_type == 'RNA':
            self.RNA[1] = score
            return None
        elif interaction_type == 'Peptide':
            self.Peptide[1] = score
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type

    def setDegree(self, interaction_type, degree):
        if interaction_type == 'Protein':
            self.Protein[0] = degree
            return None
        elif interaction_type == 'DNA':
            self.DNA[0] = degree
            return None
        elif interaction_type == 'RNA':
            self.RNA[0] = degree
            return None
        elif interaction_type == 'Peptide':
            self.Peptide[0] = degree
            return None
        else:
            return 'Unknown interactiontype: %s' % interaction_type


interchain_types = set(['Protein', 'DNA', 'RNA', 'Peptide'])

interaction_types = set(['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide'])

non_intra_interaction_types = set(['ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide'])

chain_type_interactions = set(['Protein', 'DNA', 'RNA', 'Peptide'])

small_molecule_interactions = set(['ligand', 'ion', 'metal'])


class Interaction_type_profile:
    __slots__ = ['sw_ligand', 'interchain', 'intrachain']

    def __init__(self, inter_str=None):
        if inter_str is None:
            self.sw_ligand = Ligand_types()
            self.interchain = Interchain_types()
            self.intrachain = Intrachain_types()
        else:
            lig_str, intra_str, interchain_str = inter_str.split('|')
            self.sw_ligand = Ligand_types(lig_str=lig_str)
            self.interchain = Interchain_types(interchain_str=interchain_str)
            self.intrachain = Intrachain_types(intra_str=intra_str)

    def addEdge(self, interaction_type, score):
        if interaction_type in interchain_types:
            return self.interchain.addEdge(interaction_type, score)
        elif interaction_type in intrachain_types:
            return self.intrachain.addEdge(interaction_type, score)
        elif interaction_type in ligand_types:
            return self.sw_ligand.addEdge(interaction_type, score)
        else:
            return 'Unknown interactiontype: %s' % interaction_type

    def encode(self):
        return self.sw_ligand.encode() + '|' + self.intrachain.encode() + '|' + self.interchain.encode()

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


class Chain_types:
    __slots__ = ['contact', 'h_bond', 'overlap']

    def __init__(self, chain_str=None):
        if chain_str is None:
            self.contact = Interaction_type_profile()
            self.h_bond = Interaction_type_profile()
            self.overlap = Interaction_type_profile()
        else:
            con_str, hb_str, ovl_str = chain_str.split(';')
            self.contact = Interaction_type_profile(inter_str=con_str)
            self.h_bond = Interaction_type_profile(inter_str=hb_str)
            self.overlap = Interaction_type_profile(inter_str=ovl_str)

    def addEdge(self, bondtype, interaction_type, score):
        if bondtype == 'cnt':
            return self.contact.addEdge(interaction_type, score)
        elif bondtype == 'hbond':
            return self.h_bond.addEdge(interaction_type, score)
        elif bondtype == 'ovl':
            return self.overlap.addEdge(interaction_type, score)
        else:
            return 'Unknown bondtype: %s' % bondtype

    def encode(self):
        return self.contact.encode() + ';' + self.h_bond.encode() + ';' + self.overlap.encode()

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


class Interaction_profile:
    __slots__ = ['mainchain', 'sidechain', '_class', '_simple_class', 'interacting_chains', 'interacting_ligands']

    def __init__(self, profile_str=None, interacting_chains_str=None, interacting_ligands_str=None):
        if profile_str is None:
            self.mainchain = Chain_types()
            self.sidechain = Chain_types()
            self._class = None
            self._simple_class = None
        else:
            mc_str, sc_str = profile_str.split(',')
            self.mainchain = Chain_types(chain_str=mc_str)
            self.sidechain = Chain_types(chain_str=sc_str)
            self._class = None
            self._simple_class = None

        if interacting_chains_str is not None:
            self.interacting_chains = interacting_chains_str.split(',')
        else:
            self.interacting_chains = set()

        if interacting_ligands_str is not None:
            self.interacting_ligands = set()
            for id_tuple in interacting_ligands_str.split(','):
                self.interacting_ligands.add(id_tuple)
        else:
            self.interacting_ligands = set()

    def addEdge(self, chaintype, bondtype, interaction_type, score, chain_id, res_id):
        if bondtype != 'ovl':
            if interaction_type in small_molecule_interactions:
                self.interacting_ligands.add('%s:%s' % (chain_id, res_id))
            elif interaction_type in chain_type_interactions:
                self.interacting_chains.add(chain_id)

        if chaintype == 'mc':
            return self.mainchain.addEdge(bondtype, interaction_type, score)
        elif chaintype == 'sc' or chaintype == 'ligand':
            return self.sidechain.addEdge(bondtype, interaction_type, score)
        else:
            return 'Unknown chaintype in addEdge: %s' % chaintype

    def encode(self):
        return self.mainchain.encode() + ',' + self.sidechain.encode()

    def getScore(self, chain_type, bond_type, interaction_type):
        if chain_type == 'mc':
            return self.mainchain.getScore(bond_type, interaction_type)
        elif chain_type == 'sc':
            return self.sidechain.getScore(bond_type, interaction_type)
        else:
            return 'Unknown chaintype in getScore: %s' % chain_type

    def getDegree(self, chain_type, bond_type, interaction_type):
        if chain_type == 'mc':
            return self.mainchain.getDegree(bond_type, interaction_type)
        elif chain_type == 'sc':
            return self.sidechain.getDegree(bond_type, interaction_type)
        else:
            return 'Unknown chaintype in getDegree: %s' % chain_type

    def setScore(self, chain_type, bond_type, interaction_type, score):
        if chain_type == 'mc':
            return self.mainchain.setScore(bond_type, interaction_type, score)
        elif chain_type == 'sc':
            return self.sidechain.setScore(bond_type, interaction_type, score)
        else:
            return 'Unknown chaintype in setScore: %s' % chain_type

    def setDegree(self, chain_type, bond_type, interaction_type, degree):
        if chain_type == 'mc':
            return self.mainchain.setDegree(bond_type, interaction_type, degree)
        elif chain_type == 'sc':
            return self.sidechain.setDegree(bond_type, interaction_type, degree)
        else:
            return 'Unknown chaintype in setDegree: %s' % chain_type

    def getChainSpecificCombiScore(self, chain_type, interaction_type):
        if chain_type == 'mc':
            return self.mainchain.getCombiScore(interaction_type)
        elif chain_type == 'sc':
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
            return self.mainchain.getCombiDegree(interaction_type)
        elif chain_type == 'sc':
            return self.sidechain.getCombiDegree(interaction_type)
        else:
            return None

    def computeClass(self):
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
        if self._class is not None:
            return self._class, self._simple_class
        else:
            self.computeClass()
            return self._class, self._simple_class


class Centrality_scores:
    __slots__ = ['AbsoluteCentrality', 'LengthNormalizedCentrality', 'MinMaxNormalizedCentrality', 'AbsoluteCentralityWithNegative',
                 'LengthNormalizedCentralityWithNegative', 'MinMaxNormalizedCentralityWithNegative', 'AbsoluteComplexCentrality',
                 'LengthNormalizedComplexCentrality', 'MinMaxNormalizedComplexCentrality', 'AbsoluteComplexCentralityWithNegative',
                 'LengthNormalizedComplexCentralityWithNegative', 'MinMaxNormalizedComplexCentralityWithNegative', 'cent_list']
    
    feature_names = ['AbsoluteCentrality', 'LengthNormalizedCentrality', 'MinMaxNormalizedCentrality', 'AbsoluteCentralityWithNegative',
                 'LengthNormalizedCentralityWithNegative', 'MinMaxNormalizedCentralityWithNegative', 'AbsoluteComplexCentrality',
                 'LengthNormalizedComplexCentrality', 'MinMaxNormalizedComplexCentrality', 'AbsoluteComplexCentralityWithNegative',
                 'LengthNormalizedComplexCentralityWithNegative', 'MinMaxNormalizedComplexCentralityWithNegative']

    def __init__(self, AbsoluteCentrality=None, LengthNormalizedCentrality=None, MinMaxNormalizedCentrality=None,
                 AbsoluteCentralityWithNegative=None, LengthNormalizedCentralityWithNegative=None, MinMaxNormalizedCentralityWithNegative=None,
                 AbsoluteComplexCentrality=None, LengthNormalizedComplexCentrality=None, MinMaxNormalizedComplexCentrality=None,
                 AbsoluteComplexCentralityWithNegative=None, LengthNormalizedComplexCentralityWithNegative=None,
                 MinMaxNormalizedComplexCentralityWithNegative=None, code_str=None, cent_list=None):

        if code_str is not None:
            self.str_decode(code_str)
            return

        if cent_list is not None:
            self.cent_list = cent_list
            self.setAllByCentList()
            return

        self.cent_list = []
        self.AbsoluteCentrality = AbsoluteCentrality
        self.cent_list.append(self.AbsoluteCentrality)

        self.LengthNormalizedCentrality = LengthNormalizedCentrality
        self.cent_list.append(self.LengthNormalizedCentrality)

        self.MinMaxNormalizedCentrality = MinMaxNormalizedCentrality
        self.cent_list.append(self.MinMaxNormalizedCentrality)

        self.AbsoluteCentralityWithNegative = AbsoluteCentralityWithNegative
        self.cent_list.append(self.AbsoluteCentralityWithNegative)

        self.LengthNormalizedCentralityWithNegative = LengthNormalizedCentralityWithNegative
        self.cent_list.append(self.LengthNormalizedCentralityWithNegative)

        self.MinMaxNormalizedCentralityWithNegative = MinMaxNormalizedCentralityWithNegative
        self.cent_list.append(self.MinMaxNormalizedCentralityWithNegative)

        self.AbsoluteComplexCentrality = AbsoluteComplexCentrality
        self.cent_list.append(self.AbsoluteComplexCentrality)

        self.LengthNormalizedComplexCentrality = LengthNormalizedComplexCentrality
        self.cent_list.append(self.LengthNormalizedComplexCentrality)

        self.MinMaxNormalizedComplexCentrality = MinMaxNormalizedComplexCentrality
        self.cent_list.append(self.MinMaxNormalizedComplexCentrality)

        self.AbsoluteComplexCentralityWithNegative = AbsoluteComplexCentralityWithNegative
        self.cent_list.append(self.AbsoluteComplexCentralityWithNegative)

        self.LengthNormalizedComplexCentralityWithNegative = LengthNormalizedComplexCentralityWithNegative
        self.cent_list.append(self.LengthNormalizedComplexCentralityWithNegative)

        self.MinMaxNormalizedComplexCentralityWithNegative = MinMaxNormalizedComplexCentralityWithNegative
        self.cent_list.append(self.MinMaxNormalizedComplexCentralityWithNegative)

    def str_encode(self):
        str_code = ';'.join([str(x)[:10] for x in self.cent_list])
        return str_code

    def str_decode(self, str_code):
        self.cent_list = []
        for x in str_code.split(';'):
            try:
                self.cent_list.append(float(x))
            except:
                self.cent_list.append(None)
        self.setAllByCentList()

    def setAllByCentList(self):
        self.AbsoluteCentrality = self.cent_list[0]

        self.LengthNormalizedCentrality = self.cent_list[1]

        self.MinMaxNormalizedCentrality = self.cent_list[2]

        self.AbsoluteCentralityWithNegative = self.cent_list[3]

        self.LengthNormalizedCentralityWithNegative = self.cent_list[4]

        self.MinMaxNormalizedCentralityWithNegative = self.cent_list[5]

        self.AbsoluteComplexCentrality = self.cent_list[6]

        self.LengthNormalizedComplexCentrality = self.cent_list[7]

        self.MinMaxNormalizedComplexCentrality = self.cent_list[8]

        self.AbsoluteComplexCentralityWithNegative = self.cent_list[9]

        self.LengthNormalizedComplexCentralityWithNegative = self.cent_list[10]

        self.MinMaxNormalizedComplexCentralityWithNegative = self.cent_list[11]


# called by database
# called by sdsc
def calculateAverageProfile(profiles):
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


def getIAmap(interaction_score_file):
    f = gzip.open(interaction_score_file, 'rt')
    lines = f.readlines()
    f.close()

    IAmap = {}

    for line in lines[1:]:
        # A:39:_:ASP (cnt:mc_sc) A:41:_:LEU    2.220600
        if line.count('\t') == 0:
            print(interaction_score_file)
            break
        score = float(line.split('\t')[1])
        edge = line.split('\t')[0]
        res_a, interaction_type, res_b = edge.split()
        chain_a, res_nr_a, insertioncode_a, res_name_a = res_a.split(':')
        chain_b, res_nr_b, insertioncode_b, res_name_b = res_b.split(':')
        interaction_base_type, interaction_sub_type = interaction_type[1:-1].split(':')
        interaction_type_a, interaction_type_b = interaction_sub_type.split('_')

        res_nr_a = "%s%s" % (res_nr_a, insertioncode_a.replace('_', ''))
        res_nr_b = "%s%s" % (res_nr_b, insertioncode_b.replace('_', ''))
        if res_name_a in ligands.BORING_LIGANDS and res_name_a not in residues.THREE_TO_ONE:
            continue
        if res_name_b in ligands.BORING_LIGANDS and res_name_b not in residues.THREE_TO_ONE:
            continue

        if chain_a not in IAmap:
            IAmap[chain_a] = {}
        if res_nr_a not in IAmap[chain_a]:
            IAmap[chain_a][res_nr_a] = {}
        if interaction_base_type not in IAmap[chain_a][res_nr_a]:
            IAmap[chain_a][res_nr_a][interaction_base_type] = {}
        if interaction_type_a not in IAmap[chain_a][res_nr_a][interaction_base_type]:
            IAmap[chain_a][res_nr_a][interaction_base_type][interaction_type_a] = {}
        if not chain_b in IAmap[chain_a][res_nr_a][interaction_base_type][interaction_type_a]:
            IAmap[chain_a][res_nr_a][interaction_base_type][interaction_type_a][chain_b] = {}
        IAmap[chain_a][res_nr_a][interaction_base_type][interaction_type_a][chain_b][res_nr_b] = score

        if chain_b not in IAmap:
            IAmap[chain_b] = {}
        if res_nr_b not in IAmap[chain_b]:
            IAmap[chain_b][res_nr_b] = {}
        if interaction_base_type not in IAmap[chain_b][res_nr_b]:
            IAmap[chain_b][res_nr_b][interaction_base_type] = {}
        if interaction_type_b not in IAmap[chain_b][res_nr_b][interaction_base_type]:
            IAmap[chain_b][res_nr_b][interaction_base_type][interaction_type_b] = {}
        if not chain_a in IAmap[chain_b][res_nr_b][interaction_base_type][interaction_type_b]:
            IAmap[chain_b][res_nr_b][interaction_base_type][interaction_type_b][chain_a] = {}
        IAmap[chain_b][res_nr_b][interaction_base_type][interaction_type_b][chain_a][res_nr_a] = score

    return IAmap


# Takes a dict of structure_id, chain : IAmap and calculates an aggregated RIN for the given chains.
# The chains must assign the same protein in different structures, therefor residue-position backmaps for all of the are needed.
# Called by serializedPipeline

def aggregate_IAmaps(IAmaps, backmaps):
    aggregated_IAmap_precursor = {}
    max_edge_count = {}
    for structure_id, chain in IAmaps:
        IAmap = IAmaps[(structure_id, chain)]
        backmap = backmaps[structure_id]
        if chain not in IAmap:
            continue
        for res_nr in IAmap[chain]:
            if not res_nr in backmap:
                continue
            position = backmap[res_nr]
            if position is None:
                continue
            for res_nr_2 in IAmap[chain]:
                if res_nr_2 not in backmap:
                    continue
                i_position = backmap[res_nr_2]
                if i_position is None:
                    continue
                if (position, i_position) not in max_edge_count:
                    max_edge_count[(position, i_position)] = 1
                else:
                    max_edge_count[(position, i_position)] += 1
            for interaction_base_type in IAmap[chain][res_nr]:
                for interaction_type in IAmap[chain][res_nr][interaction_base_type]:
                    for (i_chain, i_res_nr) in IAmap[chain][res_nr][interaction_base_type][interaction_type]:
                        if i_chain != chain: #ignore non-intrachain interactions
                            continue
                        if not i_res_nr in backmap:
                            continue
                        i_position = backmap[i_res_nr]
                        if i_position is None:
                            continue
                        if not (position, i_position) in aggregated_IAmap_precursor:
                            aggregated_IAmap_precursor[(position, i_position)] = {}
                        if not interaction_base_type in aggregated_IAmap_precursor[(position, i_position)]:
                            aggregated_IAmap_precursor[(position, i_position)][interaction_base_type] = {}
                        if not interaction_type in aggregated_IAmap_precursor[(position, i_position)][interaction_base_type]:
                            aggregated_IAmap_precursor[(position, i_position)][interaction_base_type][interaction_type] = []
                        aggregated_IAmap_precursor[(position, i_position)][interaction_base_type][interaction_type].append(IAmap[chain][res_nr][interaction_base_type][interaction_type][(i_chain, i_res_nr)])

    aggregate_IAmap = {}
    for (position, i_position) in aggregated_IAmap_precursor:
        max_edge = max_edge_count[(position, i_position)]
        if not (position, i_position) in aggregate_IAmap:
            aggregate_IAmap[(position, i_position)] = {}
        for interaction_base_type in aggregated_IAmap_precursor[(position, i_position)]:
            if not interaction_base_type in aggregate_IAmap[(position, i_position)]:
                aggregate_IAmap[(position, i_position)][interaction_base_type] = {}
            for interaction_type in aggregated_IAmap_precursor[(position, i_position)][interaction_base_type]:
                if not interaction_type in aggregate_IAmap[(position, i_position)][interaction_base_type]:
                    aggregate_IAmap[(position, i_position)][interaction_base_type][interaction_type] = {}
                aggregate_IAmap[(position, i_position)][interaction_base_type][interaction_type] = sum(aggregated_IAmap_precursor[(position, i_position)][interaction_base_type][interaction_type])/max_edge

    return aggregate_IAmap

def float_or_none(string):
    if string == 'None':
        return None
    else:
        return float(string)


def getCentMap(centrality_file):

    f = gzip.open(centrality_file, 'rt')
    lines = f.read().split('\n')
    f.close()

    centrality_map = {}

    for line in lines[1:]:
        # Residue        AbsoluteCentrality      LengthNormalizedCentrality      MinMaxNormalizedCentrality      AbsoluteCentralityWithNegative  LengthNormalizedCentralityWithNegative  MinMaxNormalizedCentralityWithNegative  AbsoluteComplexCentrality       LengthNormalizedComplexCentrality       MinMaxNormalizedComplexCentrality       AbsoluteComplexCentralityWithNegative   LengthNormalizedComplexCentralityWithNegative   MinMaxNormalizedComplexCentralityWithNegative
        # A:20:_:TYR      1927.0  0.24084489438820147     1.0     1947.0  0.24334458192725908     1.0     2412.0  0.07389705882352941     1.2516865594187856      2260.0  0.06924019607843138     1.1607601438109914

        # A:9:_:TYR    1036.0    0.150079675503    0.829463570857
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

            cent_scores = Centrality_scores(AbsoluteCentrality=float_or_none(words[1]),
                                            LengthNormalizedCentrality=float_or_none(words[2]),
                                            MinMaxNormalizedCentrality=float_or_none(words[3]),
                                            AbsoluteCentralityWithNegative=float_or_none(words[4]),
                                            LengthNormalizedCentralityWithNegative=float_or_none(words[5]),
                                            MinMaxNormalizedCentralityWithNegative=float_or_none(words[6]),
                                            AbsoluteComplexCentrality=float_or_none(words[7]),
                                            LengthNormalizedComplexCentrality=float_or_none(words[8]),
                                            MinMaxNormalizedComplexCentrality=float_or_none(words[9]),
                                            AbsoluteComplexCentralityWithNegative=float_or_none(words[10]),
                                            LengthNormalizedComplexCentralityWithNegative=float_or_none(words[11]),
                                            MinMaxNormalizedComplexCentralityWithNegative=float_or_none(words[12]))

            [chain, res_nr, insertioncode, res_name] = res.split(':')
            res_nr = "%s%s" % (res_nr, insertioncode.replace('_', ''))

            if chain not in centrality_map:
                centrality_map[chain] = {}

            centrality_map[chain][res_nr] = cent_scores
        except:
            # some of the rins are not updated, this leads to an error here,-> simply update all rins
            print('Error in getCentMap: ', centrality_file)
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            print('\n'.join([str(e), str(f), str(g)]))
            continue

    return centrality_map


def getProfile(interaction_map, residue, ligands, metals, ions, res_contig_map, pdb_id, chain_type_map):
    profile = Interaction_profile()
    (chain, res) = residue

    if chain not in     interaction_map:
        return profile
    if res not in interaction_map[chain]:
        return profile
    if chain not in res_contig_map:
        return profile
    if res not in res_contig_map[chain]:
        return profile

    for bondtype in interaction_map[chain][res]:
        if bondtype == 'combi':
            continue
        for chaintype in interaction_map[chain][res][bondtype]:
            for chain_b in interaction_map[chain][res][bondtype][chaintype]:
                for res_b in interaction_map[chain][res][bondtype][chaintype][chain_b]:
                    score = interaction_map[chain][res][bondtype][chaintype][chain_b][res_b]

                    if (chain_b, res_b) in ligands:
                        interaction_type = 'ligand'
                    elif (chain_b, res_b) in metals:
                        interaction_type = 'metal'
                    elif (chain_b, res_b) in ions:
                        interaction_type = 'ion'
                    elif chain != chain_b:
                        if chain_b not in chain_type_map:
                            continue
                        interaction_type = chain_type_map[chain_b]
                    else:
                        if res_b not in res_contig_map[chain]:
                            continue
                        res_dist = abs(res_contig_map[chain][res][0] - res_contig_map[chain][res_b][0])
                        if res_dist < 2:
                            interaction_type = 'neighbor'
                        elif res_dist < 6:
                            interaction_type = 'short'
                        else:
                            interaction_type = 'long'
                    error = profile.addEdge(chaintype, bondtype, interaction_type, score, chain_b, res_b)
                    if error is not None:
                        print(error, chain, res, pdb_id)
                        return None
    return profile


def calculateIAPProfiles(interaction_map, chains, ligands, metals, ions):
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}
    for chain, res in ligands:
        deg = 0
        score = 0.0
        if chain not in interaction_map:
            # This can happen if the chain consists only of ligands with at least one non-boring ligand
            # and this non-boring ligand only has bonds with boring ligands
            continue
        if res in interaction_map[chain]:
            for chain_b in interaction_map[chain][res]['combi']['all']:
                for res_b in interaction_map[chain][res]['combi']['all'][chain_b]:
                    score += interaction_map[chain][res]['combi']['all'][chain_b][res_b]
                    deg += 1
        ligand_profiles[(chain, res)] = [deg, score]

    for chain, res in metals:
        deg = 0
        score = 0.0
        if chain in interaction_map:
            if res in interaction_map[chain]:
                for chain_b in interaction_map[chain][res]['combi']['all']:
                    for res_b in interaction_map[chain][res]['combi']['all'][chain_b]:
                        score += interaction_map[chain][res]['combi']['all'][chain_b][res_b]
                        deg += 1
        metal_profiles[(chain, res)] = [deg, score]

    for chain, res in ions:
        deg = 0
        score = 0.0
        if chain in interaction_map:
            if res in interaction_map[chain]:
                for chain_b in interaction_map[chain][res]['combi']['all']:
                    for res_b in interaction_map[chain][res]['combi']['all'][chain_b]:
                        score += interaction_map[chain][res]['combi']['all'][chain_b][res_b]
                        deg += 1
        ion_profiles[(chain, res)] = [deg, score]

    for chain in chains:
        if chain in interaction_map:
            for res in interaction_map[chain]:
                if (chain, res) in ligands:
                    continue
                if (chain, res) in metals:
                    continue
                if (chain, res) in ions:
                    continue
                for chain_b in interaction_map[chain][res]['combi']['all']:
                    if chain == chain_b:
                        continue
                    for res_b in interaction_map[chain][res]['combi']['all'][chain_b]:
                        if (chain_b, res_b) in ligands:
                            continue
                        if (chain_b, res_b) in metals:
                            continue
                        if (chain_b, res_b) in ions:
                            continue
                        if not (chain, chain_b) in chain_chain_profiles:
                            chain_chain_profiles[(chain, chain_b)] = [0, 0.0]

                        score = interaction_map[chain][res]['combi']['all'][chain_b][res_b]

                        chain_chain_profiles[(chain, chain_b)][0] += 1
                        chain_chain_profiles[(chain, chain_b)][1] += score

    return ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles


# called by templateFiltering


def lookup(pdb_id, page, config, inp_residues, chains, ligands, metals, ions, res_contig_map, base_path, chain_type_map, encoded=True, model_path=None, keep_tmp_files=True):

    #TODO check if needed
    if not config.custom_db_path:
        pdb_id = pdb_id.replace('_AU', '').lower().split('.')[0]

    # set path for custom db rin calculation
    if config.custom_db_path:
        folder_path = config.temp_folder
        custom_files_list = os.listdir(config.temp_folder)     
        path_stem = "%s/%s" % (folder_path,pdb_id)

    elif model_path is None:
        folder_path = "%s/%s/%s" % (base_path, pdb_id[1:-1], pdb_id)
        path_stem = "%s/%s" % (folder_path, pdb_id)
    else:
        folder_path = config.temp_folder
        if model_path[-4:] == '.pdb':
            path_stem = model_path[:-4] # remove .pdb from the end of the model path
        else:
            path_stem = model_path[:-7] # .pdb.gz
        pdb_id = path_stem.split('/')[-1].split('.')[0]

    interaction_score_file = "%s_intsc.ea.gz" % (path_stem)
    if not os.path.isfile(interaction_score_file):
        if config.verbosity >= 3:
            print(f'Did not find RIN: {interaction_score_file}, folder path: {folder_path}' )
        rinerator_path = config.rinerator_path

        if os.path.exists(folder_path):
            rin_target_path = folder_path         
        else:
            if model_path is None:
                if os.path.exists(base_path):
                    if not os.path.exists(f'{base_path}/{pdb_id[1:-1]}'):
                        os.mkdir(f'{base_path}/{pdb_id[1:-1]}')
                    os.mkdir(f'{base_path}/{pdb_id[1:-1]}/{pdb_id}')
            if os.path.exists(folder_path):
                rin_target_path = folder_path
            else:
                rin_target_path = config.tmp_folder

        calcRIN(config, page.encode(), rin_target_path, pdb_id, rinerator_path, True, structure_path=model_path)
        path_stem = "%s/%s" % (rin_target_path, pdb_id)
        interaction_score_file = "%s_intsc.ea.gz" % (path_stem)
        
        if not os.path.isfile(interaction_score_file):
            return f'RIN not found and calculation of RIN failed for: {pdb_id}, {model_path} {interaction_score_file}'

    network_file = "%s.sif.gz" % (path_stem)
    interaction_count_file = "%s_nrint.ea.gz" % (path_stem)
    residue_file = "%s_res.txt.gz" % (path_stem)
    centrality_file = "%s_btw_cent.txt.gz" % (path_stem)

    interaction_map = getIAmap(interaction_score_file)
    if os.path.isfile(centrality_file):
        centrality_map = getCentMap(centrality_file)
    else:
        centrality_map = {}
    profiles_map = {}
    if chains is None:
        chains = list(res_contig_map.keys())
    if len(chains) == 0:
        return 'Illegal input, no chains chosen in RIN lookup of %s' % pdb_id
    for chain in chains:
        profile_map = {}
        if inp_residues is None:
            residues = res_contig_map[chain]
        else:
            residues = inp_residues
        for res in residues:
            res_name = res_contig_map[chain][res][1]
            if len(res_name) < 3:
                continue
            if chain in centrality_map:
                if res in centrality_map[chain]:
                    centrality_scores = centrality_map[chain][res]
                else:
                    centrality_scores = None
            else:
                centrality_scores = None

            profile = getProfile(interaction_map, (chain, res), ligands, metals, ions, res_contig_map, pdb_id, chain_type_map)
            if profile is None:
                print(pdb_id, res, res_name)
                continue
            if encoded:
                profile_map[res] = profile.encode(), centrality_scores
            else:
                profile_map[res] = profile, centrality_scores
        profiles_map[chain] = profile_map
    if len(profiles_map) == 0:
        return 'Profiles_map is empty, in RIN lookup of %s' % pdb_id
    ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles = calculateIAPProfiles(interaction_map, chains, ligands, metals, ions)

    if not keep_tmp_files:
        os.remove(interaction_score_file)
        os.remove(network_file)
        os.remove(interaction_count_file)
        os.remove(residue_file)
        os.remove(centrality_file)

    return profiles_map, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, interaction_map
