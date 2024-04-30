# sdsc: structman datastructures and classes
from structman.lib.rin import Interaction_profile, Centrality_scores
from structman.lib.sdsc.consts.residues import METAL_ATOMS, ION_ATOMS
from structman.lib.sdsc.sdsc_utils import rin_classify, doomsday_protocol


class Residue(object):
    __slots__ = ['res_num', 'aa', 'lig_dist_str', 'lig_dists', 'chain_dist_str', 'chain_distances', 'RSA',
                 'relative_main_chain_acc', 'relative_side_chain_acc', 'SSA',
                 'homo_dist_str', 'homomer_distances',
                 'interaction_profile', 'interaction_profile_str', 'centrality_score_str', 'centralities', 'modres', 'b_factor', 'database_id', 'stored', 'phi', 'psi',
                 'intra_ssbond', 'inter_ssbond', 'ssbond_length', 'intra_link', 'inter_link', 'link_length', 'cis_conformation', 'cis_follower',
                 'inter_chain_median_kd', 'inter_chain_dist_weighted_kd', 'inter_chain_median_rsa',
                 'inter_chain_dist_weighted_rsa', 'intra_chain_median_kd', 'intra_chain_dist_weighted_kd',
                 'inter_chain_interactions_median', 'inter_chain_interactions_dist_weighted',
                 'intra_chain_interactions_median', 'intra_chain_interactions_dist_weighted',
                 'intra_chain_median_rsa', 'intra_chain_dist_weighted_rsa', 'Class', 'simpleClass',
                 'interacting_chains_str', 'interacting_ligands_str'
                 ]

    def __init__(self, res_num, aa='X', lig_dist_str=None, lig_dists=None, chain_dist_str=None, chain_distances=None, RSA=None,
                 relative_main_chain_acc=None, relative_side_chain_acc=None,
                 SSA=None, homo_dist_str=None, homomer_distances=None, interaction_profile=None, interaction_profile_str=None,
                 centrality_score_str=None, centralities=None, modres=None, b_factor=None, database_id=None, stored=False, phi=None, psi=None,
                 intra_ssbond=None, inter_ssbond=None, ssbond_length=None, intra_link=None, inter_link = None, link_length=None, cis_conformation=None, cis_follower=None,
                 inter_chain_median_kd=None, inter_chain_dist_weighted_kd=None, inter_chain_median_rsa=None,
                 inter_chain_dist_weighted_rsa=None, intra_chain_median_kd=None, intra_chain_dist_weighted_kd=None,
                 intra_chain_median_rsa=None, intra_chain_dist_weighted_rsa=None,
                 interacting_chains_str=None, interacting_ligands_str=None, inter_chain_interactions_median=None,
                 inter_chain_interactions_dist_weighted=None, intra_chain_interactions_median=None,
                 intra_chain_interactions_dist_weighted=None):

        self.res_num = res_num
        self.aa = aa
        self.lig_dist_str = lig_dist_str
        self.lig_dists = lig_dists
        self.chain_dist_str = chain_dist_str
        self.chain_distances = chain_distances
        self.RSA = RSA  # relative surface accessible area
        self.SSA = SSA  # secondary structure assignment
        self.relative_main_chain_acc = relative_main_chain_acc
        self.relative_side_chain_acc = relative_side_chain_acc
        self.homomer_distances = homomer_distances
        self.homo_dist_str = homo_dist_str
        self.interaction_profile = interaction_profile
        self.interaction_profile_str = interaction_profile_str  # interaction profile coded into a string, see rin.py for more information
        self.interacting_chains_str = interacting_chains_str
        self.interacting_ligands_str = interacting_ligands_str
        self.centralities = centralities
        self.centrality_score_str = centrality_score_str
        self.modres = modres
        self.b_factor = b_factor
        self.database_id = database_id
        self.stored = stored
        self.phi = phi
        self.psi = psi
        self.intra_ssbond = intra_ssbond
        self.inter_ssbond = inter_ssbond
        self.ssbond_length = ssbond_length
        self.intra_link = intra_link
        self.inter_link = inter_link
        self.link_length = link_length
        self.cis_conformation = cis_conformation
        self.cis_follower = cis_follower
        self.inter_chain_median_kd = inter_chain_median_kd
        self.inter_chain_dist_weighted_kd = inter_chain_dist_weighted_kd
        self.inter_chain_median_rsa = inter_chain_median_rsa
        self.inter_chain_dist_weighted_rsa = inter_chain_dist_weighted_rsa
        self.intra_chain_median_kd = intra_chain_median_kd
        self.intra_chain_dist_weighted_kd = intra_chain_dist_weighted_kd
        self.intra_chain_median_rsa = intra_chain_median_rsa
        self.intra_chain_dist_weighted_rsa = intra_chain_dist_weighted_rsa
        self.inter_chain_interactions_median = inter_chain_interactions_median
        self.inter_chain_interactions_dist_weighted = inter_chain_interactions_dist_weighted
        self.intra_chain_interactions_median = intra_chain_interactions_median
        self.intra_chain_interactions_dist_weighted = intra_chain_interactions_dist_weighted
        self.Class = None
        self.simpleClass = None

    def deconstruct(self):
        del self.res_num
        doomsday_protocol(self)

    def set_database_id(self, value):
        self.database_id = value

    def get_database_id(self):
        return self.database_id

    def get_aa(self):
        return self.aa

    def get_chain_distances(self):
        return self.chain_distances

    def get_chain_dist_str(self):
        if self.chain_dist_str is None:
            if self.chain_distances is None:
                return None
            self.convert_chain_dist_str()
        return self.chain_dist_str

    def convert_chain_dist_str(self):
        chain_dist_strings = []
        for chain_id in self.chain_distances:
            (dist, atom_pair, min_resi) = self.chain_distances[chain_id]
            if dist is not None:
                chain_dist_strings.append("%s.%s:%1.2f:(%s-%s)" % (chain_id, min_resi, dist, str(atom_pair[0]), str(atom_pair[1])))

        self.chain_dist_str = ",".join(chain_dist_strings)

    def get_ligand_distances(self):
        return self.lig_dists

    def get_lig_dist_str(self):
        if self.lig_dist_str is None:
            if self.lig_dists is None:
                return None
            self.convert_lig_dist_str()
        return self.lig_dist_str

    def convert_lig_dist_str(self):
        lig_dist_strings = []

        for lig_id in self.lig_dists:
            (dist, atom_pair) = self.lig_dists[lig_id]
            if dist is not None:
                lig_dist_strings.append("%s:%1.2f:(%s-%s)" % (lig_id, dist, str(atom_pair[0]), str(atom_pair[1])))

        self.lig_dist_str = ",".join(lig_dist_strings)

    def get_interaction_partners(self):
        if self.interaction_profile is None:
            self.get_interaction_profile()
        if self.interaction_profile is None:
            return None, None
        return (self.interaction_profile.interacting_chains, self.interaction_profile.interacting_ligands)

    def generate_interacting_chains_str(self):
        if self.interaction_profile is None:
            self.get_interaction_profile()
        if self.interaction_profile is None:
            return
        self.interacting_chains_str = ','.join(self.interaction_profile.interacting_chains)

    def generate_interacting_ligands_str(self):
        if self.interaction_profile is None:
            self.get_interaction_profile()
        if self.interaction_profile is None:
            return
        self.interacting_ligands_str = ','.join(self.interaction_profile.interacting_ligands)

    def get_shortest_distances(self, chains):
        min_ld = None
        min_md = None
        min_id = None
        min_cd = None
        min_dd = None
        min_rd = None
        min_hd = None

        ldists = {}
        mdists = {}
        idists = {}
        if self.lig_dists is not None:
            for lig_id in self.lig_dists:
                lig_name, res, chain = lig_id.split('_')
                (dist, atom_pair) = self.lig_dists[lig_id]
                if lig_name in METAL_ATOMS:
                    mdists[dist] = lig_name, res, chain
                elif lig_name in ION_ATOMS:
                    idists[dist] = lig_name, res, chain
                else:
                    ldists[dist] = lig_name, res, chain

        min_lig = None
        min_metal = None
        min_ion = None

        if len(ldists) > 0:
            min_ld = min(ldists.keys())
            min_lig = ldists[min_ld]

        if len(mdists) > 0:
            min_md = min(mdists.keys())
            min_metal = mdists[min_md]

        if len(idists) > 0:
            min_id = min(idists.keys())
            min_ion = idists[min_id]

        cdists = {}
        ddists = {}
        rdists = {}

        iacs = {}

        if self.chain_distances is not None:
            for chain_id in self.chain_distances:
                (dist, atom_pair, min_resi) = self.chain_distances[chain_id]
                if dist is None:
                    continue
                if chain_id not in chains:
                    chaintype = 'Protein'
                else:
                    chaintype = chains[chain_id]

                if chaintype == "Protein" or chaintype == 'Peptide':
                    cdists[dist] = chain_id
                elif chaintype == "RNA":
                    rdists[dist] = chain_id
                elif chaintype == "DNA":
                    ddists[dist] = chain_id

        if len(cdists) > 0:
            min_cd = min(cdists.keys())
            iacs['Protein'] = cdists[min_cd]
        if len(rdists) > 0:
            min_rd = min(rdists.keys())
            iacs['RNA'] = rdists[min_rd]
        if len(ddists) > 0:
            min_dd = min(ddists.keys())
            iacs['DNA'] = ddists[min_dd]

        homo_dists = []
        if self.homomer_distances is not None:
            for homo_chain in self.homomer_distances:
                dist = self.homomer_distances[homo_chain]
                homo_dists.append(dist)
        if len(homo_dists) > 0:
            min_hd = min(homo_dists)

        return min_hd, min_ld, min_md, min_id, min_cd, min_rd, min_dd, min_lig, min_metal, min_ion, iacs

    def get_homomer_dists(self):
        return self.homomer_distances

    def get_homo_dist_str(self):
        if self.homo_dist_str is None:
            if self.homomer_distances is None:
                return None
            self.convert_homo_dist_str()
        return self.homo_dist_str

    def convert_homo_dist_str(self):
        homo_strs = []
        for homo_chain in self.homomer_distances:
            min_d = self.homomer_distances[homo_chain]
            homo_strs.append('%s:%1.2f' % (homo_chain, min_d))
        self.homo_dist_str = ','.join(homo_strs)

    def get_centralities(self, get_whats_there=False):
        if get_whats_there:
            if self.centralities is not None:
                return self.centralities
            return self.centrality_score_str

        if self.centrality_score_str is not None and self.centralities is None:
            self.centralities = Centrality_scores(code_str=self.centrality_score_str)
        return self.centralities

    def get_centrality_str(self):
        if self.centrality_score_str is None:
            if self.centralities is None:
                return None
            else:
                self.convert_centrality_str()
        return self.centrality_score_str

    def convert_centrality_str(self):
        if self.centralities is None:
            return
        self.centrality_score_str = self.centralities.str_encode()

    def get_modres(self):
        return self.modres

    def get_b_factor(self):
        return self.b_factor

    def get_rsa(self, splitted=False):
        if not splitted:
            return self.RSA
        else:
            return (self.RSA, self.relative_main_chain_acc, self.relative_side_chain_acc)

    def get_ssa(self):
        return self.SSA

    def get_phi(self):
        return self.phi

    def get_psi(self):
        return self.psi

    def get_angles(self):
        return self.phi, self.psi

    def get_residue_link_information(self):
        return (self.intra_ssbond, self.inter_ssbond, self.ssbond_length, self.intra_link, self.inter_link, self.link_length, self.cis_conformation, self.cis_follower)

    def get_interaction_profile(self, get_whats_there=False):
        if get_whats_there:
            if self.interaction_profile is not None:
                return self.interaction_profile
            return self.interaction_profile_str

        if self.interaction_profile_str is not None and self.interaction_profile is None:
            self.interaction_profile = Interaction_profile(profile_str=self.interaction_profile_str,
                                                               interacting_chains_str=self.interacting_chains_str,
                                                               interacting_ligands_str=self.interacting_ligands_str)
        return self.interaction_profile

    def get_interacting_chains_str(self):
        if self.interacting_chains_str is None:
             self.generate_interacting_chains_str()
        return self.interacting_chains_str

    def get_interacting_ligands_str(self):
        if self.interacting_ligands_str is None:
             self.generate_interacting_ligands_str()
        return self.interacting_ligands_str

    def get_interaction_profile_str(self):
        if self.interaction_profile_str is None:
            if self.interaction_profile is None:
                return None
            else:
                self.convert_interaction_profile_str()
        return self.interaction_profile_str

    def convert_interaction_profile_str(self):
        if self.interaction_profile is None:
            return
        self.interaction_profile_str = self.interaction_profile.encode()

    def get_milieu(self):
        return (self.inter_chain_median_kd, self.inter_chain_dist_weighted_kd,
                self.inter_chain_median_rsa, self.inter_chain_dist_weighted_rsa, self.intra_chain_median_kd,
                self.intra_chain_dist_weighted_kd, self.intra_chain_median_rsa, self.intra_chain_dist_weighted_rsa)

    def get_interface_milieu(self):
        return (self.inter_chain_interactions_median,
                self.inter_chain_interactions_dist_weighted,
                self.intra_chain_interactions_median,
                self.intra_chain_interactions_dist_weighted)

    def get_classification(self, config):
        if self.Class is None and self.simpleClass is None:
            if self.RSA is None:
                sc = None
            else:
                if self.RSA > config.surface_threshold:
                    sc = "Surface"
                else:
                    sc = "Core"
            if self.interaction_profile is None:
                self.get_interaction_profile()
            if self.interaction_profile is None:
                config.errorlog.add_warning('Interaction profile is None, classification will fail')
            rin_class, rin_simple_class = rin_classify(self.interaction_profile, sc)
            self.Class = rin_class
            self.simpleClass = rin_simple_class
        return self.Class, self.simpleClass

    def set_classification(self, Class, simpleClass):
        self.Class = Class
        self.simpleClass = simpleClass
