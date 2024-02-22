# sdsc: structman datastructures and classes
import re
import numpy as np
from collections import Counter
from structman.lib import rin
from structman.lib.sdsc.sdsc_utils import rin_classify, triple_locate, doomsday_protocol, locate


MobiDB_map = {'D_PA': 'Polyampholite', 'D_WC': 'Weak polyampholie', 'D_NPE': 'D_NPE', 'D_PPE': 'D_PPE'}

def median(l):
    n = len(l)
    l = sorted(l)
    if n == 1:
        return l[0]
    if n == 0:
        return None
    if n % 2 == 0:
        med = (l[(n // 2) - 1] + l[n // 2]) / 2.0
    else:
        med = l[(n - 1) // 2]
    return med

class Mappings:
    __slots__ = [
                    'qualities', 'covs', 'seq_ids', 'rsas', 'mc_rsas', 'sc_rsas', 'ssas', 'lig_dists',
                    'metal_dists', 'ion_dists', 'chain_dists', 'rna_dists', 'dna_dists', 'homo_dists',
                    'profiles', 'weighted_profile', 'weighted_profile_str', 'centralities', 'weighted_centralities', 'weighted_centralities_str', 'b_factors',
                    'weighted_b_factor', 'modres', 'weighted_modres', 'weighted_ssa',
                    'phis', 'weighted_phi', 'psis', 'weighted_psi', 'intra_ssbonds', 'weighted_intra_ssbond', 'weighted_inter_ssbond', 'ssbond_lengths',
                    'weighted_ssbond_length', 'intra_links', 'weighted_intra_link', 'weighted_inter_link', 'link_lengths', 'weighted_link_length',
                    'cis_conformations', 'weighted_cis_conformation', 'cis_followers', 'weighted_cis_follower',
                    'inter_chain_median_kds', 'weighted_inter_chain_median_kd', 'inter_chain_dist_weighted_kds', 'weighted_inter_chain_dist_weighted_kd',
                    'inter_chain_median_rsas', 'weighted_inter_chain_median_rsa', 'inter_chain_dist_weighted_rsas',
                    'weighted_inter_chain_dist_weighted_rsa', 'intra_chain_median_kds', 'weighted_intra_chain_median_kd',
                    'intra_chain_dist_weighted_kds', 'weighted_intra_chain_dist_weighted_kd', 'intra_chain_median_rsas',
                    'weighted_intra_chain_median_rsa', 'intra_chain_dist_weighted_rsas', 'weighted_intra_chain_dist_weighted_rsa',
                    'inter_chain_interactions_medians', 'weighted_inter_chain_interactions_median',
                    'inter_chain_interactions_dist_weighteds', 'weighted_inter_chain_interactions_dist_weighted',
                    'intra_chain_interactions_medians', 'weighted_intra_chain_interactions_median',
                    'intra_chain_interactions_dist_weighteds', 'weighted_intra_chain_interactions_dist_weighted',
                    'weighted_lig_dist', 'lig_dist_conf', 'weighted_metal_dist', 'metal_dist_conf', 'weighted_ion_dist', 'ion_dist_conf',
                    'weighted_chain_dist', 'chain_dist_conf', 'weighted_rna_dist', 'rna_dist_conf', 'weighted_dna_dist', 'dna_dist_conf',
                    'weighted_homo_dist', 'homo_dist_conf', 'weighted_surface_value', 'weighted_mainchain_surface_value', 'weighted_sidechain_surface_value',
                    'weighted_location', 'weighted_mainchain_location', 'weighted_sidechain_location', 'max_cov', 'location_conf',
                    'res_classes', 'aa_ids', 'max_seq_res', 'recommended_res', 'interaction_recommendations', 'recommendation_order',
                    'rin_class', 'rin_simple_class', 'Class', 'simple_class', 'classification_conf', 'max_cov_rsa', 'resolutions', 'res_aas', 'amount_of_structures', 'rsa_change_score', 'mc_rsa_change_score', 'sc_rsa_change_score',
                    'Backbone_RMSD_seq_id_high_weight', 'All_atom_RMSD_seq_id_high_weight', 'nof_site_residue_seq_id_high_weight', 'Site_LDDT_seq_id_high_weight',
                    'Backbone_RMSD_seq_id_low_weight', 'All_atom_RMSD_seq_id_low_weight', 'nof_site_residue_seq_id_low_weight', 'Site_LDDT_seq_id_low_weight',
                    'Backbone_RMSD_seq_id_greater_90', 'All_atom_RMSD_seq_id_greater_90', 'nof_site_residue_seq_id_greater_90', 'Site_LDDT_seq_id_greater_90',
                    'Backbone_RMSD_seq_id_between_50_and_90', 'All_atom_RMSD_seq_id_between_50_and_90', 'nof_site_residue_seq_id_between_50_and_90', 'Site_LDDT_seq_id_between_50_and_90',
                    'Backbone_RMSD_seq_id_lower_50', 'All_atom_RMSD_seq_id_lower_50', 'nof_site_residue_seq_id_lower_50', 'Site_LDDT_seq_id_lower_50',
                    'Backbone_RMSD_site_id_greater_99', 'All_atom_RMSD_site_id_greater_99', 'nof_site_residue_site_id_greater_99', 'Site_LDDT_site_id_greater_99',
                    'Backbone_RMSD_site_id_between_70_and_99', 'All_atom_RMSD_site_id_between_70_and_99', 'nof_site_residue_site_id_between_70_and_99', 'Site_LDDT_site_id_between_70_and_99',
                    'Backbone_RMSD_site_id_lower_70', 'All_atom_RMSD_site_id_lower_70', 'nof_site_residue_site_id_lower_70', 'Site_LDDT_site_id_lower_70',
                    'Backbone_RMSD_seq_id_greater_90_site_id_greater_99', 'All_atom_RMSD_seq_id_greater_90_site_id_greater_99', 'nof_site_residue_seq_id_greater_90_site_id_greater_99', 'Site_LDDT_seq_id_greater_90_site_id_greater_99',
                    'Backbone_RMSD_seq_id_between_50_and_90_site_id_greater_99', 'All_atom_RMSD_seq_id_between_50_and_90_site_id_greater_99', 'nof_site_residue_seq_id_between_50_and_90_site_id_greater_99', 'Site_LDDT_seq_id_between_50_and_90_site_id_greater_99',
                    'Backbone_RMSD_seq_id_lower_50_site_id_greater_99', 'All_atom_RMSD_seq_id_lower_50_site_id_greater_99', 'nof_site_residue_seq_id_lower_50_site_id_greater_99', 'Site_LDDT_seq_id_lower_50_site_id_greater_99',
                    'Backbone_RMSD_seq_id_greater_90_site_id_between_70_and_99', 'All_atom_RMSD_seq_id_greater_90_site_id_between_70_and_99', 'nof_site_residue_seq_id_greater_90_site_id_between_70_and_99', 'Site_LDDT_seq_id_greater_90_site_id_between_70_and_99',
                    'Backbone_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99', 'All_atom_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99', 'nof_site_residue_seq_id_between_50_and_90_site_id_between_70_and_99', 'Site_LDDT_seq_id_between_50_and_90_site_id_between_70_and_99',
                    'Backbone_RMSD_seq_id_lower_50_site_id_between_70_and_99', 'All_atom_RMSD_seq_id_lower_50_site_id_between_70_and_99', 'nof_site_residue_seq_id_lower_50_site_id_between_70_and_99', 'Site_LDDT_seq_id_lower_50_site_id_between_70_and_99',
                    'Backbone_RMSD_seq_id_greater_90_site_id_lower_70', 'All_atom_RMSD_seq_id_greater_90_site_id_lower_70', 'nof_site_residue_seq_id_greater_90_site_id_lower_70', 'Site_LDDT_seq_id_greater_90_site_id_lower_70',
                    'Backbone_RMSD_seq_id_between_50_and_90_site_id_lower_70', 'All_atom_RMSD_seq_id_between_50_and_90_site_id_lower_70', 'nof_site_residue_seq_id_between_50_and_90_site_id_lower_70', 'Site_LDDT_seq_id_between_50_and_90_site_id_lower_70',
                    'Backbone_RMSD_seq_id_lower_50_site_id_lower_70', 'All_atom_RMSD_seq_id_lower_50_site_id_lower_70', 'nof_site_residue_seq_id_lower_50_site_id_lower_70', 'Site_LDDT_seq_id_lower_50_site_id_lower_70'
                ]

    def __init__(self, raw_results=None):
        self.qualities = {}
        self.seq_ids = {}
        self.covs = {}
        self.rsas = {}
        self.mc_rsas = {}
        self.sc_rsas = {}
        self.ssas = {}
        self.lig_dists = {}
        self.metal_dists = {}
        self.ion_dists = {}
        self.chain_dists = {}
        self.rna_dists = {}
        self.dna_dists = {}
        self.homo_dists = {}
        self.profiles = {}
        self.centralities = {}
        self.phis = {}
        self.psis = {}
        self.intra_ssbonds = {}
        self.ssbond_lengths = {}
        self.intra_links = {}
        self.link_lengths = {}
        self.cis_conformations = {}
        self.cis_followers = {}
        self.inter_chain_median_kds = {}
        self.inter_chain_dist_weighted_kds = {}
        self.inter_chain_median_rsas = {}
        self.inter_chain_dist_weighted_rsas = {}
        self.intra_chain_median_kds = {}
        self.intra_chain_dist_weighted_kds = {}
        self.intra_chain_median_rsas = {}
        self.intra_chain_dist_weighted_rsas = {}
        self.inter_chain_interactions_medians = {}
        self.inter_chain_interactions_dist_weighteds = {}
        self.intra_chain_interactions_medians = {}
        self.intra_chain_interactions_dist_weighteds = {}
        self.b_factors = {}
        self.modres = {}
        self.res_classes = {}
        self.aa_ids = {}
        self.res_aas = {}
        self.resolutions = {}
        self.recommended_res = None
        self.recommendation_order = []
        self.max_seq_res = None
        self.weighted_surface_value = None
        self.weighted_mainchain_surface_value = None
        self.weighted_sidechain_surface_value = None
        self.weighted_location = None
        self.weighted_mainchain_location = None
        self.weighted_sidechain_location = None
        self.weighted_profile = None
        self.weighted_profile_str = None
        self.weighted_centralities = None
        self.weighted_centralities_str = None
        self.weighted_b_factor = None
        self.weighted_modres = None
        self.weighted_ssa = None
        self.weighted_phi = None
        self.weighted_psi = None
        self.weighted_intra_ssbond = None
        self.weighted_inter_ssbond = None
        self.weighted_ssbond_length = None
        self.weighted_intra_link = None
        self.weighted_inter_link = None
        self.weighted_link_length = None
        self.weighted_cis_conformation = None
        self.weighted_cis_follower = None
        self.weighted_inter_chain_median_kd = None
        self.weighted_inter_chain_dist_weighted_kd = None
        self.weighted_inter_chain_median_rsa = None
        self.weighted_inter_chain_dist_weighted_rsa = None
        self.weighted_intra_chain_median_kd = None
        self.weighted_intra_chain_dist_weighted_kd = None
        self.weighted_intra_chain_median_rsa = None
        self.weighted_intra_chain_dist_weighted_rsa = None
        self.weighted_inter_chain_interactions_median = None
        self.weighted_inter_chain_interactions_dist_weighted = None
        self.weighted_intra_chain_interactions_median = None
        self.weighted_intra_chain_interactions_dist_weighted = None
        self.weighted_lig_dist = None
        self.weighted_metal_dist = None
        self.weighted_ion_dist = None
        self.weighted_chain_dist = None
        self.weighted_rna_dist = None
        self.weighted_dna_dist = None
        self.weighted_homo_dist = None
        self.rin_class = None
        self.rin_simple_class = None
        self.Class = None
        self.simple_class = None
        self.interaction_recommendations = None
        self.lig_dist_conf = None
        self.metal_dist_conf = None
        self.ion_dist_conf = None
        self.chain_dist_conf = None
        self.rna_dist_conf = None
        self.dna_dist_conf = None
        self.homo_dist_conf = None
        self.location_conf = None
        self.classification_conf = None
        self.amount_of_structures = 0
        self.rsa_change_score = None
        self.mc_rsa_change_score = None
        self.sc_rsa_change_score = None

        self.Backbone_RMSD_seq_id_high_weight = None
        self.All_atom_RMSD_seq_id_high_weight = None
        self.nof_site_residue_seq_id_high_weight = None
        self.Site_LDDT_seq_id_high_weight = None
        self.Backbone_RMSD_seq_id_low_weight = None
        self.All_atom_RMSD_seq_id_low_weight = None
        self.nof_site_residue_seq_id_low_weight = None
        self.Site_LDDT_seq_id_low_weight = None
        self.Backbone_RMSD_seq_id_greater_90 = None
        self.All_atom_RMSD_seq_id_greater_90 = None
        self.nof_site_residue_seq_id_greater_90 = None
        self.Site_LDDT_seq_id_greater_90 = None
        self.Backbone_RMSD_seq_id_between_50_and_90 = None
        self.All_atom_RMSD_seq_id_between_50_and_90 = None
        self.nof_site_residue_seq_id_between_50_and_90 = None
        self.Site_LDDT_seq_id_between_50_and_90 = None
        self.Backbone_RMSD_seq_id_lower_50 = None
        self.All_atom_RMSD_seq_id_lower_50 = None
        self.nof_site_residue_seq_id_lower_50 = None
        self.Site_LDDT_seq_id_lower_50 = None
        self.Backbone_RMSD_site_id_greater_99 = None
        self.All_atom_RMSD_site_id_greater_99 = None
        self.nof_site_residue_site_id_greater_99 = None
        self.Site_LDDT_site_id_greater_99 = None
        self.Backbone_RMSD_site_id_between_70_and_99 = None
        self.All_atom_RMSD_site_id_between_70_and_99 = None
        self.nof_site_residue_site_id_between_70_and_99 = None
        self.Site_LDDT_site_id_between_70_and_99 = None
        self.Backbone_RMSD_site_id_lower_70 = None
        self.All_atom_RMSD_site_id_lower_70 = None
        self.nof_site_residue_site_id_lower_70 = None
        self.Site_LDDT_site_id_lower_70 = None
        self.Backbone_RMSD_seq_id_greater_90_site_id_greater_99 = None
        self.All_atom_RMSD_seq_id_greater_90_site_id_greater_99 = None
        self.nof_site_residue_seq_id_greater_90_site_id_greater_99 = None
        self.Site_LDDT_seq_id_greater_90_site_id_greater_99 = None
        self.Backbone_RMSD_seq_id_between_50_and_90_site_id_greater_99 = None
        self.All_atom_RMSD_seq_id_between_50_and_90_site_id_greater_99 = None
        self.nof_site_residue_seq_id_between_50_and_90_site_id_greater_99 = None
        self.Site_LDDT_seq_id_between_50_and_90_site_id_greater_99 = None
        self.Backbone_RMSD_seq_id_lower_50_site_id_greater_99 = None
        self.All_atom_RMSD_seq_id_lower_50_site_id_greater_99 = None
        self.nof_site_residue_seq_id_lower_50_site_id_greater_99 = None
        self.Site_LDDT_seq_id_lower_50_site_id_greater_99 = None
        self.Backbone_RMSD_seq_id_greater_90_site_id_between_70_and_99 = None
        self.All_atom_RMSD_seq_id_greater_90_site_id_between_70_and_99 = None
        self.nof_site_residue_seq_id_greater_90_site_id_between_70_and_99 = None
        self.Site_LDDT_seq_id_greater_90_site_id_between_70_and_99 = None
        self.Backbone_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99 = None
        self.All_atom_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99 = None
        self.nof_site_residue_seq_id_between_50_and_90_site_id_between_70_and_99 = None
        self.Site_LDDT_seq_id_between_50_and_90_site_id_between_70_and_99 = None
        self.Backbone_RMSD_seq_id_lower_50_site_id_between_70_and_99 = None
        self.All_atom_RMSD_seq_id_lower_50_site_id_between_70_and_99 = None
        self.nof_site_residue_seq_id_lower_50_site_id_between_70_and_99 = None
        self.Site_LDDT_seq_id_lower_50_site_id_between_70_and_99 = None
        self.Backbone_RMSD_seq_id_greater_90_site_id_lower_70 = None
        self.All_atom_RMSD_seq_id_greater_90_site_id_lower_70 = None
        self.nof_site_residue_seq_id_greater_90_site_id_lower_70 = None
        self.Site_LDDT_seq_id_greater_90_site_id_lower_70 = None
        self.Backbone_RMSD_seq_id_between_50_and_90_site_id_lower_70 = None
        self.All_atom_RMSD_seq_id_between_50_and_90_site_id_lower_70 = None
        self.nof_site_residue_seq_id_between_50_and_90_site_id_lower_70 = None
        self.Site_LDDT_seq_id_between_50_and_90_site_id_lower_70 = None
        self.Backbone_RMSD_seq_id_lower_50_site_id_lower_70 = None
        self.All_atom_RMSD_seq_id_lower_50_site_id_lower_70 = None
        self.nof_site_residue_seq_id_lower_50_site_id_lower_70 = None
        self.Site_LDDT_seq_id_lower_50_site_id_lower_70 = None

        if raw_results is not None:
            (self.recommendation_order, self.recommended_res, self.max_seq_res, self.weighted_surface_value, self.weighted_mainchain_surface_value,
            self.weighted_sidechain_surface_value, self.weighted_location, self.weighted_mainchain_location,
            self.weighted_sidechain_location, self.weighted_profile_str,
            self.weighted_centralities_str, self.weighted_b_factor, self.weighted_modres, self.weighted_ssa, self.weighted_phi,
            self.weighted_psi, self.weighted_intra_ssbond, self.weighted_inter_ssbond, self.weighted_ssbond_length, self.weighted_intra_link,
            self.weighted_inter_link, self.weighted_link_length, self.weighted_cis_conformation, self.weighted_cis_follower,
            self.weighted_inter_chain_median_kd, self.weighted_inter_chain_dist_weighted_kd, self.weighted_inter_chain_median_rsa,
            self.weighted_inter_chain_dist_weighted_rsa, self.weighted_intra_chain_median_kd, self.weighted_intra_chain_dist_weighted_kd,
            self.weighted_intra_chain_median_rsa, self.weighted_intra_chain_dist_weighted_rsa,
            self.weighted_inter_chain_interactions_median, self.weighted_inter_chain_interactions_dist_weighted,
            self.weighted_intra_chain_interactions_median, self.weighted_intra_chain_interactions_dist_weighted,
            self.weighted_lig_dist,
            self.weighted_metal_dist, self.weighted_ion_dist, self.weighted_chain_dist, self.weighted_rna_dist, self.weighted_dna_dist,
            self.weighted_homo_dist, self.rin_class, self.rin_simple_class, self.Class, self.simple_class, self.interaction_recommendations,
            self.lig_dist_conf, self.metal_dist_conf, self.ion_dist_conf, self.chain_dist_conf, self.rna_dist_conf, self.dna_dist_conf,
            self.homo_dist_conf, self.location_conf, self.classification_conf, self.amount_of_structures, self.rsa_change_score, self.mc_rsa_change_score, self.sc_rsa_change_score,
            self.Backbone_RMSD_seq_id_high_weight, self.All_atom_RMSD_seq_id_high_weight, self.nof_site_residue_seq_id_high_weight, self.Site_LDDT_seq_id_high_weight,
            self.Backbone_RMSD_seq_id_low_weight, self.All_atom_RMSD_seq_id_low_weight, self.nof_site_residue_seq_id_low_weight, self.Site_LDDT_seq_id_low_weight,
            self.Backbone_RMSD_seq_id_greater_90, self.All_atom_RMSD_seq_id_greater_90, self.nof_site_residue_seq_id_greater_90, self.Site_LDDT_seq_id_greater_90,
            self.Backbone_RMSD_seq_id_between_50_and_90, self.All_atom_RMSD_seq_id_between_50_and_90, self.nof_site_residue_seq_id_between_50_and_90, self.Site_LDDT_seq_id_between_50_and_90,
            self.Backbone_RMSD_seq_id_lower_50, self.All_atom_RMSD_seq_id_lower_50, self.nof_site_residue_seq_id_lower_50, self.Site_LDDT_seq_id_lower_50,
            self.Backbone_RMSD_site_id_greater_99, self.All_atom_RMSD_site_id_greater_99, self.nof_site_residue_site_id_greater_99, self.Site_LDDT_site_id_greater_99,
            self.Backbone_RMSD_site_id_between_70_and_99, self.All_atom_RMSD_site_id_between_70_and_99, self.nof_site_residue_site_id_between_70_and_99, self.Site_LDDT_site_id_between_70_and_99,
            self.Backbone_RMSD_site_id_lower_70, self.All_atom_RMSD_site_id_lower_70, self.nof_site_residue_site_id_lower_70, self.Site_LDDT_site_id_lower_70,
            self.Backbone_RMSD_seq_id_greater_90_site_id_greater_99, self.All_atom_RMSD_seq_id_greater_90_site_id_greater_99, self.nof_site_residue_seq_id_greater_90_site_id_greater_99, self.Site_LDDT_seq_id_greater_90_site_id_greater_99,
            self.Backbone_RMSD_seq_id_between_50_and_90_site_id_greater_99, self.All_atom_RMSD_seq_id_between_50_and_90_site_id_greater_99, self.nof_site_residue_seq_id_between_50_and_90_site_id_greater_99, self.Site_LDDT_seq_id_between_50_and_90_site_id_greater_99,
            self.Backbone_RMSD_seq_id_lower_50_site_id_greater_99, self.All_atom_RMSD_seq_id_lower_50_site_id_greater_99, self.nof_site_residue_seq_id_lower_50_site_id_greater_99, self.Site_LDDT_seq_id_lower_50_site_id_greater_99,
            self.Backbone_RMSD_seq_id_greater_90_site_id_between_70_and_99, self.All_atom_RMSD_seq_id_greater_90_site_id_between_70_and_99, self.nof_site_residue_seq_id_greater_90_site_id_between_70_and_99, self.Site_LDDT_seq_id_greater_90_site_id_between_70_and_99,
            self.Backbone_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99, self.All_atom_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99, self.nof_site_residue_seq_id_between_50_and_90_site_id_between_70_and_99, self.Site_LDDT_seq_id_between_50_and_90_site_id_between_70_and_99,
            self.Backbone_RMSD_seq_id_lower_50_site_id_between_70_and_99, self.All_atom_RMSD_seq_id_lower_50_site_id_between_70_and_99, self.nof_site_residue_seq_id_lower_50_site_id_between_70_and_99, self.Site_LDDT_seq_id_lower_50_site_id_between_70_and_99,
            self.Backbone_RMSD_seq_id_greater_90_site_id_lower_70, self.All_atom_RMSD_seq_id_greater_90_site_id_lower_70, self.nof_site_residue_seq_id_greater_90_site_id_lower_70, self.Site_LDDT_seq_id_greater_90_site_id_lower_70,
            self.Backbone_RMSD_seq_id_between_50_and_90_site_id_lower_70, self.All_atom_RMSD_seq_id_between_50_and_90_site_id_lower_70, self.nof_site_residue_seq_id_between_50_and_90_site_id_lower_70, self.Site_LDDT_seq_id_between_50_and_90_site_id_lower_70,
            self.Backbone_RMSD_seq_id_lower_50_site_id_lower_70, self.All_atom_RMSD_seq_id_lower_50_site_id_lower_70, self.nof_site_residue_seq_id_lower_50_site_id_lower_70, self.Site_LDDT_seq_id_lower_50_site_id_lower_70) = raw_results

    def deconstruct(self):
        doomsday_protocol(self)

    def add_mapping(self, mapping_id, mapping):
        (quality, seq_id, cov, rsa, mc_rsa, sc_rsa, ssa, lig_dist, metal_dist, ion_dist, chain_dist, rna_dist, dna_dist, homo_dist, profile, centralities,
         phi, psi, intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
         inter_chain_median_kd, inter_chain_dist_weighted_kd,
         inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
         intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
         inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
         intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
         b_factor, modres, res_class, res_simple_class,
         identical_aa, resolution, res_aa) = mapping
        if mapping_id not in self.qualities:
            self.amount_of_structures += 1
        self.qualities[mapping_id] = quality
        self.seq_ids[mapping_id] = seq_id
        self.covs[mapping_id] = cov
        self.rsas[mapping_id] = rsa
        self.mc_rsas[mapping_id] = mc_rsa
        self.sc_rsas[mapping_id] = sc_rsa
        self.ssas[mapping_id] = ssa
        self.lig_dists[mapping_id] = lig_dist
        self.metal_dists[mapping_id] = metal_dist
        self.ion_dists[mapping_id] = ion_dist
        self.chain_dists[mapping_id] = chain_dist
        self.rna_dists[mapping_id] = rna_dist
        self.dna_dists[mapping_id] = dna_dist
        self.homo_dists[mapping_id] = homo_dist
        self.profiles[mapping_id] = profile
        self.centralities[mapping_id] = centralities
        self.phis[mapping_id] = phi
        self.psis[mapping_id] = psi
        self.intra_ssbonds[mapping_id] = intra_ssbond
        self.ssbond_lengths[mapping_id] = ssbond_length
        self.intra_links[mapping_id] = intra_link
        self.link_lengths[mapping_id] = link_length
        self.cis_conformations[mapping_id] = cis_conformation
        self.cis_followers[mapping_id] = cis_follower
        self.inter_chain_median_kds[mapping_id] = inter_chain_median_kd
        self.inter_chain_dist_weighted_kds[mapping_id] = inter_chain_dist_weighted_kd
        self.inter_chain_median_rsas[mapping_id] = inter_chain_median_rsa
        self.inter_chain_dist_weighted_rsas[mapping_id] = inter_chain_dist_weighted_rsa
        self.intra_chain_median_kds[mapping_id] = intra_chain_median_kd
        self.intra_chain_dist_weighted_kds[mapping_id] = intra_chain_dist_weighted_kd
        self.intra_chain_median_rsas[mapping_id] = intra_chain_median_rsa
        self.intra_chain_dist_weighted_rsas[mapping_id] = intra_chain_dist_weighted_rsa
        self.inter_chain_interactions_medians[mapping_id] = inter_chain_interactions_median
        self.inter_chain_interactions_dist_weighteds[mapping_id] = inter_chain_interactions_dist_weighted
        self.intra_chain_interactions_medians[mapping_id] = intra_chain_interactions_median
        self.intra_chain_interactions_dist_weighteds[mapping_id] = intra_chain_interactions_dist_weighted
        self.b_factors[mapping_id] = b_factor
        self.modres[mapping_id] = modres
        self.res_classes[mapping_id] = res_class, res_simple_class
        self.aa_ids[mapping_id] = identical_aa
        self.res_aas[mapping_id] = res_aa
        self.resolutions[mapping_id[0]] = resolution

    def add_result(self, mapping_id, raw_results, quality, seq_id, cov):

        (rsa, mc_rsa, sc_rsa, profile, centralities, b_factor, modres, ssa, phi,
             psi, intra_ssbond, ssbond_length, intra_link,
             link_length, cis_conformation, cis_follower,
             inter_chain_median_kd, inter_chain_dist_weighted_kd, inter_chain_median_rsa,
             inter_chain_dist_weighted_rsa, intra_chain_median_kd, intra_chain_dist_weighted_kd,
             intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
             inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
             intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
             lig_dist, metal_dist, ion_dist, chain_dist, rna_dist, dna_dist,
             homo_dist, res_class, res_simple_class) = raw_results

        if mapping_id not in self.qualities:
            self.amount_of_structures += 1
        self.qualities[mapping_id] = quality
        self.seq_ids[mapping_id] = seq_id
        self.covs[mapping_id] = cov
        self.rsas[mapping_id] = rsa
        self.mc_rsas[mapping_id] = mc_rsa
        self.sc_rsas[mapping_id] = sc_rsa
        self.ssas[mapping_id] = ssa
        self.lig_dists[mapping_id] = lig_dist
        self.metal_dists[mapping_id] = metal_dist
        self.ion_dists[mapping_id] = ion_dist
        self.chain_dists[mapping_id] = chain_dist
        self.rna_dists[mapping_id] = rna_dist
        self.dna_dists[mapping_id] = dna_dist
        self.homo_dists[mapping_id] = homo_dist
        self.profiles[mapping_id] = profile
        self.centralities[mapping_id] = centralities
        self.phis[mapping_id] = phi
        self.psis[mapping_id] = psi
        self.intra_ssbonds[mapping_id] = intra_ssbond
        self.ssbond_lengths[mapping_id] = ssbond_length
        self.intra_links[mapping_id] = intra_link
        self.link_lengths[mapping_id] = link_length
        self.cis_conformations[mapping_id] = cis_conformation
        self.cis_followers[mapping_id] = cis_follower
        self.inter_chain_median_kds[mapping_id] = inter_chain_median_kd
        self.inter_chain_dist_weighted_kds[mapping_id] = inter_chain_dist_weighted_kd
        self.inter_chain_median_rsas[mapping_id] = inter_chain_median_rsa
        self.inter_chain_dist_weighted_rsas[mapping_id] = inter_chain_dist_weighted_rsa
        self.intra_chain_median_kds[mapping_id] = intra_chain_median_kd
        self.intra_chain_dist_weighted_kds[mapping_id] = intra_chain_dist_weighted_kd
        self.intra_chain_median_rsas[mapping_id] = intra_chain_median_rsa
        self.intra_chain_dist_weighted_rsas[mapping_id] = intra_chain_dist_weighted_rsa
        self.inter_chain_interactions_medians[mapping_id] = inter_chain_interactions_median
        self.inter_chain_interactions_dist_weighteds[mapping_id] = inter_chain_interactions_dist_weighted
        self.intra_chain_interactions_medians[mapping_id] = intra_chain_interactions_median
        self.intra_chain_interactions_dist_weighteds[mapping_id] = intra_chain_interactions_dist_weighted
        self.b_factors[mapping_id] = b_factor
        self.modres[mapping_id] = modres
        self.res_classes[mapping_id] = res_class, res_simple_class
        self.aa_ids[mapping_id] = True
        self.res_aas[mapping_id] = None

    def set_weighted_lig_dist(self, wd, conf):
        self.weighted_lig_dist = wd
        self.lig_dist_conf = conf

    def set_weighted_metal_dist(self, wd, conf):
        self.weighted_metal_dist = wd
        self.metal_dist_conf = conf

    def set_weighted_ion_dist(self, wd, conf):
        self.weighted_ion_dist = wd
        self.ion_dist_conf = conf

    def set_weighted_chain_dist(self, wd, conf):
        self.weighted_chain_dist = wd
        self.chain_dist_conf = conf

    def set_weighted_rna_dist(self, wd, conf):
        self.weighted_rna_dist = wd
        self.rna_dist_conf = conf

    def set_weighted_dna_dist(self, wd, conf):
        self.weighted_dna_dist = wd
        self.dna_dist_conf = conf

    def set_weighted_homo_dist(self, wd, conf):
        self.weighted_homo_dist = wd
        self.homo_dist_conf = conf

    def set_weighted_phi(self, value):
        self.weighted_phi = value

    def set_weighted_psi(self, value):
        self.weighted_psi = value

    def set_weighted_intra_ssbond(self, value):
        self.weighted_intra_ssbond = value

    def set_weighted_inter_ssbond(self, value):
        self.weighted_inter_ssbond = value

    def set_weighted_ssbond_length(self, value):
        self.weighted_ssbond_length = value

    def set_weighted_intra_link(self, value):
        self.weighted_intra_link = value

    def set_weighted_inter_link(self, value):
        self.weighted_inter_link = value

    def set_weighted_link_length(self, value):
        self.weighted_link_length = value

    def set_weighted_cis_conformation(self, value):
        self.weighted_cis_conformation = value

    def set_weighted_cis_follower(self, value):
        self.weighted_cis_follower = value

    def set_weighted_inter_chain_median_kd(self, value):
        self.weighted_inter_chain_median_kd = value

    def set_weighted_inter_chain_dist_weighted_kd(self, value):
        self.weighted_inter_chain_dist_weighted_kd = value

    def set_weighted_inter_chain_median_rsa(self, value):
        self.weighted_inter_chain_median_rsa = value

    def set_weighted_inter_chain_dist_weighted_rsa(self, value):
        self.weighted_inter_chain_dist_weighted_rsa = value

    def set_weighted_intra_chain_median_kd(self, value):
        self.weighted_intra_chain_median_kd = value

    def set_weighted_intra_chain_dist_weighted_kd(self, value):
        self.weighted_intra_chain_dist_weighted_kd = value

    def set_weighted_intra_chain_median_rsa(self, value):
        self.weighted_intra_chain_median_rsa = value

    def set_weighted_intra_chain_dist_weighted_rsa(self, value):
        self.weighted_intra_chain_dist_weighted_rsa = value

    def set_weighted_inter_chain_interactions_median(self, value):
        self.weighted_inter_chain_interactions_median = value

    def set_weighted_inter_chain_interactions_dist_weighted(self, value):
        self.weighted_inter_chain_interactions_dist_weighted = value

    def set_weighted_intra_chain_interactions_median(self, value):
        self.weighted_intra_chain_interactions_median = value

    def set_weighted_intra_chain_interactions_dist_weighted(self, value):
        self.weighted_intra_chain_interactions_dist_weighted = value

    def set_weighted_modres(self, value):
        self.weighted_modres = value

    def set_weighted_b_factor(self, value):
        self.weighted_b_factor = value

    def get_raw_result(self):
        return (self.recommendation_order, self.recommended_res, self.max_seq_res, self.weighted_surface_value, self.weighted_mainchain_surface_value,
                self.weighted_sidechain_surface_value, self.weighted_location, self.weighted_mainchain_location,
                self.weighted_sidechain_location, self.get_weighted_profile_str(),
                self.get_weighted_centralities_str(), self.weighted_b_factor, self.weighted_modres, self.weighted_ssa, self.weighted_phi,
                self.weighted_psi, self.weighted_intra_ssbond, self.weighted_inter_ssbond, self.weighted_ssbond_length, self.weighted_intra_link,
                self.weighted_inter_link, self.weighted_link_length, self.weighted_cis_conformation, self.weighted_cis_follower,
                self.weighted_inter_chain_median_kd, self.weighted_inter_chain_dist_weighted_kd, self.weighted_inter_chain_median_rsa,
                self.weighted_inter_chain_dist_weighted_rsa, self.weighted_intra_chain_median_kd, self.weighted_intra_chain_dist_weighted_kd,
                self.weighted_intra_chain_median_rsa, self.weighted_intra_chain_dist_weighted_rsa,
                self.weighted_inter_chain_interactions_median, self.weighted_inter_chain_interactions_dist_weighted,
                self.weighted_intra_chain_interactions_median, self.weighted_intra_chain_interactions_dist_weighted,
                self.weighted_lig_dist,
                self.weighted_metal_dist, self.weighted_ion_dist, self.weighted_chain_dist, self.weighted_rna_dist, self.weighted_dna_dist,
                self.weighted_homo_dist, self.rin_class, self.rin_simple_class, self.Class, self.simple_class, self.interaction_recommendations,
                self.lig_dist_conf, self.metal_dist_conf, self.ion_dist_conf, self.chain_dist_conf, self.rna_dist_conf, self.dna_dist_conf,
                self.homo_dist_conf, self.location_conf, self.classification_conf, self.amount_of_structures, self.rsa_change_score, self.mc_rsa_change_score, self.sc_rsa_change_score,
                self.Backbone_RMSD_seq_id_high_weight, self.All_atom_RMSD_seq_id_high_weight, self.nof_site_residue_seq_id_high_weight, self.Site_LDDT_seq_id_high_weight,
                self.Backbone_RMSD_seq_id_low_weight, self.All_atom_RMSD_seq_id_low_weight, self.nof_site_residue_seq_id_low_weight, self.Site_LDDT_seq_id_low_weight,
                self.Backbone_RMSD_seq_id_greater_90, self.All_atom_RMSD_seq_id_greater_90, self.nof_site_residue_seq_id_greater_90, self.Site_LDDT_seq_id_greater_90,
                self.Backbone_RMSD_seq_id_between_50_and_90, self.All_atom_RMSD_seq_id_between_50_and_90, self.nof_site_residue_seq_id_between_50_and_90, self.Site_LDDT_seq_id_between_50_and_90,
                self.Backbone_RMSD_seq_id_lower_50, self.All_atom_RMSD_seq_id_lower_50, self.nof_site_residue_seq_id_lower_50, self.Site_LDDT_seq_id_lower_50,
                self.Backbone_RMSD_site_id_greater_99, self.All_atom_RMSD_site_id_greater_99, self.nof_site_residue_site_id_greater_99, self.Site_LDDT_site_id_greater_99,
                self.Backbone_RMSD_site_id_between_70_and_99, self.All_atom_RMSD_site_id_between_70_and_99, self.nof_site_residue_site_id_between_70_and_99, self.Site_LDDT_site_id_between_70_and_99,
                self.Backbone_RMSD_site_id_lower_70, self.All_atom_RMSD_site_id_lower_70, self.nof_site_residue_site_id_lower_70, self.Site_LDDT_site_id_lower_70,
                self.Backbone_RMSD_seq_id_greater_90_site_id_greater_99, self.All_atom_RMSD_seq_id_greater_90_site_id_greater_99, self.nof_site_residue_seq_id_greater_90_site_id_greater_99, self.Site_LDDT_seq_id_greater_90_site_id_greater_99,
                self.Backbone_RMSD_seq_id_between_50_and_90_site_id_greater_99, self.All_atom_RMSD_seq_id_between_50_and_90_site_id_greater_99, self.nof_site_residue_seq_id_between_50_and_90_site_id_greater_99, self.Site_LDDT_seq_id_between_50_and_90_site_id_greater_99,
                self.Backbone_RMSD_seq_id_lower_50_site_id_greater_99, self.All_atom_RMSD_seq_id_lower_50_site_id_greater_99, self.nof_site_residue_seq_id_lower_50_site_id_greater_99, self.Site_LDDT_seq_id_lower_50_site_id_greater_99,
                self.Backbone_RMSD_seq_id_greater_90_site_id_between_70_and_99, self.All_atom_RMSD_seq_id_greater_90_site_id_between_70_and_99, self.nof_site_residue_seq_id_greater_90_site_id_between_70_and_99, self.Site_LDDT_seq_id_greater_90_site_id_between_70_and_99,
                self.Backbone_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99, self.All_atom_RMSD_seq_id_between_50_and_90_site_id_between_70_and_99, self.nof_site_residue_seq_id_between_50_and_90_site_id_between_70_and_99, self.Site_LDDT_seq_id_between_50_and_90_site_id_between_70_and_99,
                self.Backbone_RMSD_seq_id_lower_50_site_id_between_70_and_99, self.All_atom_RMSD_seq_id_lower_50_site_id_between_70_and_99, self.nof_site_residue_seq_id_lower_50_site_id_between_70_and_99, self.Site_LDDT_seq_id_lower_50_site_id_between_70_and_99,
                self.Backbone_RMSD_seq_id_greater_90_site_id_lower_70, self.All_atom_RMSD_seq_id_greater_90_site_id_lower_70, self.nof_site_residue_seq_id_greater_90_site_id_lower_70, self.Site_LDDT_seq_id_greater_90_site_id_lower_70,
                self.Backbone_RMSD_seq_id_between_50_and_90_site_id_lower_70, self.All_atom_RMSD_seq_id_between_50_and_90_site_id_lower_70, self.nof_site_residue_seq_id_between_50_and_90_site_id_lower_70, self.Site_LDDT_seq_id_between_50_and_90_site_id_lower_70,
                self.Backbone_RMSD_seq_id_lower_50_site_id_lower_70, self.All_atom_RMSD_seq_id_lower_50_site_id_lower_70, self.nof_site_residue_seq_id_lower_50_site_id_lower_70, self.Site_LDDT_seq_id_lower_50_site_id_lower_70)

    def get_result_for_indel_aggregation(self):
        
        return (self.weighted_surface_value, self.weighted_mainchain_surface_value, self.weighted_sidechain_surface_value,
             self.get_weighted_profile(), self.get_weighted_centralities(), self.weighted_b_factor, self.weighted_modres, self.weighted_ssa, self.weighted_phi,
             self.weighted_psi, self.weighted_intra_ssbond, self.weighted_ssbond_length, self.weighted_intra_link,
             self.weighted_link_length, self.weighted_cis_conformation, self.weighted_cis_follower,
             self.weighted_inter_chain_median_kd, self.weighted_inter_chain_dist_weighted_kd, self.weighted_inter_chain_median_rsa,
             self.weighted_inter_chain_dist_weighted_rsa, self.weighted_intra_chain_median_kd, self.weighted_intra_chain_dist_weighted_kd,
             self.weighted_intra_chain_median_rsa, self.weighted_intra_chain_dist_weighted_rsa,
             self.weighted_inter_chain_interactions_median, self.weighted_inter_chain_interactions_dist_weighted,
             self.weighted_intra_chain_interactions_median, self.weighted_intra_chain_interactions_dist_weighted,
             self.weighted_lig_dist, self.weighted_metal_dist, self.weighted_ion_dist, self.weighted_chain_dist, self.weighted_rna_dist, self.weighted_dna_dist,
             self.weighted_homo_dist, self.rin_class, self.rin_simple_class)

    def get_database_result_for_indel_aggregation(self):
        return (self.weighted_surface_value, self.weighted_mainchain_surface_value, self.weighted_sidechain_surface_value,
             self.get_weighted_profile_str(), self.get_weighted_centralities_str(), self.weighted_b_factor, self.weighted_modres, self.weighted_ssa, self.weighted_phi,
             self.weighted_psi, self.weighted_intra_ssbond, self.weighted_ssbond_length, self.weighted_intra_link,
             self.weighted_link_length, self.weighted_cis_conformation, self.weighted_cis_follower,
             self.weighted_inter_chain_median_kd, self.weighted_inter_chain_dist_weighted_kd, self.weighted_inter_chain_median_rsa,
             self.weighted_inter_chain_dist_weighted_rsa, self.weighted_intra_chain_median_kd, self.weighted_intra_chain_dist_weighted_kd,
             self.weighted_intra_chain_median_rsa, self.weighted_intra_chain_dist_weighted_rsa,
             self.weighted_inter_chain_interactions_median, self.weighted_inter_chain_interactions_dist_weighted,
             self.weighted_intra_chain_interactions_median, self.weighted_intra_chain_interactions_dist_weighted,
             self.rin_class, self.rin_simple_class)


    def get_result_copy(self):
        result_obj = Mappings()
        result_obj.recommended_res = self.recommended_res
        result_obj.max_seq_res = self.max_seq_res
        result_obj.weighted_surface_value = self.weighted_surface_value
        result_obj.weighted_mainchain_surface_value = self.weighted_mainchain_surface_value
        result_obj.weighted_sidechain_surface_value = self.weighted_sidechain_surface_value
        result_obj.weighted_location = self.weighted_location
        result_obj.weighted_mainchain_location = self.weighted_mainchain_location
        result_obj.weighted_sidechain_location = self.weighted_sidechain_location
        result_obj.weighted_profile = self.weighted_profile
        result_obj.weighted_centralities = self.weighted_centralities
        result_obj.weighted_b_factor = self.weighted_b_factor
        result_obj.weighted_modres = self.weighted_modres
        result_obj.weighted_ssa = self.weighted_ssa
        result_obj.weighted_phi = self.weighted_phi
        result_obj.weighted_psi = self.weighted_psi
        result_obj.weighted_intra_ssbond = self.weighted_intra_ssbond
        result_obj.weighted_inter_ssbond = self.weighted_inter_ssbond
        result_obj.weighted_ssbond_length = self.weighted_ssbond_length
        result_obj.weighted_intra_link = self.weighted_intra_link
        result_obj.weighted_inter_link = self.weighted_inter_link
        result_obj.weighted_link_length = self.weighted_link_length
        result_obj.weighted_cis_conformation = self.weighted_cis_conformation
        result_obj.weighted_cis_follower = self.weighted_cis_follower
        result_obj.weighted_inter_chain_median_kd = self.weighted_inter_chain_median_kd
        result_obj.weighted_inter_chain_dist_weighted_kd = self.weighted_inter_chain_dist_weighted_kd
        result_obj.weighted_inter_chain_median_rsa = self.weighted_inter_chain_median_rsa
        result_obj.weighted_inter_chain_dist_weighted_rsa = self.weighted_inter_chain_dist_weighted_rsa
        result_obj.weighted_intra_chain_median_kd = self.weighted_intra_chain_median_kd
        result_obj.weighted_intra_chain_dist_weighted_kd = self.weighted_intra_chain_dist_weighted_kd
        result_obj.weighted_intra_chain_median_rsa = self.weighted_intra_chain_median_rsa
        result_obj.weighted_intra_chain_dist_weighted_rsa = self.weighted_intra_chain_dist_weighted_rsa
        result_obj.weighted_inter_chain_interactions_median = self.weighted_inter_chain_interactions_median
        result_obj.weighted_inter_chain_interactions_dist_weighted = self.weighted_inter_chain_interactions_dist_weighted
        result_obj.weighted_intra_chain_interactions_median = self.weighted_intra_chain_interactions_median
        result_obj.weighted_intra_chain_interactions_dist_weighted = self.weighted_intra_chain_interactions_dist_weighted
        result_obj.weighted_lig_dist = self.weighted_lig_dist
        result_obj.weighted_metal_dist = self.weighted_metal_dist
        result_obj.weighted_ion_dist = self.weighted_ion_dist
        result_obj.weighted_chain_dist = self.weighted_chain_dist
        result_obj.weighted_rna_dist = self.weighted_rna_dist
        result_obj.weighted_dna_dist = self.weighted_dna_dist
        result_obj.weighted_homo_dist = self.weighted_homo_dist
        result_obj.rin_class = self.rin_class
        result_obj.rin_simple_class = self.rin_simple_class
        result_obj.Class = self.Class
        result_obj.simple_class = self.simple_class
        result_obj.interaction_recommendations = self.interaction_recommendations
        result_obj.lig_dist_conf = self.lig_dist_conf
        result_obj.metal_dist_conf = self.metal_dist_conf
        result_obj.ion_dist_conf = self.ion_dist_conf
        result_obj.chain_dist_conf = self.chain_dist_conf
        result_obj.rna_dist_conf = self.rna_dist_conf
        result_obj.dna_dist_conf = self.dna_dist_conf
        result_obj.homo_dist_conf = self.homo_dist_conf
        result_obj.location_conf = self.location_conf
        result_obj.classification_conf = self.classification_conf
        result_obj.amount_of_structures = self.amount_of_structures
        result_obj.rsa_change_score = self.rsa_change_score
        result_obj.mc_rsa_change_score = self.mc_rsa_change_score
        result_obj.sc_rsa_change_score = self.sc_rsa_change_score
        return result_obj

    def weight_all(self, config, disorder_score, disorder_region, for_indel_aggregation = False):
        dist_maps = [(self.lig_dists, self.set_weighted_lig_dist), (self.metal_dists, self.set_weighted_metal_dist),
                     (self.ion_dists, self.set_weighted_ion_dist), (self.chain_dists, self.set_weighted_chain_dist),
                     (self.rna_dists, self.set_weighted_rna_dist), (self.dna_dists, self.set_weighted_dna_dist),
                     (self.homo_dists, self.set_weighted_homo_dist)]
        for dist_map, func in dist_maps:
            weighted_value, conf = self.weight(dist_map, config, distance_weighting=True)
            func(weighted_value, conf)

        value_maps = [(self.phis, self.set_weighted_phi), (self.psis, self.set_weighted_psi), (self.ssbond_lengths, self.set_weighted_ssbond_length),
                      (self.link_lengths, self.set_weighted_link_length), (self.inter_chain_median_kds, self.set_weighted_inter_chain_median_kd),
                      (self.inter_chain_dist_weighted_kds, self.set_weighted_inter_chain_dist_weighted_kd),
                      (self.inter_chain_median_rsas, self.set_weighted_inter_chain_median_rsa),
                      (self.inter_chain_dist_weighted_rsas, self.set_weighted_inter_chain_dist_weighted_rsa),
                      (self.intra_chain_median_kds, self.set_weighted_intra_chain_median_kd),
                      (self.intra_chain_dist_weighted_kds, self.set_weighted_intra_chain_dist_weighted_kd),
                      (self.intra_chain_median_rsas, self.set_weighted_intra_chain_median_rsa),
                      (self.intra_chain_dist_weighted_rsas, self.set_weighted_intra_chain_dist_weighted_rsa),
                      (self.inter_chain_interactions_medians, self.set_weighted_inter_chain_interactions_median),
                      (self.inter_chain_interactions_dist_weighteds, self.set_weighted_inter_chain_interactions_dist_weighted),
                      (self.intra_chain_interactions_medians, self.set_weighted_intra_chain_interactions_median),
                      (self.intra_chain_interactions_dist_weighteds, self.set_weighted_intra_chain_interactions_dist_weighted),
                      (self.b_factors, self.set_weighted_b_factor)]
        for value_map, func in value_maps:
            weighted_value = self.weight(value_map, config, calc_conf=False)
            func(weighted_value)

        prop_maps = [(self.cis_conformations, self.set_weighted_cis_conformation), (self.cis_followers, self.set_weighted_cis_follower),
                     (self.modres, self.set_weighted_modres)]
        for prop_map, func in prop_maps:
            weighted_prop = self.weight_prop(prop_map)
            func(weighted_prop)

        bool_prop_maps = [(self.intra_ssbonds, self.set_weighted_intra_ssbond, self.set_weighted_inter_ssbond),
                          (self.intra_links, self.set_weighted_intra_link, self.set_weighted_inter_link)]
        for prop_map, true_func, false_func in bool_prop_maps:
            weighted_true_prop, weighted_false_prop = self.weigthed_bool_prop(prop_map)
            true_func(weighted_true_prop)
            false_func(weighted_false_prop)

        self.weighted_ssa = self.weight_majority(self.ssas)
        if len(self.covs) > 0:
            self.max_cov = max(self.covs.values())
        self.weight_structure_location(config)
        self.weight_centralities()
        self.weight_profiles()
        self.rin_based_classify()
        self.distance_classify(config, disorder_score, disorder_region)
        if not for_indel_aggregation:
            self.set_recommended_residues()  
        
        #analyse for conformational changes
        #print("analyze conformational changes")
        core_count = 0
        buried_count = 0
        surface_count = 0
        #core_list = []
        #buried_list = []
        #surface_list = []
        self.rsa_change_score = 0
        for mapping_id in self.rsas:
            #mapping_id = pdb_id, chain, res_nr
            rsa = self.rsas[mapping_id]
            location = locate(rsa, config)
            #b_factor = self.b_factors[mapping_id]
            #resolution = self.resolutions[mapping_id[0]]
            
            #criteria if residue is to be considered
            #if resolution > 4:
                #continue
            #if b_factor < 30:
                #continue

            if location == "Surface":
                surface_count +=1
                #surface_list.append((mapping_id, self.res_aas[mapping_id]))
            elif location == "Core" or location == "Buried":
                core_count += 1
                #core_list.append((mapping_id, self.res_aas[mapping_id]))
            #else:
            #    print(f"Unknown structural classification: {location}")
        
            """
            elif location == "Buried":
                buried_count += 1
                buried_list.append((mapping_id, self.res_aas[mapping_id]))
            """
        
        mc_core_count = 0
        mc_surface_count = 0
        #mc_core_list = []
        #mc_surface_list = []
        self.mc_rsa_change_score = 0
        for mapping_id in self.mc_rsas:
            mc_rsa = self.mc_rsas[mapping_id]
            mc_location = locate(mc_rsa, config)
            if mc_location == "Surface":
                mc_surface_count +=1
                #mc_surface_list.append((mapping_id, self.res_aas[mapping_id]))
            elif mc_location == "Core" or mc_location == "Buried":
                mc_core_count += 1
                #mc_core_list.append((mapping_id, self.res_aas[mapping_id]))
            #else:
            #    print(f"Unknown structural classification: {mc_location}")
        
        sc_core_count = 0
        sc_surface_count = 0
        #sc_core_list = []
        #sc_surface_list = []
        self.sc_rsa_change_score = 0
        for mapping_id in self.sc_rsas:
            sc_rsa = self.sc_rsas[mapping_id]
            sc_location = locate(sc_rsa, config)
            if sc_location == "Surface":
                sc_surface_count +=1
                #sc_surface_list.append((mapping_id, self.res_aas[mapping_id]))
            elif sc_location == "Core" or sc_location == "Buried":
                sc_core_count += 1
                #sc_core_list.append((mapping_id, self.res_aas[mapping_id]))
            #else:
            #    print(f"Unknown structural classification: {sc_location}")
        
        #criteria for possible conf change
        if core_count > 0 and surface_count > 0:
            total_sum = core_count + surface_count
            self.rsa_change_score = core_count / total_sum
            """
            #print(f"core {core_count} surface {surface_count}, score {self.rsa_change_score}")
            #sort surface and core lists by second element (amino acid)
            #core
            #sort aa alphabetically
            core_list.sort(key = lambda x: x[1])
            #frequencies of aa
            freq_core = Counter(val for key, val in core_list)
            #sort by frequencies
            sorted_core=sorted(core_list, key=lambda element: freq_core[element[1]], reverse=True)
            occurences_core = Counter([aa for (m_id, aa) in sorted_core])

            #surface
            surface_list.sort(key = lambda x:x[1])
            freq_surface = Counter(val for key, val in surface_list)
            sorted_surface = sorted(surface_list, key=lambda element: freq_surface[element[1]], reverse=True)
            occurences_surface = Counter([aa for (m_id, aa) in sorted_surface])

            occ_core = []
            for i in occurences_core.values():
                occ_core.append(i)
            occ_surface = []
            for i in occurences_surface.values():
                occ_surface.append(i)

            
            #put aa before occurrences in array[]
            #core
            curr_aa = None
            core_out = []
            for x in sorted_core:
                if x[1] == curr_aa:
                    core_out.append(x[0])
                else:
                    curr_aa=x[1]
                    core_out.append(x[1])
                    core_out.append(x[0])
        	#surface
            curr_aa = None
            surface_out = []
            for x in sorted_surface:
                if x[1] == curr_aa:
                    surface_out.append(x[0])
                else:
                    curr_aa=x[1]
                    surface_out.append(x[1])
                    surface_out.append(x[0])
            
            output = []
            output.append(" core")
            j = 0
            for i in core_out:
                if len(i) == 1:
                    output.append('\n')
                    output.append(i)
                    output.append(str(occ_core[j]))
                    j += 1
                else:
                    output.append(i[0])
            output.append('\n')
            output.append("surface")
            l = 0
            for k in surface_out:
                if len(k) == 1:
                    output.append('\n')
                    output.append(k)
                    output.append(str(occ_surface[l]))
                    l += 1
                else:
                    output.append(k[0])
        			
            what_is_s = " ".join(output)
        	#print(f"res nr.{mapping_id[2]}\n {core_count} as core, {surface_count} as surface \n {s}")
        	"""
        	
        elif core_count == 0 and surface_count > 0:
            #everything is surface (or buried, but we don"t make a difference yet)
            self.rsa_change_score = 0
            """
            #surface
            surface_list.sort(key = lambda x:x[1])
            freq_surface = Counter(val for key, val in surface_list)
            sorted_surface = sorted(surface_list, key=lambda element: freq_surface[element[1]], reverse=True)
            occurences_surface = Counter([aa for (m_id, aa) in sorted_surface])
        
            occ_surface = []
            for i in occurences_surface.values():
                occ_surface.append(i)

            curr_aa = None
            surface_out = []
            for x in sorted_surface:
                if x[1] == curr_aa:
                    surface_out.append(x[0])
                else:
                    curr_aa=x[1]
                    surface_out.append(x[1])
                    surface_out.append(x[0])
            output = []
            output.append("surface")
            l = 0
            for k in surface_out:
                if len(k) == 1:
                    output.append('\n')
                    output.append(k)
                    output.append(str(occ_surface[l]))
                    l += 1
                else:
                    output.append(k[0])
                    
            what_is_s = " ".join(output)
            """
        	
        	
        elif core_count > 0 and surface_count == 0:
            #everything is core (or buried, see above)
            self.rsa_change_score = 1
            
            """            
            core_list.sort(key = lambda x: x[1])
            freq_core = Counter(val for key, val in core_list)
            sorted_core=sorted(core_list, key=lambda element: freq_core[element[1]], reverse=True)
            occurences_core = Counter([aa for (m_id, aa) in sorted_core])
            
            occ_core = []
            for i in occurences_core.values():
                occ_core.append(i)

            curr_aa = None
            core_out = []
            for x in sorted_core:
                if x[1] == curr_aa:
                    core_out.append(x[0])
                else:
                    curr_aa=x[1]
                    core_out.append(x[1])
                    core_out.append(x[0])
                    
            output = []
            output.append(" core")
            j = 0
            for i in core_out:
                if len(i) == 1:
                    output.append('\n')
                    output.append(i)
                    output.append(str(occ_core[j]))
                    j += 1
                else:
                    output.append(i[0])
                    
            what_is_s = " ".join(output)
            """
            
        #for mainchain and sidechain only rsa change score(rcs), lists with amino acids and occurences can be added later
        #mainchain rcs
        if mc_core_count > 0 and mc_surface_count > 0:
            total_sum = mc_core_count + mc_surface_count
            self.mc_rsa_change_score = mc_core_count / total_sum
            
        elif mc_core_count == 0 and mc_surface_count > 0:
            #everything is surface
            self.mc_rsa_change_score = 0
            
        elif mc_core_count > 0 and mc_surface_count == 0:
            #everything is core
            self.mc_rsa_change_score = 1
            
        #sidechain rcs
        if sc_core_count > 0 and sc_surface_count > 0:
            total_sum = sc_core_count + sc_surface_count
            self.sc_rsa_change_score = sc_core_count / total_sum
            
        elif sc_core_count == 0 and sc_surface_count > 0:
            #everything is surface
            self.sc_rsa_change_score = 0
            
        elif sc_core_count > 0 and sc_surface_count == 0:
            #everything is core
            self.sc_rsa_change_score = 1
        
        

    def weight(self, value_map, config, distance_weighting=False, calc_conf=True, coverage_extra_weight=0):
        nom = 0.0
        denom = 0.0
        n = 0.0
        qs = []
        weighted_value = None
        conf = 0.0
        for mapping_id in value_map:
            value = value_map[mapping_id]
            if value is None:
                continue
            if distance_weighting and value > config.long_distance_threshold:
                continue
            if not isinstance(value, float):
                print('Strange value: ', value)
            qual = self.qualities[mapping_id]
            cov = self.covs[mapping_id]
            weight = qual * (cov ** coverage_extra_weight)
            nom += value * weight
            denom += weight
            n += 1.0
            qs.append(weight)
        if denom > 0.0:
            weighted_value = nom / denom
            if calc_conf:
                conf = (1.0 - 1.0 / (n + 1.0)) * (max(qs) + median(qs)) / 2
        if calc_conf:
            return weighted_value, conf
        else:
            return weighted_value

    def weight_prop(self, value_map):
        prop = 0.
        if len(value_map) == 0:
            return None
        for mapping_id in value_map:
            value = value_map[mapping_id]
            if value is None:
                continue
            prop += 1.
        weighted_prop = prop / float(len(value_map))
        return weighted_prop

    def weigthed_bool_prop(self, value_map):
        true_prop = 0.
        false_prop = 0.
        if len(value_map) == 0:
            return None, None
        for mapping_id in value_map:
            value = value_map[mapping_id]
            if value is None:
                continue
            if value:
                true_prop += 1.
            else:
                false_prop += 1.
        weighted_true_prop = true_prop / float(len(value_map))
        weighted_false_prop = false_prop / float(len(value_map))
        return weighted_true_prop, weighted_false_prop

    def weight_majority(self, value_map):
        voting = {}
        best_value = None
        for mapping_id in value_map:
            qual = self.qualities[mapping_id]
            value = value_map[mapping_id]
            if value not in voting:
                voting[value] = qual
            else:
                voting[value] += qual

        max_qual = 0.
        for value in voting:
            qual = voting[value]
            if qual > max_qual:
                max_qual = qual
                best_value = value

        return best_value

    def weight_structure_location(self, config):
        self.weighted_surface_value, conf = self.weight(self.rsas, config, coverage_extra_weight=2)
        self.weighted_mainchain_surface_value, mc_conf = self.weight(self.mc_rsas, config, coverage_extra_weight=2)
        self.weighted_sidechain_surface_value, sc_conf = self.weight(self.sc_rsas, config, coverage_extra_weight=2)

        self.location_conf = sc_conf

        (self.weighted_location, self.weighted_mainchain_location, self.weighted_sidechain_location) = triple_locate(self.weighted_surface_value,
                                                                                                                     self.weighted_mainchain_surface_value, self.weighted_sidechain_surface_value, config)

    def weight_centralities(self):
        self.weighted_centralities = None

        total_qual = 0.0

        for mapping_id in self.centralities:
            centrality_scores = self.centralities[mapping_id]
            if centrality_scores is None:
                continue
            qual = self.qualities[mapping_id]
            if self.weighted_centralities is None:
                self.weighted_centralities = [0.] * len(centrality_scores.cent_list)

            for pos, cent_score in enumerate(centrality_scores.cent_list):
                if cent_score is None:
                    continue
                self.weighted_centralities[pos] += cent_score * qual

            total_qual += qual
        if self.weighted_centralities is None:
            return

        if total_qual > 0.0:
            for pos, cent_score in enumerate(self.weighted_centralities):
                self.weighted_centralities[pos] = self.weighted_centralities[pos] / total_qual
            self.weighted_centralities = rin.Centrality_scores(cent_list=self.weighted_centralities)

    def get_weighted_centralities(self):
        if self.weighted_centralities is not None:
            return self.weighted_centralities
        if self.weighted_centralities_str is None:
            return None
        return rin.Centrality_scores(code_str=self.weighted_centralities_str)

    def get_weighted_centralities_str(self):
        if self.weighted_centralities_str is not None:
            return self.weighted_centralities_str
        if self.weighted_centralities is None:
            return None
        return self.weighted_centralities.str_encode()

    def get_weighted_profile(self):
        if self.weighted_profile is not None:
            return self.weighted_profile
        if self.weighted_profile_str is None:
            return None
        return rin.Interaction_profile(profile_str=self.weighted_profile_str)

    def get_weighted_profile_str(self):
        if self.weighted_profile_str is not None:
            return self.weighted_profile_str
        if self.weighted_profile is None:
            return None
        return self.weighted_profile.encode()

    def weight_profiles(self):
        weight_profile_tuples = []
        for mapping_id in self.profiles:
            weight_profile_tuples.append((self.qualities[mapping_id], self.profiles[mapping_id]))
        self.weighted_profile = rin.calculateAverageProfile(weight_profile_tuples)

    def set_recommended_residues(self):
        max_seq_id = 0.
        max_identical_aa_identical_class = 0.
        max_identical_aa = 0.
        max_identical_class = 0.
        max_should_not_happen = 0.
        max_qual = 0.
        max_qual_res = None
        self.max_seq_res = None
        self.recommended_res = None
        self.interaction_recommendations = {}

        identical_aa_identical_class_recs = []
        identical_class_recs = []
        identical_aa_recs = []
        other_recs = []

        for mapping_id in self.seq_ids:
            seq_id = self.seq_ids[mapping_id]
            cov = self.covs[mapping_id]
            qual = self.qualities[mapping_id]
            if qual is None:
                qual = 0.

            if seq_id > max_seq_id:
                max_seq_id = seq_id
                self.max_seq_res = qual, (mapping_id, seq_id, cov)
            elif seq_id == max_seq_id:
                if mapping_id[0] > self.max_seq_res[1][0][0]:
                    max_seq_id = seq_id
                    self.max_seq_res = qual, (mapping_id, seq_id, cov)

            if qual > max_qual:
                max_qual = qual
                max_qual_res = mapping_id, seq_id, cov
            elif qual == max_qual:
                if mapping_id[0] > max_qual_res[0][0]:
                    max_qual = qual
                    max_qual_res = mapping_id, seq_id, cov

            residue_simple_class = self.res_classes[mapping_id][1]

            if self.rin_simple_class == residue_simple_class:
                if self.aa_ids[mapping_id]:
                    identical_aa_identical_class_recs.append((qual, mapping_id))
                    if qual > max_identical_aa_identical_class:
                        max_identical_aa_identical_class = qual
                        max_identical_aa_identical_class_res = mapping_id, seq_id, cov
                    elif qual == max_identical_aa_identical_class:
                        if mapping_id[0] > max_identical_aa_identical_class_res[0][0]:
                            max_identical_aa_identical_class = qual
                            max_identical_aa_identical_class_res = mapping_id, seq_id, cov
                else:
                    identical_class_recs.append((qual, mapping_id))
                    if qual > max_identical_class:
                        max_identical_class = qual
                        max_identical_class_res = mapping_id, seq_id, cov
                    elif qual == max_identical_class:
                        if mapping_id[0] > max_identical_class_res[0][0]:
                            max_identical_class = qual
                            max_identical_class_res = mapping_id, seq_id, cov
            else:
                if self.aa_ids[mapping_id]:
                    identical_aa_recs.append((qual, mapping_id))
                    if qual > max_identical_aa:
                        max_identical_aa = qual
                        max_identical_aa_res = mapping_id, seq_id, cov
                    elif qual == max_identical_aa:
                        if mapping_id[0] > max_identical_aa_res[0][0]:
                            max_identical_aa = qual
                            max_identical_aa_res = mapping_id, seq_id, cov
                else:
                    other_recs.append((qual, mapping_id))
                    if qual > max_should_not_happen:
                        max_should_not_happen = qual
                        max_should_not_happen_res = mapping_id, seq_id, cov
                    elif qual == max_should_not_happen:
                        if mapping_id[0] > max_should_not_happen_res[0][0]:
                            max_should_not_happen = qual
                            max_should_not_happen_res = mapping_id, seq_id, cov

            if residue_simple_class is not None:
                if residue_simple_class.count('Interaction') > 0:
                    if residue_simple_class not in self.interaction_recommendations:
                        self.interaction_recommendations[residue_simple_class] = (qual, mapping_id, seq_id, cov)
                    elif qual > self.interaction_recommendations[residue_simple_class][0]:
                        self.interaction_recommendations[residue_simple_class] = (qual, mapping_id, seq_id, cov)

        recommendation_order = sorted(identical_aa_identical_class_recs, key = lambda x:x[0], reverse = True)
        recommendation_order += sorted(identical_class_recs, key = lambda x:x[0], reverse = True)
        recommendation_order += sorted(identical_aa_recs, key = lambda x:x[0], reverse = True)
        recommendation_order += sorted(other_recs, key = lambda x:x[0], reverse = True)
        #print(f'Recommendation order in mappings class: {recommendation_order}')
        self.recommendation_order = recommendation_order

        if max_identical_aa_identical_class > 0.:
            self.recommended_res = max_identical_aa_identical_class, max_identical_aa_identical_class_res
        elif max_identical_class > 0.:
            self.recommended_res = max_identical_class, max_identical_class_res
        elif max_identical_aa > 0.:
            self.recommended_res = max_identical_aa, max_identical_aa_res
        else:
            # If this case is reached, then max_qual_res should be identical to max_should_not_happen_res
            # This would mean there is not mapped residue with identical simple class or identical aa, which should not happen (thus the name)
            self.recommended_res = max_qual, max_qual_res

        if self.recommended_res is not None:
            if self.recommended_res[1] is not None:
                qual, (recommended_res, seq_id, cov) = self.recommended_res
                pdb_id, chain, res_nr = recommended_res
                resolution = self.resolutions[pdb_id]
                res_aa = self.res_aas[(pdb_id, chain, res_nr)]
                self.recommended_res = '%s:%s %s:%s;%1.2f;%1.2f;%1.1f' % (pdb_id, chain, res_nr, res_aa, seq_id, cov, resolution)
            else:
                self.recommended_res = None

        if self.max_seq_res is not None:
            if self.max_seq_res[1] is not None:
                qual, (recommended_res, seq_id, cov) = self.max_seq_res
                pdb_id, chain, res_nr = recommended_res
                resolution = self.resolutions[pdb_id]
                res_aa = self.res_aas[(pdb_id, chain, res_nr)]
                self.max_seq_res = '%s:%s %s:%s;%1.2f;%1.2f;%1.1f' % (pdb_id, chain, res_nr, res_aa, seq_id, cov, resolution)
            else:
                self.max_seq_res = None

    def get_recommended_res_str(self):
        return self.recommended_res

    def get_max_seq_structure_res_str(self):
        return self.max_seq_res

    def distance_classify(self, config, disorder_score, disorder_region):
        self.Class, self.classification_conf = self.getWeightedClass(config)
        if self.Class is None:
            self.Class, self.classification_conf = self.disorderCheck(config, disorder_score, disorder_region)
            self.simple_class = self.Class
            self.rin_class = self.Class
            self.rin_simple_class = self.Class
        else:
            self.simple_class = self.simplifyClass(self.Class, self.weighted_location)

    def rin_based_classify(self):
        rin_class, rin_simple_class = rin_classify(self.weighted_profile, self.weighted_sidechain_location)
        self.rin_class = rin_class
        self.rin_simple_class = rin_simple_class

    def simplifyClass(self, c, sc):
        if c == "Surface" or c == "Core" or c == 'Disorder' or c == 'Buried' or c is None:
            return c
        interactions = re.sub(r' Interaction$', '', re.sub(r'^[^:]*: ', '', c)).split(' and ')
        non_far_interactions = set([x for x in interactions if ' far' not in x])

        priorities = ['Metal', 'Ligand', 'DNA', 'RNA', 'Protein', 'Ion']

        if len(non_far_interactions) == 0:
            return sc

        for interaction in priorities:
            if interaction in non_far_interactions:
                return interaction + ' Interaction'

        print('Unknown Class: %s' % c)
        return sc

    def disorderCheck(self, config, disorder_score, region):
        if region is None or disorder_score is None:
            return None, None
        glob = True
        if region == 'disorder':
            glob = False
        elif region == 'globular':
            glob = True
        elif region in MobiDB_map:
            glob = False
        elif config.verbosity >= 3:
            print('Unknown disorder region: ', region)
        if not glob:
            return 'Disorder', disorder_score
        else:
            return None, None

    def getWeightedClass(self, config):
        dt = config.short_distance_threshold

        contacts = []
        confs = []

        # Determine the near macromolecule with the highest conf.
        # If there is no near macromolecule, also consider far ones.
        # In case of multiple maximum confs use the first one added to macros.
        macros = []
        macro_confs = []

        if self.weighted_chain_dist and self.weighted_chain_dist < dt:
            macros.append('Protein')
            macro_confs.append(self.chain_dist_conf)

        if self.weighted_dna_dist and self.weighted_dna_dist < dt:
            macros.append('DNA')
            macro_confs.append(self.dna_dist_conf)

        if self.weighted_rna_dist and self.weighted_rna_dist < dt:
            macros.append('RNA')
            macro_confs.append(self.rna_dist_conf)

        contact_macro = None
        conf_macro = None

        if len(macros) > 0:
            best_macro_idx = np.min(np.argmax(macro_confs))
            contacts.append(macros[best_macro_idx])
            confs.append(macro_confs[best_macro_idx])

        # Additional interactions
        if self.weighted_lig_dist and self.weighted_lig_dist < dt:
            contacts.append('Ligand')
            confs.append(self.lig_dist_conf)

        if self.weighted_metal_dist and self.weighted_metal_dist < dt:
            contacts.append('Metal')
            confs.append(self.metal_dist_conf)

        if self.weighted_ion_dist and self.weighted_ion_dist < dt:
            contacts.append('Ion')
            confs.append(self.ion_dist_conf)

        if len(contacts) == 0:
            return self.weighted_location, self.location_conf

        clas = " and ".join(contacts)
        if len(contacts) == 4:
            clas = "Quadruple Interaction: " + clas
        elif len(contacts) == 3:
            clas = "Triple Interaction: " + clas
        elif len(contacts) == 2:
            clas = "Double Interaction: " + clas
        elif len(contacts) == 1:
            clas = clas + " Interaction"

        conf = sum(confs) / len(confs)

        return clas, conf

    def set_mm_feature(self, feat_value, feature_name):
        setattr(self, feature_name.replace(' ','_'), feat_value)