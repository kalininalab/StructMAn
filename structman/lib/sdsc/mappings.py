import statistics

# sdsc: structman datastructures and classes

from structman.lib import rin
from structman.base_utils.base_utils import pack
from structman.lib.sdsc.sdsc_utils import classify, triple_locate, doomsday_protocol, locate

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


def weight(value_map, config, qualities, covs, distance_weighting=False, calc_conf=True, coverage_extra_weight=0):
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
        qual = qualities[mapping_id]
        cov = covs[mapping_id]
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


def weight_propensity(value_map):
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


def weigthed_bool_propensity(value_map):
    true_prop = 0.

    if len(value_map) == 0:
        return None
    for mapping_id in value_map:
        value = value_map[mapping_id]
        if value is None:
            continue
        if value:
            true_prop += 1.

    weighted_true_prop = true_prop / float(len(value_map))
    return weighted_true_prop


def weight_majority(value_map, qualities):
    voting = {}
    best_value = None
    for mapping_id in value_map:
        qual = qualities[mapping_id]
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

class Feature_set:
    def __init__(self):
        for feature_name in self.__slots__:
            setattr(self, feature_name, None)
        return
    
    def set_values(self, values):
        if values is None:
            return
        for feat_pos, feature_name in enumerate(self.__slots__):
            setattr(self, feature_name, values[feat_pos])

    def set_value_by_name(self, feature_name, value):
        setattr(self, feature_name, value)

    def get_raw_list(self):
        values = []
        for feature_name in self.__slots__:
            values.append(getattr(self, feature_name))
        return values
    
    def get_feature_names(self):
        return self.__slots__

    def add_to_output_object(self, output_object):
        for feature_name in self.__slots__:
            output_object.add_value(feature_name, getattr(self, feature_name))

    def create_unweighted_feature_dicts(self):
        unweighted_feature_dicts = {}
        for feature_name in self.__slots__:
            unweighted_feature_dict = {}
            unweighted_feature_dicts[feature_name] = unweighted_feature_dict
        return unweighted_feature_dicts

class Microminer_features(Feature_set):
    __slots__ = [
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

class Structural_features(Feature_set):
    __slots__ = [
        'b_factor', 'modres', 'ssa', 'phi', 'psi', 'intra_ssbond', 'inter_ssbond', 'ssbond_length', 'intra_link', 'inter_link', 'link_length',
        'cis_conformation', 'cis_follower', 'inter_chain_median_kd', 'inter_chain_dist_weighted_kd',
        'inter_chain_median_rsa', 'inter_chain_dist_weighted_rsa', 'intra_chain_median_kd', 'intra_chain_dist_weighted_kd', 'intra_chain_median_rsa', 'intra_chain_dist_weighted_rsa',
        'inter_chain_interactions_median', 'inter_chain_interactions_dist_weighted', 'intra_chain_interactions_median', 'intra_chain_interactions_dist_weighted',
        'lig_dist', 'metal_dist', 'ion_dist', 'chain_dist', 'rna_dist', 'dna_dist', 'homo_dist', 'surface_value', 'mainchain_surface_value', 'sidechain_surface_value',
        'rin_class', 'rin_simple_class',
    ]

    distance_based_features = set([
        'lig_dist', 'metal_dist', 'ion_dist', 'chain_dist', 'rna_dist', 'dna_dist', 'homo_dist'
    ])

    float_features = set([
        'b_factor', 'phi', 'psi', 'ssbond_length', 'link_length', 'inter_chain_median_kd', 'inter_chain_dist_weighted_kd',
        'inter_chain_median_rsa', 'inter_chain_dist_weighted_rsa', 'intra_chain_median_kd', 'intra_chain_dist_weighted_kd', 'intra_chain_median_rsa', 'intra_chain_dist_weighted_rsa',
        'inter_chain_interactions_median', 'inter_chain_interactions_dist_weighted', 'intra_chain_interactions_median', 'intra_chain_interactions_dist_weighted'
    ])

    probability_features = set([
        'modres', 'cis_conformation', 'cis_follower'
    ])

    boolean_features = set([
        'intra_ssbond', 'inter_ssbond', 'intra_link', 'inter_link',
    ])

    categorical_features = set([
        'ssa'
    ])

    surface_value_features = set([
        'surface_value', 'mainchain_surface_value', 'sidechain_surface_value'
    ])

    def weight_function(self, feature_name, value_map, qualities, covs, config):
        if feature_name in Structural_features.distance_based_features:
            return weight(value_map, config, qualities, covs, distance_weighting=True, calc_conf=False)
        if feature_name in Structural_features.float_features:
            return weight(value_map, config, qualities, covs, calc_conf=False)
        if feature_name in Structural_features.probability_features:
            return weight_propensity(value_map)
        if feature_name in Structural_features.boolean_features:
            return weigthed_bool_propensity(value_map)
        if feature_name in Structural_features.categorical_features:
            return weight_majority(value_map, qualities)
        if feature_name in Structural_features.surface_value_features:
            return weight(value_map, config, qualities, covs, calc_conf=False, coverage_extra_weight=2)
        return None

class Integrated_features(Feature_set):
    __slots__ = [
        'location', 'mainchain_location', 'sidechain_location',
        'structural_classification', 'simple_class',
        'rsa_change_score', 'mc_rsa_change_score', 'sc_rsa_change_score', 'sc_rsa_std'
    ]

    def generate(self, structural_features, config, profile, disorder_score, disorder_region, surface_values, mc_surface_values, sc_surface_values):
        (self.location, self.mainchain_location, self.sidechain_location) = triple_locate(structural_features.surface_value, structural_features.mainchain_surface_value, structural_features.sidechain_surface_value, config)
        
        self.structural_classification, self.simple_class = classify(config, profile, self.sidechain_location, disorder_score, disorder_region)

        #analyse for conformational changes
        core_count = 0
        surface_count = 0
        #core_list = []
        #buried_list = []
        #surface_list = []
        self.rsa_change_score = None
        for mapping_id in surface_values:
            rsa = surface_values[mapping_id]
            location = locate(rsa, config)

            if location == "Surface":
                surface_count +=1
            elif location == "Core" or location == "Buried":
                core_count += 1
        
        mc_core_count = 0
        mc_surface_count = 0

        self.mc_rsa_change_score = None
        for mapping_id in mc_surface_values:
            mc_rsa = mc_surface_values[mapping_id]
            mc_location = locate(mc_rsa, config)
            if mc_location == "Surface":
                mc_surface_count +=1
            elif mc_location == "Core" or mc_location == "Buried":
                mc_core_count += 1

        
        sc_core_count = 0
        sc_surface_count = 0
        #sc_core_list = []
        sc_surface_list = []
        self.sc_rsa_change_score = None
        for mapping_id in sc_surface_values:
            sc_rsa = sc_surface_values[mapping_id]
            if sc_rsa is not None:
                sc_surface_list.append(sc_rsa)
            sc_location = locate(sc_rsa, config)
            if sc_location == "Surface":
                sc_surface_count +=1
            elif sc_location == "Core" or sc_location == "Buried":
                sc_core_count += 1


        if len(sc_surface_list) > 1:
            self.sc_rsa_std = statistics.stdev(sc_surface_list)
        elif len(sc_surface_list) == 1:
            self.sc_rsa_std = 0.
        else:
            self.sc_rsa_std = None

        #criteria for possible conf change
        if core_count > 0 and surface_count > 0:
            total_sum = core_count + surface_count
            self.rsa_change_score = core_count / total_sum
        	
        elif core_count == 0 and surface_count > 0:
            #everything is surface (or buried, but we don"t make a difference yet)
            self.rsa_change_score = 0
        	
        elif core_count > 0 and surface_count == 0:
            #everything is core (or buried, see above)
            self.rsa_change_score = 1
            
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

class RIN_based_features(Feature_set):
    __slots__ = [
        'profile', 'centralities', 'profile_str', 'centralities_str'
    ]

    def get_profile(self):
        if self.profile is not None:
            return self.profile
        if self.profile_str is None:
            return None
        return rin.Interaction_profile(profile_str=self.profile_str)

    def get_centralities(self):
        if self.centralities is not None:
            return self.centralities
        if self.centralities_str is None:
            return None
        return rin.Centrality_scores(code_str=self.centralities_str)

    def get_centralities_str(self):
        if self.centralities_str is not None:
            return self.centralities_str
        if self.centralities is None:
            return None
        return self.centralities.str_encode()

    def get_profile_str(self):
        if self.profile_str is not None:
            return self.profile_str
        if self.profile is None:
            return None
        return self.profile.encode()

    def get_raw_list(self):
        raw_list = [self.get_profile_str(), self.get_centralities_str()]
        return raw_list
    
    def set_values(self, raw_rin_based_features):
        if raw_rin_based_features is None:
            return
        self.profile_str, self.centrality_str = raw_rin_based_features
        self.centralities = rin.Centrality_scores(code_str=self.centrality_str)
        if self.profile_str is not None:
            self.profile = rin.Interaction_profile(profile_str=self.profile_str)

    def add_to_output_object(self, output_object):
        for feature_name in rin.Centrality_scores.feature_names:
            try:
                output_object.add_value(feature_name, getattr(self.centralities, feature_name))
            except:
                output_object.add_value(feature_name, None)

        for chaintype in ['mc', 'sc']:
            for interaction_type in ['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']:
                feature_name = '%s %s score' % (chaintype, interaction_type)
                if self.profile is None:
                    value = '-'
                else:
                    value = self.profile.getChainSpecificCombiScore(chaintype, interaction_type)
                output_object.add_value(feature_name, value)

                feature_name = '%s %s degree' % (chaintype, interaction_type)
                if self.profile is None:
                    value = '-'
                else:
                    value = self.profile.getChainSpecificCombiDegree(chaintype, interaction_type)
                output_object.add_value(feature_name, value)

                feature_name = '%s %s H-bond score' % (chaintype, interaction_type)
                if self.profile is None:
                    value = '-'
                else:
                    value = self.profile.getScore(chaintype, 'hbond', interaction_type)
                output_object.add_value(feature_name, value)

    def weight(self, unweighted_dicts, qualities):
        unweighted_centralities = unweighted_dicts['centralities']
        unweighted_profiles = unweighted_dicts['profile']
        self.centralities = None

        total_qual = 0.0

        for mapping_id in unweighted_centralities:
            centrality_scores = unweighted_centralities[mapping_id]
            if centrality_scores is None:
                continue
            qual = qualities[mapping_id]
            if self.centralities is None:
                self.centralities = [0.] * len(centrality_scores.cent_list)

            for pos, cent_score in enumerate(centrality_scores.cent_list):
                if cent_score is None:
                    continue
                self.centralities[pos] += cent_score * qual

            total_qual += qual
        if self.centralities is None:
            return

        if total_qual > 0.0:
            for pos, cent_score in enumerate(self.centralities):
                self.centralities[pos] = self.centralities[pos] / total_qual
            self.centralities = rin.Centrality_scores(cent_list=self.centralities)

        weight_profile_tuples = []
        for mapping_id in unweighted_profiles:
            weight_profile_tuples.append((qualities[mapping_id], unweighted_profiles[mapping_id]))
        self.profile = rin.calculateAverageProfile(weight_profile_tuples)

class Mappings:
    __slots__ = [
                    'qualities', 'covs', 'seq_ids',
                    'aa_ids', 'max_seq_res', 'recommended_res', 'interaction_recommendations', 'recommendation_order',
                    'resolutions', 'res_aas', 'amount_of_structures',
                    'structural_features',
                    'structural_unweighted_feature_dicts',
                    'rin_based_features',
                    'rin_based_unweighted_feature_dicts',
                    'integrated_features',
                    'microminer_features'
                ]


    def __init__(self, raw_results=None):
        self.qualities = {}
        self.seq_ids = {}
        self.covs = {}

        self.amount_of_structures = 0
        self.aa_ids = {}
        self.res_aas = {}
        self.resolutions = {}
        self.recommended_res = None
        self.recommendation_order = []
        self.max_seq_res = None

        self.interaction_recommendations = None

        self.structural_features = Structural_features()
        self.structural_unweighted_feature_dicts = self.structural_features.create_unweighted_feature_dicts()

        self.microminer_features = Microminer_features()

        self.rin_based_features = RIN_based_features()
        self.rin_based_unweighted_feature_dicts = self.rin_based_features.create_unweighted_feature_dicts()

        self.integrated_features = Integrated_features()

        if raw_results is not None:
            (self.recommendation_order, self.recommended_res, self.max_seq_res, self.interaction_recommendations, self.amount_of_structures,
            raw_structural_features, raw_microminer_features, raw_rin_features, raw_integrated_features
            ) = raw_results

            if isinstance(raw_structural_features, Structural_features):
                self.structural_features = raw_structural_features
            else:
                self.structural_features.set_values(raw_structural_features)

            if isinstance(raw_microminer_features, Microminer_features):
                self.microminer_features = raw_microminer_features
            else:
                self.microminer_features.set_values(raw_microminer_features)

            if isinstance(raw_rin_features, RIN_based_features):
                self.rin_based_features = raw_rin_features
            else:
                profile_str, centralities_str = raw_rin_features
                self.rin_based_features.profile_str = profile_str
                self.rin_based_features.centralities_str = centralities_str

            if isinstance(raw_integrated_features, Integrated_features):
                self.integrated_features = raw_integrated_features
            else:
                self.integrated_features.set_values(raw_integrated_features)

    def deconstruct(self):
        doomsday_protocol(self)

    def add_mapping(self, mapping_id, mapping):
        quality, seq_id, cov, structural_feature_dict, rin_based_feature_dict, identical_aa, resolution, res_aa = mapping
        if mapping_id not in self.qualities:
            self.amount_of_structures += 1
        self.qualities[mapping_id] = quality
        self.seq_ids[mapping_id] = seq_id
        self.covs[mapping_id] = cov
        for feature_name in structural_feature_dict:
            self.structural_unweighted_feature_dicts[feature_name][mapping_id] = structural_feature_dict[feature_name]
        for feature_name in rin_based_feature_dict:
            self.rin_based_unweighted_feature_dicts[feature_name][mapping_id] = rin_based_feature_dict[feature_name]
        
        self.aa_ids[mapping_id] = identical_aa
        self.res_aas[mapping_id] = res_aa
        self.resolutions[mapping_id[0]] = resolution

    def add_result(self, mapping_id, raw_results, quality, seq_id, cov):

        (recommendation_order, recommended_res, max_seq_res, interaction_recommendations, amount_of_structures, structural_features, microminer_features, rin_based_features, integrated_features) = raw_results

        if mapping_id not in self.qualities:
            self.amount_of_structures += 1
        self.qualities[mapping_id] = quality
        self.seq_ids[mapping_id] = seq_id
        self.covs[mapping_id] = cov

        structural_feature_values = structural_features.get_raw_list()
        for feat_pos, feat_name in enumerate(structural_features.get_feature_names()):
            self.structural_unweighted_feature_dicts[feat_name][mapping_id] = structural_feature_values[feat_pos]

        self.rin_based_unweighted_feature_dicts['profile'][mapping_id] = rin_based_features.get_profile()
        self.rin_based_unweighted_feature_dicts['centralities'][mapping_id] = rin_based_features.get_centralities()

        self.aa_ids[mapping_id] = True
        self.res_aas[mapping_id] = None

    def pack_features(self):
        packed_features = pack((
            str(self.interaction_recommendations),
            self.amount_of_structures,
            self.structural_features.get_raw_list(),
            self.rin_based_features.get_raw_list(),
            self.microminer_features.get_raw_list(),
            self.integrated_features.get_raw_list()
        ))
        return packed_features


    def get_raw_result(self):
        return (self.recommendation_order, self.recommended_res, self.max_seq_res, self.interaction_recommendations, self.amount_of_structures, self.structural_features, self.microminer_features, self.rin_based_features, self.integrated_features)

 
    def weight_all(self, config, disorder_score, disorder_region, for_indel_aggregation = False):

        for feature_name in self.structural_unweighted_feature_dicts:
            weighted_value = self.structural_features.weight_function(feature_name, self.structural_unweighted_feature_dicts[feature_name], self.qualities, self.covs, config)
            self.structural_features.set_value_by_name(feature_name, weighted_value)

        self.rin_based_features.weight(self.rin_based_unweighted_feature_dicts, self.qualities)
        self.integrated_features.generate(self.structural_features, config, self.rin_based_features.profile, disorder_score, disorder_region,
                                          self.structural_unweighted_feature_dicts['surface_value'],  self.structural_unweighted_feature_dicts['mainchain_surface_value'],  self.structural_unweighted_feature_dicts['sidechain_surface_value'])

        if not for_indel_aggregation:
            self.set_recommended_residues()  

        
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

            residue_simple_class = self.structural_unweighted_feature_dicts['rin_simple_class'][mapping_id]

            if self.integrated_features.simple_class == residue_simple_class:
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


