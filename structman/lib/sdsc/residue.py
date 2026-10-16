# sdsc: structman datastructures and classes
from structman.base_utils.base_utils import Slotted_obj
from structman.lib.sdsc.sdsc_utils import classify, doomsday_protocol
from structman.lib.rin import (
    Interaction_profile, Centrality_scores, Inflow_values, simple_class_codes_from_columns
)
from structman.lib.sdsc.consts.residues import THREE_TO_ONE

import numpy as np


modres_keys = sorted(THREE_TO_ONE.keys())
modres_dict = {x:p for p,x in enumerate(modres_keys)}

string_encodes = {
    3 : {'H': 0, ' ': 1, 'B': 2, 'E': 3, 'T': 4, 'S': 5, 'G': 6, 'I': 7},
    34: modres_dict,
    35 : {
    'I': 0,
    'V': 1,
    'L': 2,
    'F': 3,
    'C': 4,
    'M': 5,
    'A': 6,
    'G': 7,
    'T': 8,
    'S': 9,
    'W': 10,
    'Y': 11,
    'P': 12,
    'H': 13,
    'E': 14,
    'Q': 15,
    'D': 16,
    'N': 17,
    'K': 18,
    'R': 19,
    'X': 20,
    'B': 21,
    'Z': 22,
    'U': 23
}
}

string_decodes = {
    3 : ['H', ' ', 'B', 'E', 'T', 'S', 'G', 'I'],
    34 : modres_keys,
    35 : ['I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'S', 'W', 'Y', 'P', 'H', 'E', 'Q', 'D', 'N', 'K', 'R', 'X', 'B', 'Z', 'U']
}


class Residue(Slotted_obj):
    __slots__ = [
        'res_num', 'aa', 
        'RSA', 'relative_main_chain_acc',
        'relative_side_chain_acc', 'SSA',
        'interaction_profile', 'centralities', 'inflow_map',
        'modres', 'b_factor', 'pLDDT', 'database_id',
        'stored', 'phi', 'psi', 'tco', 'kappa', 'alpha',
        'intra_ssbond', 'inter_ssbond', 'ssbond_length',
        'intra_link', 'inter_link', 'link_length',
        'cis_conformation', 'cis_follower', 'inter_chain_median_kd',
        'inter_chain_dist_weighted_kd', 'inter_chain_median_rsa', 'inter_chain_dist_weighted_rsa',
        'intra_chain_median_kd', 'intra_chain_dist_weighted_kd', 'inter_chain_interactions_median',
        'inter_chain_interactions_dist_weighted', 'intra_chain_interactions_median', 'intra_chain_interactions_dist_weighted',
        'intra_chain_median_rsa', 'intra_chain_dist_weighted_rsa', 'Class',
        'simpleClass',
        'lig_dist', 'metal_dist', 'ion_dist', 'chain_dist', 'rna_dist', 'dna_dist',
        'intra_chain_ppi_degree', 'intra_chain_hom_degree',
        'inter_chain_ppi_degree', 'inter_chain_hom_degree'
        ]

    data_slots = [
        'RSA', 'relative_main_chain_acc', 'relative_side_chain_acc', 'SSA',
        'interaction_profile', 'centralities', 'inflow_map',
        'phi', 'psi', 'tco', 'kappa', 'alpha',
        'intra_ssbond', 'inter_ssbond', 'ssbond_length',
        'intra_link', 'inter_link', 'link_length', 'cis_conformation', 'cis_follower',
        'inter_chain_median_kd', 'inter_chain_dist_weighted_kd',
        'inter_chain_median_rsa', 'inter_chain_dist_weighted_rsa', 'intra_chain_median_kd',
        'intra_chain_dist_weighted_kd', 'intra_chain_median_rsa', 'intra_chain_dist_weighted_rsa',
        'inter_chain_interactions_median', 'inter_chain_interactions_dist_weighted',
        'intra_chain_interactions_median', 'intra_chain_interactions_dist_weighted',
        'b_factor', 'pLDDT', 'modres', 'aa',
        'lig_dist', 'metal_dist', 'ion_dist', 'chain_dist', 'rna_dist', 'dna_dist',
        'intra_chain_ppi_degree', 'intra_chain_hom_degree',
        'inter_chain_ppi_degree', 'inter_chain_hom_degree'
    ]

    values_per_slot = [
        1, 1, 1, 1,
        132, 12, 8,
        1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1,
        1, 1
    ]

    def __init__(self):
        pass

    def deconstruct(self):
        del self.res_num
        doomsday_protocol(self)

    def get_interaction_partners(self):
        if self.interaction_profile is None:
            return None, None
        return (self.interaction_profile.interacting_chains, self.interaction_profile.interacting_ligands)

    def get_rsa(self, splitted=False):
        if not splitted:
            return self.RSA
        else:
            return (self.RSA, self.relative_main_chain_acc, self.relative_side_chain_acc)

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
                config.errorlog.add_warning('Interaction profile is None, classification will fail')
            rin_class, rin_simple_class = classify(self.interaction_profile, sc)
            self.Class = rin_class
            self.simpleClass = rin_simple_class
        return self.Class, self.simpleClass

    def set_classification(self, Class, simpleClass):
        self.Class = Class
        self.simpleClass = simpleClass

    def checked_get_attribute(self, attribute_name: str):
        try:
            value = self.__getattribute__(attribute_name)
        except AttributeError:
            value = None
        return value

    def get_res_info(self):
        res_info = []
        for slot in self.data_slots:
            res_info.append(self.checked_get_attribute(slot))

        return res_info
    
    def set_res_info(self, res_data):
        if res_data is None:
            for slot in self.data_slots:
                self.__setattr__(slot, None)
            return
        for index, data in enumerate(res_data):
            self.__setattr__(self.data_slots[index], data)

    def restore_data_from_numpy(self, encoded: np.ndarray):
        left = 0
        for slot_index, nr_of_values in enumerate(self.values_per_slot):
            slot = self.data_slots[slot_index]
            if nr_of_values == 1:
                self.__setattr__(slot, encoded[left])
            else:
                if slot_index == 4:
                    #self.interaction_profile = Interaction_profile()
                    #self.interaction_profile.set_data_from_list(encoded[left:(left+nr_of_values)])
                    self.interaction_profile: np.ndarray = encoded[left:(left+nr_of_values)]
                elif slot_index == 5:
                    self.centralities = Centrality_scores()
                    self.centralities.set_data_from_list(encoded[left:(left+nr_of_values)])
                elif slot_index == 6:
                    self.inflow_map = Inflow_values()
                    self.inflow_map.set_data_from_list(encoded[left:(left+nr_of_values)])
            left += nr_of_values

    def encode_to_numpy(self) -> np.ndarray[np.float32]:
        # float32, not float64: nothing in the encoded row needs 64 bits.
        # Accessibilities, angles, distances and the small interaction counts all
        # carry far fewer significant digits than float32's 7, and the integer
        # codes for SSA/aa/modres are exact well past their largest value. Doing
        # it here rather than on read halves the stored blob as well as the
        # in-memory matrix, and makes the narrowing cast in retrieval a no-op for
        # newly written rows. Rows written before this change stay float64 and
        # are still narrowed on read, so both dtypes remain readable.
        raw_data = self.get_res_info()

        float_data: list[int | float] = []

        for slot_index, entry in enumerate(raw_data):
            # np.floating, not just float: np.float64 subclasses float but np.float32
            # does not, and every dist-weighted slot arrives as float32 - the
            # distances come out of exact_top20_search as a float32 array, and
            # `kd / dist` and the running weight total stay float32 from there under
            # NEP 50. Those slots therefore hit the TypeError below on any residue
            # with a populated neighbour sphere.
            if isinstance(entry, (int, float, np.integer, np.floating)):
                float_data.append(entry)
            elif isinstance(entry, str):
                try:
                    float_data.append(string_encodes[slot_index][entry])
                except (IndexError, KeyError) as err:
                    if slot_index == 34:
                        float_data.append(string_encodes[34]['UNK'])
                    else:
                        raise IndexError(f'{slot_index=} {self.data_slots[slot_index]=} {entry=} {err=}')
            elif isinstance(entry, np.ndarray):
                # An already-encoded row for one of the object slots below.
                # Interaction_profile.to_list walks a four level object tree to
                # rebuild the same 132 floats, and the aggregation path already has
                # them as a matrix row, so it hands that over instead.
                float_data += entry.tolist()
            elif isinstance(entry, Interaction_profile):
                float_data += entry.to_list()
            elif isinstance(entry, Centrality_scores):
                float_data += entry.to_list()
            elif isinstance(entry, Inflow_values):
                float_data += entry.to_list()
            elif entry is None:
                nr_of_nones = self.values_per_slot[slot_index]
                if nr_of_nones == 1:
                    float_data.append(None)
                else:
                    nones = [None] * nr_of_nones
                    float_data += nones
            else:
                raise TypeError(f'{slot_index=} {self.data_slots[slot_index]=} {entry=} {type(entry)=}')

        #print(f'{len(float_data)=}')
        return np.array(float_data, dtype=np.float32)


# The aggregation hot loop decodes millions of residues and only ever needs the
# get_res_info() list back out - the Residue instance itself is a throwaway
# intermediary. Round-tripping every data slot through __setattr__ and then
# __getattribute__ by name costs ~92 by-name attribute operations per residue,
# so the layout below lets the list be built straight off the encoded array.
# Offsets are derived from values_per_slot, so they cannot drift from the class.
def _build_res_info_layout() -> tuple[list[tuple[int, int]], dict[int, tuple[int, int]], int]:
    single_value_offsets: list[tuple[int, int]] = []
    multi_value_spans: dict[int, tuple[int, int]] = {}
    left = 0
    for slot_index, nr_of_values in enumerate(Residue.values_per_slot):
        if nr_of_values == 1:
            single_value_offsets.append((slot_index, left))
        else:
            multi_value_spans[slot_index] = (left, left + nr_of_values)
        left += nr_of_values
    return single_value_offsets, multi_value_spans, left


_RES_SINGLE_OFFSETS, _RES_MULTI_SPANS, _RES_ENCODED_LEN = _build_res_info_layout()
_RES_INFO_LEN = len(Residue.data_slots)

# Public view of the same layout, keyed by data slot name. The aggregation hot
# loop feeds whole encoded rows into column kernels instead of decoding every
# residue into python objects, so it needs to address the encoded array by
# feature name. Derived from the class definition, so it cannot drift either.
ENCODED_COLUMN_OF_SLOT: dict[str, int] = {
    Residue.data_slots[slot_index]: offset for slot_index, offset in _RES_SINGLE_OFFSETS
}
ENCODED_SPAN_OF_SLOT: dict[str, tuple[int, int]] = {
    Residue.data_slots[slot_index]: span for slot_index, span in _RES_MULTI_SPANS.items()
}
ENCODED_ROW_LEN: int = _RES_ENCODED_LEN

# Amino acid letter -> the code the encoded residue stores in its 'aa' slot.
# string_encodes is keyed by data slot index, so look that up rather than
# repeating the number.
AA_STRING_ENCODES: dict[str, int] = string_encodes[Residue.data_slots.index('aa')]

# The three names the aggregation loop needs individually, resolved once.
AA_COLUMN: int = ENCODED_COLUMN_OF_SLOT['aa']
PLDDT_COLUMN: int = ENCODED_COLUMN_OF_SLOT['pLDDT']
RSA_COLUMN: int = ENCODED_COLUMN_OF_SLOT['RSA']
MC_RSA_COLUMN: int = ENCODED_COLUMN_OF_SLOT['relative_main_chain_acc']
SC_RSA_COLUMN: int = ENCODED_COLUMN_OF_SLOT['relative_side_chain_acc']
PROFILE_SPAN: tuple[int, int] = ENCODED_SPAN_OF_SLOT['interaction_profile']
CENTRALITIES_SPAN: tuple[int, int] = ENCODED_SPAN_OF_SLOT['centralities']
INFLOW_SPAN: tuple[int, int] = ENCODED_SPAN_OF_SLOT['inflow_map']


def simple_class_codes_for_matrix(
        residue_matrix: np.ndarray, surface_threshold: float, buried_threshold: float
        ) -> np.ndarray[np.int8]:
    """The simple class of every residue of a chain, as the int codes rin defines.

    The simple class is not one of the data_slots, so it never reaches Residue_Data;
    every reader recomputes it from two parts of the encoded row. Storing it beside
    the matrix lets the alignment output read one byte per residue instead of the
    whole 195 float32 row.

    Lives here rather than in rin because it is the encoded row layout - PROFILE_SPAN
    and SC_RSA_COLUMN - that it needs, and rin cannot import this module back.
    """
    # An empty chain is packed as np.array([], dtype=np.float32), which is 1-D and
    # cannot be sliced by column.
    if residue_matrix.ndim < 2 or residue_matrix.shape[0] == 0:
        return np.empty(0, dtype=np.int8)

    profile_left, profile_right = PROFILE_SPAN
    return simple_class_codes_from_columns(
        residue_matrix[:, profile_left:profile_right],
        residue_matrix[:, SC_RSA_COLUMN],
        surface_threshold, buried_threshold)


def res_info_from_numpy(encoded: np.ndarray) -> list:
    # Equivalent to Residue().restore_data_from_numpy(encoded) followed by
    # get_res_info(), without materialising the Residue. Single-value slots stay
    # numpy scalars, exactly as the slot round-trip produced them.
    res_info: list = [None] * _RES_INFO_LEN

    for slot_index, offset in _RES_SINGLE_OFFSETS:
        res_info[slot_index] = encoded[offset]

    left, right = _RES_MULTI_SPANS[4]
    res_info[4] = encoded[left:right]

    left, right = _RES_MULTI_SPANS[5]
    centralities = Centrality_scores()
    centralities.set_data_from_list(encoded[left:right])
    res_info[5] = centralities

    left, right = _RES_MULTI_SPANS[6]
    inflow_map = Inflow_values()
    inflow_map.set_data_from_list(encoded[left:right])
    res_info[6] = inflow_map

    return res_info
