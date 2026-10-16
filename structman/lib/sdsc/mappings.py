import numpy as np
import time

# sdsc: structman datastructures and classes
from numba import njit
from structman.lib import rin
from structman.base_utils.base_utils import Slotted_obj, add_to_times
from structman.lib.sdsc.sdsc_utils import classify, triple_locate, doomsday_protocol, locate
from structman.lib.sdsc.residue import (
    ENCODED_COLUMN_OF_SLOT, CENTRALITIES_SPAN, INFLOW_SPAN, PROFILE_SPAN,
)

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

#@njit(cache=True)
def weight_array(
    matrix: np.ndarray,
    qualities: np.ndarray,
    subset_threshold = None,
    na_value = np.nan
    ) -> np.ndarray:
    
    nom_vec: np.ndarray = matrix * qualities

    if subset_threshold is None:
        weighted_val = np.nansum(nom_vec) / np.sum(qualities)
    else:
        indeces = np.flatnonzero(nom_vec > subset_threshold)
        if len(indeces) == 0:
            return na_value
        sub_vector = nom_vec[indeces]

        sub_quals = qualities[indeces]

        weighted_val = np.sum(sub_vector, axis=0) / np.sum(sub_quals)

    return weighted_val

@njit(cache=True)
def _count_surface_core(values: np.ndarray, surface_threshold: float) -> tuple[int, int]:
    # locate() only distinguishes Surface from {Core, Buried} for these counts, so the
    # buried threshold never matters here. NaN fails the comparison and lands in the
    # non-surface bucket, which is exactly what locate() does with it ("Core").
    surface_count: int = 0
    non_surface_count: int = 0

    for index in range(values.shape[0]):
        if values[index] > surface_threshold:
            surface_count += 1
        else:
            non_surface_count += 1

    return surface_count, non_surface_count


def location_counts(values, surface_threshold: float) -> tuple[int, int]:
    values_array = np.asarray(values, dtype=np.float64)

    surface_count, non_surface_count = _count_surface_core(values_array, surface_threshold)

    # locate() returns None for a None rsa, which the original loops counted as neither
    # surface nor core. numpy collapses None into NaN on conversion, so those entries
    # land in the non-surface bucket above and have to be taken back out.
    if isinstance(values, list):
        non_surface_count -= values.count(None)

    return surface_count, non_surface_count


def not_none_float_array(values) -> np.ndarray:
    # Explicitly float64: the encoded residue data is float32, and letting np.std
    # accumulate in float32 both loses precision against the previous behaviour and
    # produces an np.float32 that msgpack cannot serialise.
    if isinstance(values, list) and values.count(None) > 0:
        return np.array([value for value in values if value is not None], dtype=np.float64)
    return np.asarray(values, dtype=np.float64)


@njit(cache=True)
def weight_profile_columns(
    profiles: np.ndarray,
    qualities: np.ndarray,
    subset_threshold: float,
    na_value: float,
    ) -> np.ndarray:
    # Column-wise equivalent of weight_array(..., subset_threshold, na_value) over a
    # whole (mappings x features) profile matrix. The per-column Python loop this
    # replaces made ~6 numpy calls per column for 132 columns per position, where the
    # per-call overhead dwarfed the actual arithmetic on these short vectors.
    n_rows, n_cols = profiles.shape
    weighted = np.empty(n_cols, dtype=np.float64)

    for col in range(n_cols):
        nom: float = 0.0
        denom: float = 0.0
        any_above_threshold: bool = False

        for row in range(n_rows):
            # NaNs compare False here, matching the np.flatnonzero(nom_vec > thresh) filter.
            value = profiles[row, col] * qualities[row]
            if value > subset_threshold:
                nom += value
                denom += qualities[row]
                any_above_threshold = True

        if any_above_threshold:
            weighted[col] = nom / denom
        else:
            weighted[col] = na_value

    return weighted


@njit(cache=True)
def weight_profile_columns_span(
    matrix: np.ndarray,
    left: int,
    right: int,
    qualities: np.ndarray,
    subset_threshold: float,
    na_value: float,
    ) -> np.ndarray:
    # weight_profile_columns over a column span of the encoded-row matrix, so the
    # profile never has to be sliced out into its own array.
    n_rows = matrix.shape[0]
    weighted = np.empty(right - left, dtype=np.float64)

    for col in range(left, right):
        nom: float = 0.0
        denom: float = 0.0
        any_above_threshold: bool = False

        for row in range(n_rows):
            value = matrix[row, col] * qualities[row]
            if value > subset_threshold:
                nom += value
                denom += qualities[row]
                any_above_threshold = True

        if any_above_threshold:
            weighted[col - left] = nom / denom
        else:
            weighted[col - left] = na_value

    return weighted


@njit(cache=True)
def weight_span_shared_denominator(
    matrix: np.ndarray,
    left: int,
    right: int,
    qualities: np.ndarray,
    ) -> tuple:
    # Mirrors the centrality / inflow aggregation in RIN_based_features.weight:
    # NaN entries are skipped per value, but the denominator is the total quality
    # over all rows regardless, so a NaN in one feature does not change the
    # divisor of the others.
    n_rows = matrix.shape[0]
    accumulated = np.zeros(right - left, dtype=np.float64)
    total_qual: float = 0.0

    for row in range(n_rows):
        qual = qualities[row]
        for col in range(left, right):
            value = matrix[row, col]
            if value != value:  # NaN
                continue
            accumulated[col - left] += value * qual
        total_qual += qual

    return accumulated, total_qual


@njit(cache=True)
def weight_columns_float(
    matrix: np.ndarray,
    columns: np.ndarray,
    qualities: np.ndarray,
    ) -> np.ndarray:
    # Column-wise weight_float_no_conf. One call replaces one python-level call
    # plus one np.array(list) conversion per feature.
    n_rows = matrix.shape[0]
    out = np.empty(columns.shape[0], dtype=np.float64)

    for i in range(columns.shape[0]):
        col = columns[i]
        nom: float = 0.0
        denom: float = 0.0
        for row in range(n_rows):
            value = matrix[row, col]
            if value != value:
                continue
            weight = qualities[row]
            nom += value * weight
            denom += weight
        if denom > 0.0:
            out[i] = nom / denom
        else:
            out[i] = np.nan

    return out


@njit(cache=True)
def weight_columns_dist(
    matrix: np.ndarray,
    columns: np.ndarray,
    qualities: np.ndarray,
    threshold: float,
    ) -> np.ndarray:
    # Column-wise weight_dist_float_no_conf.
    n_rows = matrix.shape[0]
    out = np.empty(columns.shape[0], dtype=np.float64)

    for i in range(columns.shape[0]):
        col = columns[i]
        nom: float = 0.0
        denom: float = 0.0
        for row in range(n_rows):
            value = matrix[row, col]
            if value != value:
                continue
            if value > threshold:
                continue
            weight = qualities[row]
            nom += value * weight
            denom += weight
        if denom > 0.0:
            out[i] = nom / denom
        else:
            out[i] = np.nan

    return out


@njit(cache=True)
def weight_columns_cov(
    matrix: np.ndarray,
    columns: np.ndarray,
    qualities: np.ndarray,
    covs: np.ndarray,
    coverage_extra_weight: float,
    ) -> np.ndarray:
    # Column-wise weight_float_no_conf_cov_weight.
    n_rows = matrix.shape[0]
    out = np.empty(columns.shape[0], dtype=np.float64)

    for i in range(columns.shape[0]):
        col = columns[i]
        nom: float = 0.0
        denom: float = 0.0
        for row in range(n_rows):
            value = matrix[row, col]
            if value != value:
                continue
            weight = qualities[row] * (covs[row] ** coverage_extra_weight)
            nom += value * weight
            denom += weight
        if denom > 0.0:
            out[i] = nom / denom
        else:
            out[i] = np.nan

    return out


@njit(cache=True)
def weight_columns_propensity(
    matrix: np.ndarray,
    columns: np.ndarray,
    ) -> np.ndarray:
    # Column-wise weight_propensity: the fraction of mappings that actually
    # carry the feature. These slots hold a measured angle that is only set for
    # residues in a cis conformation, so a missing one encodes as NaN.
    n_rows = matrix.shape[0]
    out = np.empty(columns.shape[0], dtype=np.float64)

    for i in range(columns.shape[0]):
        col = columns[i]
        present: float = 0.0
        for row in range(n_rows):
            value = matrix[row, col]
            if value != value:
                continue
            present += 1.0
        out[i] = present / n_rows

    return out


@njit(cache=True)
def weight_columns_bool_propensity(
    matrix: np.ndarray,
    columns: np.ndarray,
    ) -> np.ndarray:
    # Column-wise weigthed_bool_propensity. That function skips None and counts
    # everything truthy, and an encoded value is never None, so the test is
    # 'not zero' - which NaN passes, exactly as `if value:` did.
    n_rows = matrix.shape[0]
    out = np.empty(columns.shape[0], dtype=np.float64)

    for i in range(columns.shape[0]):
        col = columns[i]
        true_prop: float = 0.0
        for row in range(n_rows):
            if matrix[row, col] != 0.0:
                true_prop += 1.0
        out[i] = true_prop / n_rows

    return out


@njit(cache=True)
def weight_columns_majority(
    matrix: np.ndarray,
    columns: np.ndarray,
    qualities: np.ndarray,
    ) -> np.ndarray:
    # Column-wise weight_majority. The candidate values are small integer codes
    # (SSA, modres), so a linear scan over the distinct values seen so far is
    # cheap and - unlike bucketing by code - preserves the original dict's
    # first-seen ordering, which decides ties under the strict > below.
    n_rows = matrix.shape[0]
    out = np.empty(columns.shape[0], dtype=np.float64)

    seen = np.empty(n_rows, dtype=np.float64)
    weights = np.empty(n_rows, dtype=np.float64)

    for i in range(columns.shape[0]):
        col = columns[i]
        n_seen = 0

        for row in range(n_rows):
            value = matrix[row, col]
            if value != value:
                continue
            qual = qualities[row]
            found = -1
            for s in range(n_seen):
                if seen[s] == value:
                    found = s
                    break
            if found < 0:
                seen[n_seen] = value
                weights[n_seen] = qual
                n_seen += 1
            else:
                weights[found] += qual

        best_value = np.nan
        max_qual: float = 0.0
        for s in range(n_seen):
            if weights[s] > max_qual:
                max_qual = weights[s]
                best_value = seen[s]
        out[i] = best_value

    return out


# --- streaming forms of the weighting kernels ---------------------------------------
#
# A classification task used to decode every structure one protein maps to before the
# position loop started, so it held the whole residue matrix - up to 8.6 GiB on a
# whole-PDB protein, which is what capped the stage at ~70 of 255 cores. Batching the
# annotation list means each position is visited once per batch, so every reduction
# weight_all performs has to accumulate across batches instead of seeing all the
# mappings at once.
#
# Each kernel below mirrors its one-shot twin exactly, statement for statement, and
# carries the state the twin kept in locals. The twins stay: the unbatched path still
# uses them, and dev_tests/test_classification_batching.py pins the two against each
# other. Where the one-shot divides by its own row count, the streaming finalize divides
# by the total across batches, which is the same number - and where the one-shot would
# divide by zero, so does the finalize, deliberately, rather than papering over a case
# that would otherwise change behaviour.

# Distinct values a majority column can hold before the state cannot track first-seen
# order any more. The categorical features are ssa and modres - a secondary structure
# code and a modified-residue flag - so the real count is under ten; this is only here so
# that an unexpected encoding fails loudly instead of silently dropping candidates.
_MAX_MAJORITY_VALUES: int = 64


@njit(cache=True)
def accumulate_columns_float(
    matrix: np.ndarray,
    columns: np.ndarray,
    qualities: np.ndarray,
    nom: np.ndarray,
    denom: np.ndarray,
    ) -> None:
    # Streaming weight_columns_float.
    for i in range(columns.shape[0]):
        col = columns[i]
        for row in range(matrix.shape[0]):
            value = matrix[row, col]
            if value != value:
                continue
            weight = qualities[row]
            nom[i] += value * weight
            denom[i] += weight


@njit(cache=True)
def accumulate_columns_dist(
    matrix: np.ndarray,
    columns: np.ndarray,
    qualities: np.ndarray,
    threshold: float,
    nom: np.ndarray,
    denom: np.ndarray,
    ) -> None:
    # Streaming weight_columns_dist.
    for i in range(columns.shape[0]):
        col = columns[i]
        for row in range(matrix.shape[0]):
            value = matrix[row, col]
            if value != value:
                continue
            if value > threshold:
                continue
            weight = qualities[row]
            nom[i] += value * weight
            denom[i] += weight


@njit(cache=True)
def accumulate_columns_cov(
    matrix: np.ndarray,
    columns: np.ndarray,
    qualities: np.ndarray,
    covs: np.ndarray,
    coverage_extra_weight: float,
    nom: np.ndarray,
    denom: np.ndarray,
    ) -> None:
    # Streaming weight_columns_cov.
    for i in range(columns.shape[0]):
        col = columns[i]
        for row in range(matrix.shape[0]):
            value = matrix[row, col]
            if value != value:
                continue
            weight = qualities[row] * (covs[row] ** coverage_extra_weight)
            nom[i] += value * weight
            denom[i] += weight


@njit(cache=True)
def finalize_columns_nom_denom(nom: np.ndarray, denom: np.ndarray) -> np.ndarray:
    # The tail of weight_columns_float / _dist / _cov, which share it.
    out = np.empty(nom.shape[0], dtype=np.float64)
    for i in range(nom.shape[0]):
        if denom[i] > 0.0:
            out[i] = nom[i] / denom[i]
        else:
            out[i] = np.nan
    return out


@njit(cache=True)
def accumulate_columns_propensity(
    matrix: np.ndarray,
    columns: np.ndarray,
    present: np.ndarray,
    ) -> None:
    # Streaming weight_columns_propensity; the row total is counted by the caller, since
    # it is shared with the boolean columns.
    for i in range(columns.shape[0]):
        col = columns[i]
        for row in range(matrix.shape[0]):
            value = matrix[row, col]
            if value != value:
                continue
            present[i] += 1.0


@njit(cache=True)
def accumulate_columns_bool_propensity(
    matrix: np.ndarray,
    columns: np.ndarray,
    true_prop: np.ndarray,
    ) -> None:
    # Streaming weight_columns_bool_propensity. NaN is not zero, so it counts as true -
    # exactly as `if value:` did.
    for i in range(columns.shape[0]):
        col = columns[i]
        for row in range(matrix.shape[0]):
            if matrix[row, col] != 0.0:
                true_prop[i] += 1.0


@njit(cache=True)
def finalize_columns_over_rows(counts: np.ndarray, n_rows: int) -> np.ndarray:
    # The tail of both propensity kernels: a count over the row total. Divides by zero
    # for an empty matrix, which is what the one-shot form does too.
    out = np.empty(counts.shape[0], dtype=np.float64)
    for i in range(counts.shape[0]):
        out[i] = counts[i] / n_rows
    return out


@njit(cache=True)
def accumulate_columns_majority(
    matrix: np.ndarray,
    columns: np.ndarray,
    qualities: np.ndarray,
    seen: np.ndarray,
    weights: np.ndarray,
    n_seen: np.ndarray,
    ) -> int:
    # Streaming weight_columns_majority. The one-shot keeps its candidates in
    # first-seen order because the strict > in the finalize resolves ties to the earlier
    # one; carrying that order across batches is what makes the batched result identical,
    # and it is why batches have to be walked in annotation order.
    #
    # Returns how many values did not fit, so the caller can fail loudly rather than
    # quietly dropping a candidate and changing the winner.
    overflow = 0
    max_values = seen.shape[1]

    for i in range(columns.shape[0]):
        col = columns[i]
        for row in range(matrix.shape[0]):
            value = matrix[row, col]
            if value != value:
                continue
            qual = qualities[row]
            found = -1
            for s in range(n_seen[i]):
                if seen[i, s] == value:
                    found = s
                    break
            if found < 0:
                if n_seen[i] >= max_values:
                    overflow += 1
                    continue
                seen[i, n_seen[i]] = value
                weights[i, n_seen[i]] = qual
                n_seen[i] += 1
            else:
                weights[i, found] += qual

    return overflow


@njit(cache=True)
def finalize_columns_majority(
    seen: np.ndarray,
    weights: np.ndarray,
    n_seen: np.ndarray,
    ) -> np.ndarray:
    # The tail of weight_columns_majority: the heaviest candidate, ties to the earlier.
    out = np.empty(seen.shape[0], dtype=np.float64)
    for i in range(seen.shape[0]):
        best_value = np.nan
        max_qual: float = 0.0
        for s in range(n_seen[i]):
            if weights[i, s] > max_qual:
                max_qual = weights[i, s]
                best_value = seen[i, s]
        out[i] = best_value
    return out


@njit(cache=True)
def accumulate_span_shared_denominator(
    matrix: np.ndarray,
    left: int,
    right: int,
    qualities: np.ndarray,
    accumulated: np.ndarray,
    total_qual: np.ndarray,
    ) -> None:
    # Streaming weight_span_shared_denominator. total_qual is a one-element array so a
    # nopython function can carry it.
    for row in range(matrix.shape[0]):
        qual = qualities[row]
        for col in range(left, right):
            value = matrix[row, col]
            if value != value:
                continue
            accumulated[col - left] += value * qual
        total_qual[0] += qual


@njit(cache=True)
def accumulate_profile_columns_span(
    matrix: np.ndarray,
    left: int,
    right: int,
    qualities: np.ndarray,
    subset_threshold: float,
    nom: np.ndarray,
    denom: np.ndarray,
    any_above: np.ndarray,
    ) -> None:
    # Streaming weight_profile_columns_span.
    for col in range(left, right):
        index = col - left
        for row in range(matrix.shape[0]):
            value = matrix[row, col] * qualities[row]
            if value > subset_threshold:
                nom[index] += value
                denom[index] += qualities[row]
                any_above[index] = 1


@njit(cache=True)
def finalize_profile_columns_span(
    nom: np.ndarray,
    denom: np.ndarray,
    any_above: np.ndarray,
    na_value: float,
    ) -> np.ndarray:
    weighted = np.empty(nom.shape[0], dtype=np.float64)
    for index in range(nom.shape[0]):
        if any_above[index] != 0:
            weighted[index] = nom[index] / denom[index]
        else:
            weighted[index] = na_value
    return weighted


@njit(cache=True)
def weight_float_no_conf(
    values: np.ndarray,
    qualities: np.ndarray,
    ) -> float:
    nom: float = 0.0
    denom: float = 0.0

    for mapping_id, value in enumerate(values):
        if np.isnan(value):
            continue
        
        weight: float = qualities[mapping_id]
        nom += value * weight
        denom += weight

    if denom > 0.0:
        weighted_value = nom / denom
        return weighted_value
    else:
        return np.nan
    
@njit(cache=True)
def weight_dist_float_no_conf(
    values: np.ndarray,
    qualities: np.ndarray,
    threshold: float,
    ) -> float:
    nom: float = 0.0
    denom: float = 0.0

    for mapping_id, value in enumerate(values):
        if np.isnan(value):
            continue
        if value > threshold:
            continue
        weight: float = qualities[mapping_id]
        nom += value * weight
        denom += weight

    if denom > 0.0:
        weighted_value = nom / denom
        return weighted_value
    else:
        return np.nan

@njit(cache=True)
def weight_float_no_conf_cov_weight(
    values: np.ndarray,
    qualities: np.ndarray,
    covs: np.ndarray,
    coverage_extra_weight: float,
    ) -> float:
    nom: float = 0.0
    denom: float = 0.0

    for mapping_id, value in enumerate(values):
        if np.isnan(value):
            continue

        qual: float = qualities[mapping_id]
        cov: float = covs[mapping_id]
        weight: float = qual * (cov ** coverage_extra_weight)
        nom += value * weight
        denom += weight

    if denom > 0.0:
        weighted_value = nom / denom
        return weighted_value
    else:
        return np.nan

def weight(
        value_map: list,
        config,
        qualities: list[float],
        covs:list[float],
        distance_weighting=False,
        calc_conf=True,
        coverage_extra_weight=0):
    
    if not distance_weighting and not calc_conf and coverage_extra_weight == 0:
        return weight_float_no_conf(np.array(value_map), np.array(qualities))
    
    if not distance_weighting and not calc_conf and coverage_extra_weight != 0:
        return weight_float_no_conf_cov_weight(np.array(value_map), np.array(qualities), np.array(covs), coverage_extra_weight)

    if distance_weighting and not calc_conf and coverage_extra_weight == 0:
        return weight_dist_float_no_conf(np.array(value_map), np.array(qualities), config.long_distance_threshold)

    nom = 0.0
    denom = 0.0
    n = 0.0
    qs = []
    weighted_value = None
    conf = 0.0
    for mapping_id, value in enumerate(value_map):
        if value is None or np.isnan(value):
            continue
        if distance_weighting and value > config.long_distance_threshold:
            continue
        if not isinstance(value, float):
            raise ValueError(f'Strange value: {value=}')
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
    for value in value_map:
        if value is None:
            continue
        # A missing value survives the numpy encoding as NaN rather than None,
        # so testing only for None counted every mapping and made this constant
        # 1.0. These slots hold an angle that is only set for residues in a cis
        # conformation, so NaN means absent and must not count.
        if value != value:
            continue
        prop += 1.
    weighted_prop = prop / float(len(value_map))
    return weighted_prop


def weigthed_bool_propensity(value_map):
    true_prop = 0.

    if len(value_map) == 0:
        return None
    for value in value_map:
        if value is None:
            continue
        if value:
            true_prop += 1.

    weighted_true_prop = true_prop / float(len(value_map))
    return weighted_true_prop


def weight_majority(value_map, qualities):
    voting = {}
    best_value = None
    for mapping_id, value in enumerate(value_map):
        if value is None or np.isnan(value):
            continue
        qual = qualities[mapping_id]
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

class Feature_set(Slotted_obj):
    def __init__(self, init_none = False):
        if init_none:
            for feature_name in self.__slots__:
                setattr(self, feature_name, None)
        return
    
    """
    def __serialize__(self) -> list[any]:
        serialized_obj: list[any] = []
        for attribute_name in self.__slots__:
            serialized_obj.append(self.__getattribute__(attribute_name))

        return serialized_obj
    """

    def set_values(self, values):
        if values is None:
            return
        for feat_pos, feature_name in enumerate(self.__slots__):
            setattr(self, feature_name, values[feat_pos])

    def set_value_by_name(self, feature_name, value):
        setattr(self, feature_name, value)
    
    def get_feature_names(self):
        return self.__slots__

    def add_to_output_object(self, output_object):
        for feature_name in self.__slots__:
            output_object.add_value(feature_name, getattr(self, feature_name))

    def create_unweighted_feature_dicts(self) -> dict[str, list]:
        unweighted_feature_dicts: dict[str, list] = {}
        for feature_name in self.__slots__:
            unweighted_feature_dicts[feature_name] = []
        return unweighted_feature_dicts
    
    def create_unweighted_feature_lists(self) -> dict[str, list]:
        unweighted_feature_dicts: list[list] = []
        for _ in self.__slots__:
            unweighted_feature_dicts.append([])
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
        'b_factor', 'pLDDT', 'modres', 'ssa', 'phi', 'psi', 'tco', 'kappa', 'alpha',
        'intra_ssbond', 'inter_ssbond', 'ssbond_length', 'intra_link', 'inter_link', 'link_length',
        'cis_conformation', 'cis_follower', 'inter_chain_median_kd', 'inter_chain_dist_weighted_kd',
        'inter_chain_median_rsa', 'inter_chain_dist_weighted_rsa', 'intra_chain_median_kd', 'intra_chain_dist_weighted_kd', 'intra_chain_median_rsa', 'intra_chain_dist_weighted_rsa',
        'inter_chain_interactions_median', 'inter_chain_interactions_dist_weighted', 'intra_chain_interactions_median', 'intra_chain_interactions_dist_weighted',
        'lig_dist', 'metal_dist', 'ion_dist', 'chain_dist', 'rna_dist', 'dna_dist', 'surface_value', 'mainchain_surface_value', 'sidechain_surface_value',
        'rin_class', 'rin_simple_class',
        'intra_chain_ppi_degree', 'intra_chain_hom_degree',
        'inter_chain_ppi_degree', 'inter_chain_hom_degree'
    ]

    distance_based_features = set([
        'lig_dist', 'metal_dist', 'ion_dist', 'chain_dist', 'rna_dist', 'dna_dist'
    ])

    float_features = set([
        'b_factor', 'pLDDT', 'phi', 'psi', 'tco', 'kappa', 'alpha',
        'ssbond_length', 'link_length', 'inter_chain_median_kd', 'inter_chain_dist_weighted_kd',
        'inter_chain_median_rsa', 'inter_chain_dist_weighted_rsa', 'intra_chain_median_kd', 'intra_chain_dist_weighted_kd', 'intra_chain_median_rsa', 'intra_chain_dist_weighted_rsa',
        'inter_chain_interactions_median', 'inter_chain_interactions_dist_weighted', 'intra_chain_interactions_median', 'intra_chain_interactions_dist_weighted',
        'intra_chain_ppi_degree', 'intra_chain_hom_degree',
        'inter_chain_ppi_degree', 'inter_chain_hom_degree'
    ])

    probability_features = set([
        'cis_conformation', 'cis_follower'
    ])

    boolean_features = set([
        'intra_ssbond', 'inter_ssbond', 'intra_link', 'inter_link',
    ])

    categorical_features = set([
        'ssa', 'modres'
    ])

    surface_value_features = set([
        'surface_value', 'mainchain_surface_value', 'sidechain_surface_value'
    ])

    def weight_function(self, feature_name: str, value_map: list, qualities, covs, config):
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

    def set_values_from_encoded(self, matrix: np.ndarray, qualities: np.ndarray, covs: np.ndarray, config) -> None:
        # Matrix equivalent of calling weight_function() once per feature. Every
        # numeric structural feature is a single column of the encoded residue
        # row, so the whole feature set is six kernel calls over one matrix
        # instead of 44 python calls each preceded by an np.array(list).
        layout = _STRUCTURAL_COLUMN_LAYOUT

        weighted = weight_columns_dist(matrix, layout.distance_columns, qualities, config.long_distance_threshold)
        for i, feature_name in enumerate(layout.distance_names):
            self.__setattr__(feature_name, weighted[i])

        weighted = weight_columns_float(matrix, layout.float_columns, qualities)
        for i, feature_name in enumerate(layout.float_names):
            self.__setattr__(feature_name, weighted[i])

        weighted = weight_columns_cov(matrix, layout.surface_columns, qualities, covs, 2)
        for i, feature_name in enumerate(layout.surface_names):
            self.__setattr__(feature_name, weighted[i])

        weighted = weight_columns_bool_propensity(matrix, layout.boolean_columns)
        for i, feature_name in enumerate(layout.boolean_names):
            self.__setattr__(feature_name, weighted[i])

        weighted = weight_columns_majority(matrix, layout.categorical_columns, qualities)
        for i, feature_name in enumerate(layout.categorical_names):
            value = weighted[i]
            # weight_majority returns None when nothing outvoted a zero quality.
            self.__setattr__(feature_name, None if value != value else value)

        if matrix.shape[0] > 0:
            weighted = weight_columns_propensity(matrix, layout.probability_columns)
            for i, feature_name in enumerate(layout.probability_names):
                self.__setattr__(feature_name, weighted[i])
        else:
            for feature_name in layout.probability_names:
                self.__setattr__(feature_name, None)

        # weight_function returns None for anything in no category.
        for feature_name in layout.uncategorised_names:
            self.__setattr__(feature_name, None)


class _StructuralColumnLayout:
    """Encoded-column indices for each weighting category, resolved once."""

    __slots__ = [
        'distance_names', 'distance_columns', 'float_names', 'float_columns',
        'surface_names', 'surface_columns', 'boolean_names', 'boolean_columns',
        'categorical_names', 'categorical_columns', 'probability_names',
        'probability_columns', 'uncategorised_names',
    ]

    # Structural feature names that do not match their Residue.data_slots name.
    NAME_ALIASES = {
        'ssa': 'SSA',
        'surface_value': 'RSA',
        'mainchain_surface_value': 'relative_main_chain_acc',
        'sidechain_surface_value': 'relative_side_chain_acc',
    }

    def __init__(self):
        def columns_for(feature_names):
            names = [n for n in Structural_features.__slots__ if n in feature_names]
            columns = np.array(
                [ENCODED_COLUMN_OF_SLOT[self.NAME_ALIASES.get(n, n)] for n in names],
                dtype=np.int64)
            return names, columns

        self.distance_names, self.distance_columns = columns_for(Structural_features.distance_based_features)
        self.float_names, self.float_columns = columns_for(Structural_features.float_features)
        self.surface_names, self.surface_columns = columns_for(Structural_features.surface_value_features)
        self.boolean_names, self.boolean_columns = columns_for(Structural_features.boolean_features)
        self.categorical_names, self.categorical_columns = columns_for(Structural_features.categorical_features)
        self.probability_names, self.probability_columns = columns_for(Structural_features.probability_features)

        categorised = set()
        for names in (self.distance_names, self.float_names, self.surface_names,
                      self.boolean_names, self.categorical_names, self.probability_names):
            categorised.update(names)
        self.uncategorised_names = [n for n in Structural_features.__slots__ if n not in categorised]


_STRUCTURAL_COLUMN_LAYOUT = _StructuralColumnLayout()


class Integrated_features(Feature_set):
    __slots__ = [
        'location', 'mainchain_location', 'sidechain_location',
        'structural_classification', 'simple_class',
        'rsa_change_score', 'mc_rsa_change_score', 'sc_rsa_change_score', 'sc_rsa_std'
    ]

    def generate(self, structural_features, config, profile, surface_values, mc_surface_values, sc_surface_values, encoded_profile=None):
        (self.location, self.mainchain_location, self.sidechain_location) = triple_locate(structural_features.surface_value, structural_features.mainchain_surface_value, structural_features.sidechain_surface_value, config)

        if encoded_profile is None:
            # profile is the stored flat list now, and compute_class_from_encoded only
            # reads scores, which the stored form carries unchanged. classify() is kept
            # for the callers that still hand in an Interaction_profile object.
            self.structural_classification, self.simple_class = rin.classify_encoded(profile, self.sidechain_location)
        else:
            self.structural_classification, self.simple_class = rin.classify_encoded(encoded_profile, self.sidechain_location)

        #analyse for conformational changes
        # These three counts used to be python loops calling locate() per mapping, which
        # is thousands of calls per position - each returning a string that was then
        # compared against three literals - for what is a single threshold test.
        self.rsa_change_score = None
        self.mc_rsa_change_score = None
        self.sc_rsa_change_score = None

        surface_threshold: float = config.surface_threshold

        surface_count, core_count = location_counts(surface_values, surface_threshold)
        mc_surface_count, mc_core_count = location_counts(mc_surface_values, surface_threshold)
        sc_surface_count, sc_core_count = location_counts(sc_surface_values, surface_threshold)

        sc_surface_list = not_none_float_array(sc_surface_values)

        if len(sc_surface_list) > 1:
            self.sc_rsa_std = np.std(sc_surface_list)
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
        'profile', 'centralities', 'inflow_values'
    ]
    
    def add_to_output_object(self, output_object):
        #print(f'{self.centralities=} {self.inflow_values}')
        if isinstance(self.centralities, rin.Centrality_scores):
            for feature_name in rin.Centrality_scores.__slots__:
                try:
                    output_object.add_value(feature_name, getattr(self.centralities, feature_name))
                except:
                    output_object.add_value(feature_name, None)
        else:
            for index, feature_name in enumerate(rin.Centrality_scores.__slots__):
                output_object.add_value(feature_name, None)

        for feature_name in rin.Inflow_values.__slots__:
            try:
                output_object.add_value(feature_name, getattr(self.inflow_values, feature_name))
            except:
                output_object.add_value(feature_name, None)

        # The profile is the flat list rin.stored_profile_from_encoded produces, so
        # the 66 values below are indexed straight out of it. They used to walk four
        # levels of string dispatch into the object tree per value, which was 86% of
        # this method and ran once per SNV.
        if self.profile is None:
            for feature_name in rin.PROFILE_OUTPUT_FEATURE_NAMES:
                output_object.add_value(feature_name, '-')
        else:
            # as_stored_profile, because a row packed by an older version still
            # decodes into the object tree.
            profile = rin.as_stored_profile(self.profile)
            for (score_name, degree_name, hbond_name,
                 score_indices, degree_indices, hbond_index) in rin.PROFILE_OUTPUT_LAYOUT:

                # getCombiScore nansummed these three, which is a plain sum with the
                # NaNs dropped - and 0.0 when all three are missing, as nansum gives.
                combi_score = 0.0
                for index in score_indices:
                    value = profile[index]
                    if value == value:
                        combi_score += value
                output_object.add_value(score_name, combi_score)

                # Degrees were truncated to int on the way in, so this stays an int
                # sum and formats without a decimal point, the way it always did.
                output_object.add_value(
                    degree_name,
                    profile[degree_indices[0]] + profile[degree_indices[1]] + profile[degree_indices[2]])

                output_object.add_value(hbond_name, profile[hbond_index])

    def weight_from_encoded(self, matrix: np.ndarray, qualities: np.ndarray):
        # Matrix equivalent of weight(). The centralities and inflow values live
        # in the encoded row already, so aggregating them straight out of it
        # avoids materialising a Centrality_scores / Inflow_values object per
        # mapping only to read every slot back by name.
        times = []
        ti = time.time()

        self.centralities = None
        if matrix.shape[0] > 0:
            left, right = CENTRALITIES_SPAN
            accumulated, total_qual = weight_span_shared_denominator(matrix, left, right, qualities)
            self.centralities = list(accumulated)
            if total_qual > 0.0:
                agg_cent_obj = rin.Centrality_scores()
                for pos, feat_name in enumerate(agg_cent_obj.__slots__):
                    self.centralities[pos] = self.centralities[pos] / total_qual
                    agg_cent_obj.__setattr__(feat_name, self.centralities[pos])
                self.centralities = agg_cent_obj

        if self.centralities is None:
            # Same early exit as weight(): no mapping contributed centralities,
            # so the profile stays unset and weight_all raises AttributeError.
            # Still has to match the (times, encoded_profile) contract of the
            # normal exit, since weight_all unpacks the result unconditionally.
            return times, None

        ti = add_to_times(times, ti, 'aggregate centralities')

        self.inflow_values = None
        left, right = INFLOW_SPAN
        accumulated, total_qual = weight_span_shared_denominator(matrix, left, right, qualities)
        self.inflow_values = list(accumulated)
        if total_qual > 0.0:
            agg_inflow_obj = rin.Inflow_values()
            for pos, feat_name in enumerate(agg_inflow_obj.__slots__):
                self.inflow_values[pos] = self.inflow_values[pos] / total_qual
                agg_inflow_obj.__setattr__(feat_name, self.inflow_values[pos])
            self.inflow_values = agg_inflow_obj

        ti = add_to_times(times, ti, 'aggregate inflow values')

        left, right = PROFILE_SPAN
        weighted_encoded_profile = weight_profile_columns_span(matrix, left, right, qualities, 0., 0.)

        self.profile = rin.stored_profile_from_encoded(weighted_encoded_profile)

        ti = add_to_times(times, ti, f'weight profile columns {(matrix.shape[0], right - left)=}')

        return times, weighted_encoded_profile

    def weight(self, unweighted_dicts, qualities):
        times = []
        ti = time.time()

        unweighted_centralities = unweighted_dicts[1]
        unweighted_profiles: list[np.ndarray] = unweighted_dicts[0]
        unweighted_inflow_values = unweighted_dicts[2]

        self.centralities = None

        total_qual: float = 0.0

        for mapping_id, centrality_scores in enumerate(unweighted_centralities):
            if centrality_scores is None:
                continue
            qual = qualities[mapping_id]
            if self.centralities is None:
                self.centralities = [0.] * len(centrality_scores.__slots__)

            for pos, feat_name in enumerate(centrality_scores.__slots__):
                cent_score = centrality_scores.__getattribute__(feat_name)
                if cent_score is None or np.isnan(cent_score):
                    continue
                self.centralities[pos] += cent_score * qual

            total_qual += qual
            
        if self.centralities is None:
            # Matches the (times, encoded_profile) contract of the normal exit below;
            # weight_all unpacks this unconditionally. self.profile stays unset, so
            # weight_all still raises the AttributeError this path always produced.
            return times, None

        if total_qual > 0.0:
            agg_cent_obj = rin.Centrality_scores()
            for pos, cent_score in enumerate(self.centralities):
                self.centralities[pos] = self.centralities[pos] / total_qual
            for pos, feat_name in enumerate(agg_cent_obj.__slots__):
                agg_cent_obj.__setattr__(feat_name, self.centralities[pos])
            self.centralities = agg_cent_obj

        ti = add_to_times(times, ti, 'aggregate centralities')

        self.inflow_values = None

        total_qual: float = 0.0

        for mapping_id, inflow_values in enumerate(unweighted_inflow_values):
            if inflow_values is None:
                continue
            qual = qualities[mapping_id]
            if self.inflow_values is None:
                self.inflow_values = [0.] * len(inflow_values.__slots__)

            for pos, feat_name in enumerate(inflow_values.__slots__):
                measure = inflow_values.__getattribute__(feat_name)
                if measure is None or np.isnan(measure):
                    continue
                try:
                    self.inflow_values[pos] += measure * qual
                except RuntimeWarning:
                    self.inflow_values[pos] = float('nan')

            total_qual += qual

        if total_qual > 0.0:
            agg_inflow_obj = rin.Inflow_values()
            for pos, cent_score in enumerate(self.inflow_values):
                self.inflow_values[pos] = self.inflow_values[pos] / total_qual
            for pos, feat_name in enumerate(agg_inflow_obj.__slots__):
                agg_inflow_obj.__setattr__(feat_name, self.inflow_values[pos])
            self.inflow_values = agg_inflow_obj


        ti = add_to_times(times, ti, 'aggregate inflow values')

        weight_profiles:list[np.ndarray] = []
        weight_quals:list[float] = []
        for mapping_id, unweighted_profile in enumerate(unweighted_profiles):
            if unweighted_profile is None:
                continue
            weight_quals.append(qualities[mapping_id])
            weight_profiles.append(unweighted_profile)

            #print(f'{unweighted_profile=}')

        weight_profiles_mat = np.array(weight_profiles, dtype=np.float64)

        # Keeps the IndexError the per-column loop used to raise when no profile
        # survived the filter above, rather than failing later inside the kernel.
        if weight_profiles_mat.ndim != 2:
            raise IndexError(f'no interaction profiles to weight: {weight_profiles_mat.shape=}')

        weighted_encoded_profile = weight_profile_columns(
            weight_profiles_mat, np.array(weight_quals, dtype=np.float64), 0., 0.)
        #print(f'{len(unweighted_profiles)=} {weight_quals=} {weighted_encoded_profile=}')

        self.profile = rin.stored_profile_from_encoded(weighted_encoded_profile)

        ti = add_to_times(times, ti, f'weight profile columns {weight_profiles_mat.shape=}')

        return times, weighted_encoded_profile

        #print(f'{self.profile.to_list()=}')

# --- recommended-residue selection -------------------------------------------
#
# set_recommended_residues walks every mapping of every position, which is the
# same len(seq) x len(annotation_list) count the aggregation loop runs. The
# selection is six running maxima, a four-way bucketing and a per-simple-class
# best, all of which are numeric once two things are turned into integers:
#
#   * the tie-breaks compare mapping_id[0], the structure id string. The caller
#     passes a lexicographic rank per mapping instead, computed once per protein,
#     so a tie is an integer comparison.
#   * the simple class is already an integer code, and the code space is small and
#     fixed, so the interaction recommendations become one slot per code.

_CLASS_CODE_MIN: int = min(rin.SIMPLE_CLASS_STRING_OF_CODE)
_N_CLASS_CODES: int = max(rin.SIMPLE_CLASS_STRING_OF_CODE) - _CLASS_CODE_MIN + 1


def _interaction_class_flags() -> np.ndarray:
    """Which simple class codes are interaction classes.

    Derived from the class strings rather than from the code ranges, so it keeps
    tracking `residue_simple_class.count('interaction') > 0` if the class names
    ever change.
    """
    flags = np.zeros(_N_CLASS_CODES, dtype=np.uint8)
    for code, text in rin.SIMPLE_CLASS_STRING_OF_CODE.items():
        if text is not None and text.count('interaction') > 0:
            flags[code - _CLASS_CODE_MIN] = 1
    return flags


_INTERACTION_CLASS_FLAGS: np.ndarray = _interaction_class_flags()

# Slots of the running maxima the kernel keeps.
_MAX_SEQ_ID = 0
_MAX_QUAL = 1
_MAX_BUCKET_BASE = 2       # + bucket, for the four recommendation buckets


# inline='always', so this never becomes a call into a second dispatcher.
#
# Two of the three call sites below pass a module level int, which numba types as
# Literal[int] and specialises on, while the third passes the runtime int64 of the
# bucket slot - so as a separate dispatcher this needed three signatures where every
# other nested kernel here needs one. Dispatcher.get_call_template only compiles a
# missing signature while _can_compile is set; when it is not, typing falls back to
# whatever overloads already exist, and in the ray workers the int64 one was the one
# that did not get compiled:
#   Invalid use of type(CPUDispatcher(_track_maximum)) with parameters (int64, ...)
#   Known signatures: (Literal[int](0), ...) (Literal[int](1), ...)
# Inlined there is no second dispatcher to type, no per-literal specialisation, and no
# cache entry of its own - the body is compiled as part of _rank_mappings.
@njit(cache=True, inline='always')
def _track_maximum(slot, value, index, structure_ranks, best_value, best_index):
    """One running maximum, ties broken by the structure id's rank.

    The python version indexed the current winner's mapping id unconditionally in
    the tie branch, so an empty tracker (no winner yet, and a value equal to the
    0.0 the maxima start at) raised TypeError or NameError. Here an empty tracker
    simply takes the mapping; the value is still 0.0 afterwards, so the `> 0.`
    gates that pick between the buckets are unaffected.
    """
    if value > best_value[slot]:
        best_value[slot] = value
        best_index[slot] = index
    elif value == best_value[slot]:
        if best_index[slot] == -1 or structure_ranks[index] > structure_ranks[best_index[slot]]:
            best_value[slot] = value
            best_index[slot] = index


@njit(cache=True)
def _rank_mappings(qualities, seq_ids, aa_ids, class_codes, structure_ranks,
                   aggregated_class_code, interaction_flags, class_code_min,
                   buckets, best_value, best_index,
                   interaction_best_qual, interaction_best_index):
    for index in range(qualities.shape[0]):
        qual = qualities[index]

        _track_maximum(_MAX_SEQ_ID, seq_ids[index], index, structure_ranks,
                       best_value, best_index)
        _track_maximum(_MAX_QUAL, qual, index, structure_ranks, best_value, best_index)

        code = class_codes[index]
        identical_aa = aa_ids[index]
        if aggregated_class_code == code:
            bucket = 0 if identical_aa else 1
        else:
            bucket = 2 if identical_aa else 3
        buckets[index] = bucket
        _track_maximum(_MAX_BUCKET_BASE + bucket, qual, index, structure_ranks,
                       best_value, best_index)

        class_slot = code - class_code_min
        if 0 <= class_slot < interaction_flags.shape[0] and interaction_flags[class_slot] == 1:
            # First occurrence wins, later ones only on a strictly better quality.
            if interaction_best_index[class_slot] == -1 or qual > interaction_best_qual[class_slot]:
                interaction_best_index[class_slot] = index
                interaction_best_qual[class_slot] = qual


class EncodedWeightingState:
    """One position's aggregation, accumulated batch by batch instead of in one shot.

    weight_all's encoded path reduces a (mappings x 195) matrix down to ~44 scalars, a
    centrality and inflow vector and a 132-wide profile. Every one of those reductions is
    a running (numerator, denominator) or (count, total), so a protein's annotations can
    be decoded in batches and folded in as they come - which is the point: holding the
    whole matrix cost up to 8.6 GiB on a whole-PDB protein.

    The state is sized by columns, never by mappings: about 5 KB per position, so even a
    1800-residue protein keeps under 10 MB of it. What is *not* here is the handful of
    per-mapping values the recommendation needs; Mappings retains those as thin columns,
    because set_recommended_residues runs after the aggregated class is known and so
    cannot be folded into the accumulation.
    """

    __slots__ = ('n_rows',
                 'dist_nom', 'dist_denom', 'float_nom', 'float_denom',
                 'surface_nom', 'surface_denom', 'bool_true', 'prob_present',
                 'cat_seen', 'cat_weights', 'cat_n_seen', 'cat_overflow',
                 'cent_acc', 'cent_qual', 'inflow_acc', 'inflow_qual',
                 'profile_nom', 'profile_denom', 'profile_any')

    def __init__(self):
        layout = _STRUCTURAL_COLUMN_LAYOUT
        self.n_rows = 0

        def pair(n):
            return np.zeros(n, dtype=np.float64), np.zeros(n, dtype=np.float64)

        self.dist_nom, self.dist_denom = pair(layout.distance_columns.shape[0])
        self.float_nom, self.float_denom = pair(layout.float_columns.shape[0])
        self.surface_nom, self.surface_denom = pair(layout.surface_columns.shape[0])
        self.bool_true = np.zeros(layout.boolean_columns.shape[0], dtype=np.float64)
        self.prob_present = np.zeros(layout.probability_columns.shape[0], dtype=np.float64)

        n_cat = layout.categorical_columns.shape[0]
        self.cat_seen = np.zeros((n_cat, _MAX_MAJORITY_VALUES), dtype=np.float64)
        self.cat_weights = np.zeros((n_cat, _MAX_MAJORITY_VALUES), dtype=np.float64)
        self.cat_n_seen = np.zeros(n_cat, dtype=np.int64)
        self.cat_overflow = 0

        cent_left, cent_right = CENTRALITIES_SPAN
        self.cent_acc = np.zeros(cent_right - cent_left, dtype=np.float64)
        self.cent_qual = np.zeros(1, dtype=np.float64)

        inflow_left, inflow_right = INFLOW_SPAN
        self.inflow_acc = np.zeros(inflow_right - inflow_left, dtype=np.float64)
        self.inflow_qual = np.zeros(1, dtype=np.float64)

        profile_left, profile_right = PROFILE_SPAN
        span = profile_right - profile_left
        self.profile_nom = np.zeros(span, dtype=np.float64)
        self.profile_denom = np.zeros(span, dtype=np.float64)
        self.profile_any = np.zeros(span, dtype=np.uint8)

    def accumulate(self, matrix: np.ndarray, qualities: np.ndarray, covs: np.ndarray,
                   config) -> None:
        """Fold one batch's rows in. Batches must arrive in annotation order, because
        weight_columns_majority resolves ties to the first-seen candidate."""
        if matrix.shape[0] == 0:
            return

        layout = _STRUCTURAL_COLUMN_LAYOUT
        self.n_rows += matrix.shape[0]

        accumulate_columns_dist(matrix, layout.distance_columns, qualities,
                                config.long_distance_threshold, self.dist_nom, self.dist_denom)
        accumulate_columns_float(matrix, layout.float_columns, qualities,
                                 self.float_nom, self.float_denom)
        accumulate_columns_cov(matrix, layout.surface_columns, qualities, covs, 2,
                               self.surface_nom, self.surface_denom)
        accumulate_columns_bool_propensity(matrix, layout.boolean_columns, self.bool_true)
        self.cat_overflow += accumulate_columns_majority(
            matrix, layout.categorical_columns, qualities,
            self.cat_seen, self.cat_weights, self.cat_n_seen)
        accumulate_columns_propensity(matrix, layout.probability_columns, self.prob_present)

        cent_left, cent_right = CENTRALITIES_SPAN
        accumulate_span_shared_denominator(matrix, cent_left, cent_right, qualities,
                                           self.cent_acc, self.cent_qual)
        inflow_left, inflow_right = INFLOW_SPAN
        accumulate_span_shared_denominator(matrix, inflow_left, inflow_right, qualities,
                                           self.inflow_acc, self.inflow_qual)
        profile_left, profile_right = PROFILE_SPAN
        accumulate_profile_columns_span(matrix, profile_left, profile_right, qualities, 0.,
                                        self.profile_nom, self.profile_denom, self.profile_any)

    def apply_structural(self, structural_features) -> None:
        """Set what set_values_from_encoded sets, from the accumulated state."""
        if self.cat_overflow:
            raise ValueError(
                f'{self.cat_overflow} categorical values did not fit the majority state '
                f'(cap {_MAX_MAJORITY_VALUES}); the winner could differ from the '
                f'unbatched result, so this is raised rather than guessed at')

        layout = _STRUCTURAL_COLUMN_LAYOUT

        weighted = finalize_columns_nom_denom(self.dist_nom, self.dist_denom)
        for i, feature_name in enumerate(layout.distance_names):
            structural_features.__setattr__(feature_name, weighted[i])

        weighted = finalize_columns_nom_denom(self.float_nom, self.float_denom)
        for i, feature_name in enumerate(layout.float_names):
            structural_features.__setattr__(feature_name, weighted[i])

        weighted = finalize_columns_nom_denom(self.surface_nom, self.surface_denom)
        for i, feature_name in enumerate(layout.surface_names):
            structural_features.__setattr__(feature_name, weighted[i])

        weighted = finalize_columns_majority(self.cat_seen, self.cat_weights, self.cat_n_seen)
        for i, feature_name in enumerate(layout.categorical_names):
            value = weighted[i]
            # weight_majority returns None when nothing outvoted a zero quality.
            structural_features.__setattr__(feature_name, None if value != value else value)

        if self.n_rows > 0:
            weighted = finalize_columns_over_rows(self.prob_present, self.n_rows)
            for i, feature_name in enumerate(layout.probability_names):
                structural_features.__setattr__(feature_name, weighted[i])

            weighted = finalize_columns_over_rows(self.bool_true, self.n_rows)
            for i, feature_name in enumerate(layout.boolean_names):
                structural_features.__setattr__(feature_name, weighted[i])
        else:
            for feature_name in layout.probability_names:
                structural_features.__setattr__(feature_name, None)

            for i, feature_name in enumerate(layout.boolean_names):
                structural_features.__setattr__(feature_name, None)

        # weight_function returns None for anything in no category.
        for feature_name in layout.uncategorised_names:
            structural_features.__setattr__(feature_name, None)

    def apply_rin(self, rin_based_features) -> np.ndarray | None:
        """Set what weight_from_encoded sets; returns the encoded profile, or None when
        no mapping contributed centralities - the same early exit weight() has."""
        rin_based_features.centralities = None
        if self.n_rows > 0:
            total_qual = self.cent_qual[0]
            rin_based_features.centralities = list(self.cent_acc)
            if total_qual > 0.0:
                agg_cent_obj = rin.Centrality_scores()
                for pos, feat_name in enumerate(agg_cent_obj.__slots__):
                    rin_based_features.centralities[pos] = rin_based_features.centralities[pos] / total_qual
                    agg_cent_obj.__setattr__(feat_name, rin_based_features.centralities[pos])
                rin_based_features.centralities = agg_cent_obj

        if rin_based_features.centralities is None:
            return None

        rin_based_features.inflow_values = None
        total_qual = self.inflow_qual[0]
        rin_based_features.inflow_values = list(self.inflow_acc)
        if total_qual > 0.0:
            agg_inflow_obj = rin.Inflow_values()
            for pos, feat_name in enumerate(agg_inflow_obj.__slots__):
                rin_based_features.inflow_values[pos] = rin_based_features.inflow_values[pos] / total_qual
                agg_inflow_obj.__setattr__(feat_name, rin_based_features.inflow_values[pos])
            rin_based_features.inflow_values = agg_inflow_obj

        weighted_encoded_profile = finalize_profile_columns_span(
            self.profile_nom, self.profile_denom, self.profile_any, 0.)
        rin_based_features.profile = rin.stored_profile_from_encoded(weighted_encoded_profile)
        return weighted_encoded_profile


class MappingIds:
    """(pdb_id, chain, res_nr) per mapping, built on demand from two int32 columns.

    The aggregation loop knows a mapping as an annotation index and a residue contig
    index; the triple is what the recommendation output and calc_mm_features read.
    Materialising it per mapping cost about 92 B times positions times annotations -
    0.95 GiB for a protein with 37k annotations, held for that protein's whole
    position loop - while only the first few entries of each position's
    quality-ordered list are ever looked at.

    structure_ids and chains are the protein's own per-annotation columns
    (aggregation._gather_structure_columns), held by reference, so every position of
    the protein shares one copy of the strings.

    Indexable and sized like the list it replaces, so every reader - including
    calc_mm_features walking recommendation_order - is unchanged.
    """

    __slots__ = ('anno_index', 'res_contig', 'structure_ids', 'chains')

    def __init__(self, anno_index: np.ndarray, res_contig: np.ndarray,
                 structure_ids: list[str], chains: list[str]):
        self.anno_index = anno_index
        self.res_contig = res_contig
        self.structure_ids = structure_ids
        self.chains = chains

    def __len__(self) -> int:
        return self.anno_index.shape[0]

    def __getitem__(self, index: int) -> tuple[str, str, int]:
        anno = self.anno_index[index]
        return (self.structure_ids[anno], self.chains[anno], int(self.res_contig[index]))

    def __iter__(self):
        for index in range(self.anno_index.shape[0]):
            yield self[index]

    def reorder(self, order: np.ndarray) -> 'MappingIds':
        """The same mappings in a new order, still without materialising a triple."""
        return MappingIds(self.anno_index[order], self.res_contig[order],
                          self.structure_ids, self.chains)

    def __repr__(self) -> str:
        # Only ever reached from the verbosity >= 6 log lines, which want the triples.
        return repr(list(self))


class Mappings(Slotted_obj):
    __slots__ = [
                    'qualities',                            'covs',                 'seq_ids',
                    'aa_ids',                               'max_seq_res',          'recommended_res',
                    'interaction_recommendations',          'recommendation_order', 'resolutions',
                    'res_aas',                              'amount_of_structures', 'structural_features',
                    'structural_unweighted_feature_dicts',  'rin_based_features',   'rin_based_unweighted_feature_dicts',
                    'integrated_features',                  'mapping_ids',          'encoded_rows',
                    'simple_class_codes',                   'structure_ranks',
                    'weighting_state',                      'thin_columns'
                ]

    # New slots go at the end with a False mask, so the serialised layout - and
    # therefore everything already packed into the database - is unchanged.
    slot_mask = [
        False, False, False,
        False, True, True,
        False, True, False,
        False, True, True,
        False, True, False,
        True,  False, False,
        False, False,
        False, False
    ]

    def __init__(self, init_none = False):
        self.mapping_ids = []
        self.encoded_rows = []
        # Set only on the batched path (see accumulate_encoded_mappings); None means the
        # position was aggregated in one shot and weight_all applies as before.
        self.weighting_state = None
        self.thin_columns = None
        self.simple_class_codes = []
        self.structure_ranks = None
        self.qualities = []
        self.seq_ids = []
        self.covs = []

        self.amount_of_structures = 0
        self.aa_ids = []
        self.res_aas = []
        self.resolutions = []

        self.recommendation_order: list[tuple[str, str, int | str]] = []

        self.structural_features = Structural_features(init_none=init_none)
        self.structural_unweighted_feature_dicts: dict[str, list] = self.structural_features.create_unweighted_feature_dicts()


        self.rin_based_features = RIN_based_features(init_none=init_none)
        self.rin_based_unweighted_feature_dicts = self.rin_based_features.create_unweighted_feature_lists()

        self.integrated_features = Integrated_features(init_none=init_none)

    def deconstruct(self):
        doomsday_protocol(self)

    def add_encoded_mapping(self, mapping_id, encoded_row, quality, seq_id, cov,
                            identical_aa, resolution, res_aa, simple_class_code):
        # The whole encoded residue row is kept as-is. Every numeric structural
        # feature, the interaction profile, the centralities and the inflow
        # values are columns of it, so nothing has to be unpacked into a
        # 44-key dict and spread over 44 python lists per mapping - the columns
        # are read once, in the kernels, at weighting time.
        self.amount_of_structures += 1
        self.mapping_ids.append(mapping_id)
        self.encoded_rows.append(encoded_row)
        self.qualities.append(quality)
        self.seq_ids.append(seq_id)
        self.covs.append(cov)
        self.aa_ids.append(identical_aa)
        self.res_aas.append(res_aa)
        self.resolutions.append(resolution)
        self.simple_class_codes.append(simple_class_code)

    def set_encoded_mappings(self, mapping_ids, encoded_rows, qualities, seq_ids, covs,
                             aa_ids, resolutions, res_aas, simple_class_codes,
                             structure_ranks=None):
        """Install a whole position's mappings as columns rather than row by row.

        The aggregation loop collects into numba-filled arrays now, so the nine
        per-mapping appends of add_encoded_mapping become nine assignments, and
        encoded_rows arrives as the (mappings x columns) matrix that weight_all
        used to rebuild from a list of row views on every position. None of these
        slots is serialised (see slot_mask), so the packed representation is
        untouched by them being arrays.

        structure_ranks is the lexicographic rank of each mapping's structure id
        among this protein's structures, which is what set_recommended_residues
        needs to break quality ties without comparing strings per mapping.
        """
        self.structure_ranks = structure_ranks
        self.mapping_ids = mapping_ids
        self.encoded_rows = encoded_rows
        self.qualities = qualities
        self.seq_ids = seq_ids
        self.covs = covs
        self.aa_ids = aa_ids
        self.resolutions = resolutions
        self.res_aas = res_aas
        self.simple_class_codes = simple_class_codes
        self.amount_of_structures = len(mapping_ids)

    def accumulate_encoded_mappings(self, matrix, qualities, seq_ids, covs, aa_ids,
                                    resolutions, res_aas, simple_class_codes,
                                    structure_ranks, anno_index, res_contig,
                                    surface_values, mc_surface_values, sc_surface_values,
                                    config):
        """One batch of this position's mappings, folded in instead of installed.

        The wide (mappings x 195) matrix is reduced into the accumulator and dropped with
        the batch; only the thin per-mapping values the recommendation needs are kept, at
        61 bytes a mapping against 780. Batches must arrive in annotation order - see
        EncodedWeightingState.accumulate.

        The three surface columns are kept as well: integrated_features.generate counts
        locations over them, and takes np.std of the side-chain one, which is a two-pass
        algorithm and so cannot be streamed without changing the result.
        """
        if self.weighting_state is None:
            self.weighting_state = EncodedWeightingState()
            self.thin_columns = {name: [] for name in (
                'qualities', 'seq_ids', 'covs', 'aa_ids', 'resolutions', 'res_aas',
                'simple_class_codes', 'structure_ranks', 'anno_index', 'res_contig',
                'surface_values', 'mc_surface_values', 'sc_surface_values')}

        self.weighting_state.accumulate(matrix, qualities, covs, config)

        columns = self.thin_columns
        columns['qualities'].append(qualities)
        columns['seq_ids'].append(seq_ids)
        columns['covs'].append(covs)
        columns['aa_ids'].append(aa_ids)
        columns['resolutions'].append(resolutions)
        columns['res_aas'].append(res_aas)
        columns['simple_class_codes'].append(simple_class_codes)
        columns['structure_ranks'].append(structure_ranks)
        columns['anno_index'].append(anno_index)
        columns['res_contig'].append(res_contig)
        columns['surface_values'].append(surface_values)
        columns['mc_surface_values'].append(mc_surface_values)
        columns['sc_surface_values'].append(sc_surface_values)

    def finalize_encoded(self, config, structure_ids, chains, for_indel_aggregation=False):
        """What weight_all does, from the accumulated state rather than from a matrix.

        Deliberately the same sequence of calls in the same order - the structural
        features, then the rin features, then integrated_features.generate, then the
        recommendation - so that the batched result is the unbatched one.
        """
        times = []
        ti = time.time()

        if self.weighting_state is None:
            # No batch ever contributed, so there are no mappings at this position. Same
            # shape of answer weight_all gives for an empty matrix.
            return self.weight_all(config, for_indel_aggregation=for_indel_aggregation)

        def joined(name, dtype):
            # The pieces are released as they are merged. Keeping them until the end
            # would put a second copy of every retained column on the peak, which is the
            # one thing batching is supposed to avoid; each name is joined exactly once,
            # so popping it is safe.
            parts = self.thin_columns.pop(name)
            if len(parts) == 1:
                return parts[0]
            if not parts:
                return np.empty(0, dtype=dtype)
            merged = np.concatenate(parts)
            del parts[:]
            return merged

        self.qualities = joined('qualities', np.float64)
        self.seq_ids = joined('seq_ids', np.float64)
        self.covs = joined('covs', np.float64)
        self.aa_ids = joined('aa_ids', np.bool_)
        self.resolutions = joined('resolutions', np.float64)
        self.res_aas = joined('res_aas', np.float32)
        self.simple_class_codes = joined('simple_class_codes', np.int64)
        self.structure_ranks = joined('structure_ranks', np.int64)
        self.mapping_ids = MappingIds(joined('anno_index', np.int32),
                                      joined('res_contig', np.int32),
                                      structure_ids, chains)
        self.amount_of_structures = len(self.mapping_ids)

        state = self.weighting_state
        state.apply_structural(self.structural_features)
        ti = add_to_times(times, ti, 'weight structural features')

        encoded_profile = state.apply_rin(self.rin_based_features)
        try:
            profile = self.rin_based_features.profile
            if profile is None:
                raise AttributeError('profile is None in finalize_encoded')
        except AttributeError as err:
            raise AttributeError(f'{err=} profile not set in finalize_encoded')

        ti = add_to_times(times, ti, 'weight rin_based features')

        self.integrated_features.generate(
            self.structural_features, config, profile,
            joined('surface_values', np.float32),
            joined('mc_surface_values', np.float32),
            joined('sc_surface_values', np.float32),
            encoded_profile=encoded_profile)

        ti = add_to_times(times, ti, 'generate integrated features')

        if not for_indel_aggregation:
            self.set_recommended_residues()

        ti = add_to_times(times, ti, 'set_recommended_residues')

        # The state and the batch columns are dead once the features are set, and this
        # object is about to be packed.
        self.weighting_state = None
        self.thin_columns = None

        return None, times

    def add_mapping(self, mapping_id, mapping):
        quality, seq_id, cov, structural_feature_dict, rin_based_features, identical_aa, resolution, res_aa = mapping
        self.amount_of_structures += 1
        self.mapping_ids.append(mapping_id)
        self.qualities.append(quality)
        self.seq_ids.append(seq_id)
        self.covs.append(cov)
        for feature_name in structural_feature_dict:
            self.structural_unweighted_feature_dicts[feature_name].append(structural_feature_dict[feature_name])

        self.rin_based_unweighted_feature_dicts[0].append(rin_based_features[0])
        self.rin_based_unweighted_feature_dicts[1].append(rin_based_features[1])
        self.rin_based_unweighted_feature_dicts[2].append(rin_based_features[2])

        self.aa_ids.append(identical_aa)
        self.res_aas.append(res_aa)
        self.resolutions.append(resolution)

    def add_result(self, mapping_id: int | str, structural_features: Structural_features, rin_based_features: RIN_based_features, quality: float, seq_id: float, cov: float) -> None:

        if mapping_id not in self.qualities:
            self.amount_of_structures += 1
        self.qualities[mapping_id] = quality
        self.seq_ids[mapping_id] = seq_id
        self.covs[mapping_id] = cov

        for feat_name in structural_features.get_feature_names():
            try:
                self.structural_unweighted_feature_dicts[feat_name][mapping_id] = structural_features.__getattribute__(feat_name)
            except AttributeError:
                self.structural_unweighted_feature_dicts[feat_name][mapping_id] = None

        try:
            self.rin_based_unweighted_feature_dicts[0][mapping_id] = rin.as_stored_profile(rin_based_features.profile)
        except AttributeError:
            self.rin_based_unweighted_feature_dicts[0][mapping_id] = None
        try:
            self.rin_based_unweighted_feature_dicts[1][mapping_id] = rin_based_features.centralities
        except AttributeError:
            self.rin_based_unweighted_feature_dicts[1][mapping_id] = None
        try:
            self.rin_based_unweighted_feature_dicts[2][mapping_id] = rin_based_features.inflow_values
        except AttributeError:
            self.rin_based_unweighted_feature_dicts[2][mapping_id] = None

        self.aa_ids[mapping_id] = True
        self.res_aas[mapping_id] = None


    def weight_all(self, config, for_indel_aggregation = False) -> None | str:

        times = []
        ti = time.time()

        error = None

        # Two storage modes: the aggregation loop feeds whole encoded residue
        # rows (matrix mode), while the indel path still builds Structural_features
        # objects and goes through add_result (object mode).
        encoded_mode = len(self.encoded_rows) > 0

        if encoded_mode:
            # asarray, not array: set_encoded_mappings hands over a matrix and two
            # float64 columns that are already in the right shape and dtype, and
            # none of the weighting kernels writes to them. The list form (from
            # add_encoded_mapping) still converts here as before.
            matrix = np.asarray(self.encoded_rows, dtype=np.float32)
            qualities = np.asarray(self.qualities, dtype=np.float64)
            covs = np.asarray(self.covs, dtype=np.float64)

            self.structural_features.set_values_from_encoded(matrix, qualities, covs, config)

            ti = add_to_times(times, ti, 'weight structural features')

            rin_weighting_times, encoded_profile = self.rin_based_features.weight_from_encoded(matrix, qualities)
        else:
            for feature_name in self.structural_unweighted_feature_dicts:
                weighted_value = self.structural_features.weight_function(feature_name, self.structural_unweighted_feature_dicts[feature_name], self.qualities, self.covs, config)
                self.structural_features.set_value_by_name(feature_name, weighted_value)

            ti = add_to_times(times, ti, 'weight structural features')

            rin_weighting_times, encoded_profile = self.rin_based_features.weight(self.rin_based_unweighted_feature_dicts, self.qualities)

        try:
            profile = self.rin_based_features.profile
            if profile is None:
                raise AttributeError(f'profile is None in weight_all {len(self.rin_based_unweighted_feature_dicts)=}')
        except AttributeError as err:
            raise AttributeError(f'{err=} profile not set in weight_all')

        times.append(rin_weighting_times)

        ti = add_to_times(times, ti, 'weight rin_based features')
            
        if encoded_mode:
            surface_values = matrix[:, ENCODED_COLUMN_OF_SLOT['RSA']]
            mc_surface_values = matrix[:, ENCODED_COLUMN_OF_SLOT['relative_main_chain_acc']]
            sc_surface_values = matrix[:, ENCODED_COLUMN_OF_SLOT['relative_side_chain_acc']]
        else:
            surface_values = self.structural_unweighted_feature_dicts['surface_value']
            mc_surface_values = self.structural_unweighted_feature_dicts['mainchain_surface_value']
            sc_surface_values = self.structural_unweighted_feature_dicts['sidechain_surface_value']

        self.integrated_features.generate(self.structural_features, config, profile,
                                          surface_values, mc_surface_values, sc_surface_values,
                                          encoded_profile=encoded_profile)

        ti = add_to_times(times, ti, 'generate integrated features')

        if not for_indel_aggregation:
            self.set_recommended_residues()  

        ti = add_to_times(times, ti, 'set_recommended_residues')

        return error, times
        
    def _residue_tuple(self, index):
        return (self.mapping_ids[index], self.seq_ids[index], self.covs[index],
                self.resolutions[index], self.res_aas[index])

    def _set_recommended_residues_ranked(self, aggregated_class_code):
        """The column form of set_recommended_residues; see _rank_mappings."""
        qualities = self.qualities
        n_mappings = len(qualities)

        buckets = np.empty(n_mappings, dtype=np.int64)
        best_value = np.zeros(_MAX_BUCKET_BASE + 4, dtype=np.float64)
        best_index = np.full(_MAX_BUCKET_BASE + 4, -1, dtype=np.int64)
        interaction_best_qual = np.zeros(_N_CLASS_CODES, dtype=np.float64)
        interaction_best_index = np.full(_N_CLASS_CODES, -1, dtype=np.int64)

        _rank_mappings(
            qualities, self.seq_ids, self.aa_ids,
            np.asarray(self.simple_class_codes, dtype=np.int64), self.structure_ranks,
            np.int64(aggregated_class_code), _INTERACTION_CLASS_FLAGS,
            np.int64(_CLASS_CODE_MIN),
            buckets, best_value, best_index,
            interaction_best_qual, interaction_best_index)

        # Buckets in order, each by descending quality, ties in mapping order -
        # which is what sorted(..., reverse=True) gives, since it stays stable.
        # The trailing index key reproduces that stability.
        order = np.lexsort((np.arange(n_mappings), -qualities, buckets))
        if isinstance(self.mapping_ids, MappingIds):
            # Permuted columns, not a list of triples: calc_mm_features reads only
            # the first few of these and the list cost about a gigabyte per heavy
            # protein (see MappingIds).
            self.recommendation_order = self.mapping_ids.reorder(order)
        else:
            self.recommendation_order = [self.mapping_ids[index] for index in order.tolist()]

        self.interaction_recommendations = {}
        for class_slot in range(_N_CLASS_CODES):
            index = interaction_best_index[class_slot]
            if index == -1:
                continue
            self.interaction_recommendations[
                rin.SIMPLE_CLASS_STRING_OF_CODE[class_slot + _CLASS_CODE_MIN]] = (
                    qualities[index], self.mapping_ids[index],
                    self.seq_ids[index], self.covs[index])

        self.max_seq_res = None
        if best_index[_MAX_SEQ_ID] != -1:
            winner = best_index[_MAX_SEQ_ID]
            self.max_seq_res = qualities[winner], self._residue_tuple(winner)

        for bucket in range(3):
            slot = _MAX_BUCKET_BASE + bucket
            if best_value[slot] > 0.:
                self.recommended_res = best_value[slot], self._residue_tuple(best_index[slot])
                break
        else:
            # No mapping with an identical simple class or identical aa, which
            # should not happen - fall back to the best quality overall.
            winner = best_index[_MAX_QUAL]
            self.recommended_res = (
                best_value[_MAX_QUAL],
                self._residue_tuple(winner) if winner != -1 else None)

    def set_recommended_residues(self):
        # A streamed position has no encoded_rows - finalize_encoded reduced them
        # batch by batch - but it is the same encoded data and wants the same path.
        encoded_mode = len(self.encoded_rows) > 0 or self.structure_ranks is not None
        if encoded_mode:
            # -999 is deliberately outside the code space, so an unexpected
            # aggregated class matches nothing - as an unknown string did.
            aggregated_simple_class_code = rin.SIMPLE_CLASS_CODE_OF_STRING.get(
                self.integrated_features.simple_class, -999)

            # The column path needs the structure ranks the aggregation loop
            # supplies, and it orders by lexsort, which puts a NaN quality last
            # rather than wherever timsort's comparisons happened to leave it. A
            # NaN quality needs a NaN coverage or sequence identity in the
            # database, so this is insurance rather than a live case.
            if self.structure_ranks is not None and not np.isnan(self.qualities).any():
                self._set_recommended_residues_ranked(aggregated_simple_class_code)
                self._format_recommendations()
                return

        max_seq_id = 0.
        max_identical_aa_identical_class = 0.
        max_identical_aa = 0.
        max_identical_class = 0.
        max_should_not_happen = 0.
        max_qual = 0.
        max_qual_res = None
        # Initialised so the tie branches below can tell "no winner yet" from a
        # real one. They used to index the winner unconditionally, which raised
        # TypeError (max_seq_res, max_qual_res) or NameError (the four bucket
        # winners) as soon as the first mapping's value equalled the 0.0 the
        # maxima start at - a zero sequence identity or quality was enough.
        max_identical_aa_identical_class_res = None
        max_identical_class_res = None
        max_identical_aa_res = None
        max_should_not_happen_res = None
        self.max_seq_res = None
        self.recommended_res = None
        self.interaction_recommendations = {}

        identical_aa_identical_class_recs = []
        identical_class_recs = []
        identical_aa_recs = []
        other_recs = []

        for mapping_nr, seq_id in enumerate(self.seq_ids):
            mapping_id = self.mapping_ids[mapping_nr]
            seq_id = self.seq_ids[mapping_nr]
            cov = self.covs[mapping_nr]
            qual = self.qualities[mapping_nr]
            if qual is None:
                qual = 0.
            resolution = self.resolutions[mapping_nr]
            res_aa = self.res_aas[mapping_nr]

            if seq_id > max_seq_id:
                max_seq_id = seq_id
                self.max_seq_res = qual, (mapping_id, seq_id, cov, resolution, res_aa)
            elif seq_id == max_seq_id:
                if self.max_seq_res is None or mapping_id[0] > self.max_seq_res[1][0][0]:
                    max_seq_id = seq_id
                    self.max_seq_res = qual, (mapping_id, seq_id, cov, resolution, res_aa)

            if qual > max_qual:
                max_qual = qual
                max_qual_res = mapping_id, seq_id, cov, resolution, res_aa
            elif qual == max_qual:
                if max_qual_res is None or mapping_id[0] > max_qual_res[0][0]:
                    max_qual = qual
                    max_qual_res = mapping_id, seq_id, cov, resolution, res_aa

            # In matrix mode the per-mapping simple class is an integer code, so
            # the position's aggregated simple class is translated into the same
            # code space once, above the loop, instead of comparing strings that
            # had to be built for every single mapping.
            if encoded_mode:
                code = self.simple_class_codes[mapping_nr]
                matches_class = aggregated_simple_class_code == code
                # Only the interaction_recommendations block below still wants
                # the string, and it is one dict lookup per mapping rather than
                # a formatted string built during the aggregation loop.
                residue_simple_class = rin.SIMPLE_CLASS_STRING_OF_CODE.get(code)
            else:
                residue_simple_class = self.structural_unweighted_feature_dicts['rin_simple_class'][mapping_nr]
                matches_class = self.integrated_features.simple_class == residue_simple_class

            if matches_class:
                if self.aa_ids[mapping_nr]:
                    identical_aa_identical_class_recs.append((qual, mapping_id))
                    if qual > max_identical_aa_identical_class:
                        max_identical_aa_identical_class = qual
                        max_identical_aa_identical_class_res = mapping_id, seq_id, cov, resolution, res_aa
                    elif qual == max_identical_aa_identical_class:
                        if max_identical_aa_identical_class_res is None or mapping_id[0] > max_identical_aa_identical_class_res[0][0]:
                            max_identical_aa_identical_class = qual
                            max_identical_aa_identical_class_res = mapping_id, seq_id, cov, resolution, res_aa
                else:
                    identical_class_recs.append((qual, mapping_id))
                    if qual > max_identical_class:
                        max_identical_class = qual
                        max_identical_class_res = mapping_id, seq_id, cov, resolution, res_aa
                    elif qual == max_identical_class:
                        if max_identical_class_res is None or mapping_id[0] > max_identical_class_res[0][0]:
                            max_identical_class = qual
                            max_identical_class_res = mapping_id, seq_id, cov, resolution, res_aa
            else:
                if self.aa_ids[mapping_nr]:
                    identical_aa_recs.append((qual, mapping_id))
                    if qual > max_identical_aa:
                        max_identical_aa = qual
                        max_identical_aa_res = mapping_id, seq_id, cov, resolution, res_aa
                    elif qual == max_identical_aa:
                        if max_identical_aa_res is None or mapping_id[0] > max_identical_aa_res[0][0]:
                            max_identical_aa = qual
                            max_identical_aa_res = mapping_id, seq_id, cov, resolution, res_aa
                else:
                    other_recs.append((qual, mapping_id))
                    if qual > max_should_not_happen:
                        max_should_not_happen = qual
                        max_should_not_happen_res = mapping_id, seq_id, cov, resolution, res_aa
                    elif qual == max_should_not_happen:
                        if max_should_not_happen_res is None or mapping_id[0] > max_should_not_happen_res[0][0]:
                            max_should_not_happen = qual
                            max_should_not_happen_res = mapping_id, seq_id, cov, resolution, res_aa

            if residue_simple_class is not None:
                # Every simple class this can hold is either a location or
                # '<interaction type> interaction', all lower case, so the old
                # capitalised test never matched and this map stayed empty.
                if residue_simple_class.count('interaction') > 0:
                    if residue_simple_class not in self.interaction_recommendations:
                        self.interaction_recommendations[residue_simple_class] = (qual, mapping_id, seq_id, cov)
                    elif qual > self.interaction_recommendations[residue_simple_class][0]:
                        self.interaction_recommendations[residue_simple_class] = (qual, mapping_id, seq_id, cov)

        recommendation_order = [x[1] for x in sorted(identical_aa_identical_class_recs, key = lambda x:x[0], reverse = True)]
        recommendation_order += [x[1] for x in sorted(identical_class_recs, key = lambda x:x[0], reverse = True)]
        recommendation_order += [x[1] for x in sorted(identical_aa_recs, key = lambda x:x[0], reverse = True)]
        recommendation_order += [x[1] for x in sorted(other_recs, key = lambda x:x[0], reverse = True)]
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

        self._format_recommendations()

    def _format_recommendations(self):
        """Turn the two selected (qual, residue) pairs into their stored strings."""
        if self.recommended_res is not None:
            if self.recommended_res[1] is not None:
                qual, (recommended_res, seq_id, cov, resolution, res_aa) = self.recommended_res
                pdb_id, chain, res_nr = recommended_res
                self.recommended_res = '%s:%s %s:%s;%1.2f;%1.2f;%1.1f' % (pdb_id, chain, res_nr, res_aa, seq_id, cov, resolution)
            else:
                self.recommended_res = None

        if self.max_seq_res is not None:
            if self.max_seq_res[1] is not None:
                qual, (recommended_res, seq_id, cov, resolution, res_aa) = self.max_seq_res
                pdb_id, chain, res_nr = recommended_res
                self.max_seq_res = '%s:%s %s:%s;%1.2f;%1.2f;%1.1f' % (pdb_id, chain, res_nr, res_aa, seq_id, cov, resolution)
            else:
                self.max_seq_res = None

    def get_recommended_res_str(self):
        return self.recommended_res

    def get_max_seq_structure_res_str(self):
        return self.max_seq_res
