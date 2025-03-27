import sys
import re
import traceback

from bisect import bisect_left
from more_itertools import consecutive_groups

from structman.lib.sdsc.consts import codons, residues
from structman.lib.sdsc import indel as indel_package
from structman.lib.sdsc import mutations as mutation_package
from structman.lib.sdsc import position as position_package
from structman.lib.globalAlignment import call_biopython_alignment


def extract_dedicated_tags(tags):

    # Parses tags and extracts tags with dedicated tag identifiers
    # Currently only 'gene' is a valid dedicated tag identifier
    #
    # Input:
    #   tags - either None or a set of tags
    #
    # Output:
    #   extracted_tags - tags witgout the dedicated tags
    #   gene_id        - gene identifier

    gene_id = None
    if tags is None:
        tags = set()

    extracted_tags = set()
    for tag in tags:
        tag_tuple = tag.split(':')
        if len(tag_tuple) != 2:
            extracted_tags.add(tag)
            continue
        tag_id, tag_content = tag_tuple
        if tag_id == 'gene':
            gene_id = tag_content
        else:
            extracted_tags.add(tag)

    return extracted_tags, gene_id


def three_to_one(three_letter):
    three_letter = three_letter.upper()
    one_letter = residues.THREE_TO_ONE[three_letter][0]
    return one_letter


def processHgvsPro(aachange):
    aachange = aachange[2:] #get rid of 'p.'
    aa1_three_letter = aachange[:3]
    aa2_three_letter = aachange[-3:]
    pos = int(aachange[3:-3])
    aa1 = three_to_one(aa1_three_letter)
    aa2 = three_to_one(aa2_three_letter)
    aachange = f'{aa1}{pos}{aa2}'
    return aachange, aa1, aa2, pos



def processAAChange(aachange, pdb_style=False):
    if aachange[:2] == 'p.':
        return processHgvsPro(aachange)

    if ord(aachange[0]) > 47 and ord(aachange[0]) < 58:  # if first char is a number
        aa1 = 'X'
        aachange = "X%s" % aachange
    else:
        aa1 = aachange[0]

    if ord(aachange[-1]) > 47 and ord(aachange[-1]) < 58:  # if last char is a number
        pos = int(aachange[1:])
        aa2 = None
    elif pdb_style and aachange.count('_') == 1:
        aachange = aachange.replace('_', '')
        aa2 = aachange[-1]
        pdb_res_nr = aachange[1:-1]
        return aachange, aa1, aa2, pdb_res_nr
    elif pdb_style:
        aa2 = None
        pdb_res_nr = aachange[1:]
        return aachange, aa1, aa2, pdb_res_nr
    else:
        aa2 = aachange[-1]
        pos = int(aachange[1:-1])

    return aachange, aa1, aa2, pos



def process_variants(config, variants):
    
    multi_mutations = []
    positions = []
    
    skipped_aachanges = 0
    pos_set = set()

    for variant in variants:
        if variant == '':
            skipped_aachanges += 1
            continue
        if variant.count('delins') == 1:
            indel = indel_package.Substitution(raw_str=variant)
            multi_mutations.append(indel)

        elif variant.count('del') == 1:
            indel = indel_package.Deletion(raw_str=variant)
            multi_mutations.append(indel)

        elif variant.count('ins') == 1:
            indel = indel_package.Insertion(raw_str=variant)
            multi_mutations.append(indel)

        else:
            indel = None
            try:
                aachange, aa1, aa2, pos = processAAChange(variant) 
            except:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                config.errorlog.add_warning("Position/Mutation String Format Error: %s\n%s\n%s\n%s" % (variant, str(e), str(f), str(g)))
                return None

            if aa2 is None:
                if aachange not in pos_set:
                    position = position_package.Position(pos=pos, wt_aa=aa1)
                    positions.append(position)
                    pos_set.add(aachange)
            else:
                if aachange not in pos_set:
                    position = position_package.Position(pos=pos, wt_aa=aa1, mut_aas=set(aa2))
                    positions.append(position)
                    pos_set.add(aachange)
                else:
                    position = pos
                multi_mutations.append((position, aa2))

    return positions, multi_mutations


def process_mutations_str(config, mutation_str, tags, pos_set = None, pdb_style=False):
    if isinstance(tags, str):
        tags = set(tags.split(','))
    tags, gene_id = extract_dedicated_tags(tags)
    multi_mutations = []
    positions = []
    aachanges = mutation_str.split(',')
    skipped_aachanges = 0
    if pos_set is None:
        pos_set = set()

    for aachange in aachanges:
        if aachange == '':
            skipped_aachanges += 1
            continue
        if aachange.count('delins') == 1:
            if len(aachanges) == 1:
                indel = indel_package.Substitution(raw_str=aachange, tags=tags)
            else:
                indel = indel_package.Substitution(raw_str=aachange)
            multi_mutations.append(indel)

        elif aachange.count('del') == 1:
            if len(aachanges) == 1:
                indel = indel_package.Deletion(raw_str=aachange, tags=tags)
            else:
                indel = indel_package.Deletion(raw_str=aachange)
            multi_mutations.append(indel)

        elif aachange.count('ins') == 1:
            if len(aachanges) == 1:
                indel = indel_package.Insertion(raw_str=aachange, tags=tags)
            else:
                indel = indel_package.Insertion(raw_str=aachange)
            multi_mutations.append(indel)

        else:
            indel = None
            try:
                aachange, aa1, aa2, pos = processAAChange(aachange, pdb_style=pdb_style)  # if pdb_style is True, then pos is actually pdb res_id
            except:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                config.errorlog.add_warning("Position/Mutation String Format Error: %s\n%s\n%s\n%s" % (mutation_str, str(e), str(f), str(g)))
                return None

            if pdb_style:
                if aa2 is None:
                    position = position_package.Position(pdb_res_nr=pos, wt_aa=aa1, tags=tags)
                    positions.append(position)
                else:
                    if len(aachanges) == 1:
                        position = position_package.Position(pdb_res_nr=pos, wt_aa=aa1, mut_aas=set(aa2), mut_tags_map={aa2: tags})
                    else:
                        position = position_package.Position(pdb_res_nr=pos, wt_aa=aa1, mut_aas=set(aa2))
                    multi_mutations.append((position, aa2))
            else:
                if aa2 is None:
                    if aachange not in pos_set:
                        position = position_package.Position(pos=pos, wt_aa=aa1, tags=tags)
                        positions.append(position)
                        pos_set.add(aachange)
                else:
                    if len(aachanges) == 1:
                        position = position_package.Position(pos=pos, wt_aa=aa1, mut_aas=set(aa2), mut_tags_map={aa2: tags})
                        positions.append(position)
                        pos_set.add(aachange)
                    elif aachange not in pos_set:
                        position = position_package.Position(pos=pos, wt_aa=aa1, mut_aas=set(aa2))
                        positions.append(position)
                        pos_set.add(aachange)
                    else:
                        position = pos
                    multi_mutations.append((position, aa2))
            

    if len(aachanges) == skipped_aachanges:
        protein_specific_tags = tags
    else:
        protein_specific_tags = None

    if len(multi_mutations) > 1:
        multi_mutation_tags = tags
    else:
        multi_mutation_tags = None

    return positions, multi_mutations, gene_id, protein_specific_tags, multi_mutation_tags, pos_set


def fuse_multi_mutations(list_of_multi_mutations):
    fused_snvs = {}
    indels_to_process = []
    for multi_mutation in list_of_multi_mutations:
        for pos in multi_mutation.snvs:
            fused_snvs[pos] = multi_mutation.snvs[pos]
        for indel in multi_mutation.indels:
            indels_to_process.append(indel)
    fused_deletions = []
    fused_insertions = []
    fused_substitutions = []
    while len(indels_to_process) > 0:
        #print(f'fused_deletions:')
        #for f_del in fused_deletions:
        #    print(f_del.get_notation())

        indel = indels_to_process.pop()

        #print(f'Indel to process: {indel.get_notation()}')

        indel_type, _, _ = indel.get_type()
        if indel_type == 'Deletion':
            fused_one = False
            already_in = False
            for f_del_pos, f_del in enumerate(fused_deletions):
                fuse = False
                if indel.left_flank_pos < f_del.left_flank_pos and indel.right_flank_pos > f_del.left_flank_pos:
                    fuse = True
                    new_left_flank = indel.left_flank
                else:
                    new_left_flank = f_del.left_flank
                if indel.right_flank_pos > f_del.right_flank_pos and indel.left_flank_pos < f_del.right_flank_pos:
                    fuse = True
                    new_right_flank = indel.right_flank
                else:
                    new_right_flank = f_del.right_flank

                if fuse:
                    new_f_del = indel_package.Deletion(left_flank = new_left_flank, right_flank = new_right_flank)
                    new_f_del.wt_prot = indel.wt_prot
                    fused_one = True
                    break
                if indel.left_flank_pos >= f_del.left_flank_pos and indel.right_flank_pos <= f_del.right_flank_pos:
                    already_in = True

            if fused_one:    
                fused_deletions.pop(f_del_pos)
                indels_to_process.append(new_f_del)
                continue

            del_list = []
            for f_ins_pos, f_ins in enumerate(fused_insertions):
                if f_ins.left_flank_pos >= indel.left_flank_pos and f_ins.right_flank_pos <= indel.right_flank_pos:
                    del_list.append(f_ins_pos)
            del_list = sorted(del_list, reverse = True)
            for pos in del_list:
                fused_insertions.pop(pos)

            fused_one = False
            for f_sub_pos, f_sub in enumerate(fused_substitutions):
                fuse = False
                if indel.left_flank_pos < f_sub.left_flank_pos and indel.right_flank_pos > f_sub.left_flank_pos:
                    fuse = True
                    new_left_flank = indel.left_flank
                else:
                    new_left_flank = f_sub.left_flank
                if indel.right_flank_pos > f_sub.right_flank_pos and indel.left_flank_pos < f_sub.right_flank_pos:
                    fuse = True
                    new_right_flank = indel.right_flank
                else:
                    new_right_flank = f_sub.right_flank
                if fuse:
                    new_f_sub = indel_package.Substitution(left_flank = new_left_flank, right_flank = new_right_flank, inserted_sequence = f_sub.inserted_sequence)
                    new_f_sub.wt_prot = indel.wt_prot
                    fused_one = True
                    break
            if fused_one:    
                fused_substitutions.pop(f_sub_pos)
                indels_to_process.append(new_f_sub)
                continue

            #This point gets only reached in the loop, if no fusion happened
            if not already_in:
                fused_deletions.append(indel)

        elif indel_type == 'Insertion':
            good_to_add = True

            for f_del in fused_insertions:
                if indel.left_flank_pos >= f_del.left_flank_pos and indel.right_flank_pos <= f_del.right_flank_pos:
                    good_to_add = False
                    break

            need_to_remove = False
            for f_ins_pos, f_ins in enumerate(fused_insertions):
                if f_ins.left_flank_pos == indel.left_flank_pos:
                    if len(indel.inserted_sequence) > len(f_ins.inserted_sequence):
                        need_to_remove = True
                    else:
                        good_to_add = False
                    break
            if need_to_remove:
                fused_insertions.pop(f_ins_pos)

            for f_sub in fused_substitutions:
                if indel.left_flank_pos >= f_sub.left_flank_pos and indel.right_flank_pos <= f_sub.right_flank_pos:
                    good_to_add = False
                    break
            if good_to_add:
                fused_insertions.append(indel)

        else:
            fused_one = False
            for f_del_pos, f_del in enumerate(fused_deletions):
                fuse = False
                if indel.left_flank_pos < f_del.left_flank_pos and indel.right_flank_pos > f_del.left_flank_pos:
                    fuse = True
                    new_left_flank = indel.left_flank
                else:
                    new_left_flank = f_del.left_flank
                if indel.right_flank_pos > f_del.right_flank_pos and indel.left_flank_pos < f_del.right_flank_pos:
                    fuse = True
                    new_right_flank = indel.right_flank
                else:
                    new_right_flank = f_del.right_flank
                if fuse:
                    new_f_sub = indel_package.Substitution(left_flank = new_left_flank, right_flank = new_right_flank, inserted_sequence = indel.inserted_sequence)
                    new_f_sub.wt_prot = indel.wt_prot
                    fused_one = True
                    break
            if fused_one:    
                fused_deletions.pop(f_del_pos)
                indels_to_process.append(new_f_sub)
                continue

            del_list = []
            for f_ins_pos, f_ins in enumerate(fused_insertions):
                if f_ins.left_flank_pos >= indel.left_flank_pos and f_ins.right_flank_pos <= indel.right_flank_pos:
                    del_list.append(f_ins_pos)
            del_list = sorted(del_list, reverse = True)
            for pos in del_list:
                fused_insertions.pop(pos)

            fused_one = False
            already_in = False
            for f_sub_pos, f_sub in enumerate(fused_substitutions):
                fuse = False
                if indel.left_flank_pos < f_sub.left_flank_pos and indel.right_flank_pos > f_sub.left_flank_pos:
                    fuse = True
                    new_left_flank = indel.left_flank
                else:
                    new_left_flank = f_sub.left_flank
                if indel.right_flank_pos > f_sub.right_flank_pos and indel.left_flank_pos < f_sub.right_flank_pos:
                    fuse = True
                    new_right_flank = indel.right_flank
                else:
                    new_right_flank = f_sub.right_flank
                if fuse:
                    if len(indel.inserted_sequence) > len(f_sub.inserted_sequence):
                        new_inserted_sequence = indel.inserted_sequence
                    else:
                        new_inserted_sequence = f_sub.inserted_sequence
                    new_f_sub = indel_package.Substitution(left_flank = new_left_flank, right_flank = new_right_flank, inserted_sequence = new_inserted_sequence)
                    new_f_sub.wt_prot = indel.wt_prot
                    fused_one = True
                    break
                if indel.left_flank_pos >= f_sub.left_flank_pos and indel.right_flank_pos <= f_sub.right_flank_pos:
                    already_in = True
            if fused_one:    
                fused_substitutions.pop(f_sub_pos)
                indels_to_process.append(new_f_sub)
                continue
            #This point gets only reached in the loop, if no fusion happened
            if not already_in:
                fused_substitutions.append(indel)

    fused_multi_mutations = []
    filtered_fused_snvs = {}
    for pos in fused_snvs:
        good_to_add = True
        for f_del in fused_deletions:
            if pos >= f_del.left_flank_pos and pos <= f_del.right_flank_pos:
                good_to_add = False
                break
        for f_sub in fused_substitutions:
            if pos >= f_sub.left_flank_pos and pos <= f_sub.right_flank_pos:
                good_to_add = False
                break       
        if good_to_add:
            filtered_fused_snvs[pos] =  fused_snvs[pos]
    for f_ins in fused_insertions:
        fused_multi_mutations.append(f_ins)

    #print(f'final fused_deletions:')
    #for f_del in fused_deletions:
    #    print(f_del.get_notation())

    for f_del in fused_deletions:
        fused_multi_mutations.append(f_del)

    for f_sub in fused_substitutions:
        fused_multi_mutations.append(f_sub)

    fused_mm_obj = mutation_package.MultiMutation(multi_mutation.wt_prot, f'{multi_mutation.wt_prot}_fused_multi_mutation', fused_multi_mutations)
    fused_mm_obj.snvs = filtered_fused_snvs

    return fused_mm_obj

def generate_multi_mutation_list(wt_seq, mut_seq, aligner_class = None):
    try:
        (aligned_wt_seq, aligned_mut_seq) = call_biopython_alignment(wt_seq, mut_seq, aligner_class = aligner_class)
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()

        return None, f'Alignment failed:\n{list(wt_seq)}\n{list(mut_seq)}\n{e}\n{f}\n{g}'
    variants = extract_variants(aligned_wt_seq, aligned_mut_seq)
    return 0, variants

def generate_multi_mutation(wt_seq, mut_seq, config, aligner_class = None):
    
    mutation_str_results = generate_multi_mutation_list(wt_seq, mut_seq, aligner_class = aligner_class)
    if mutation_str_results[0] is None:
        return mutation_str_results[1]

    positions, multi_mutation = process_variants(config, mutation_str_results[1])

    return positions, multi_mutation


def extract_variants(ref_align, query_align):
    """Idneitifes variants in the given pairwise sequence alignment
    Follows HGVS nomeclature (https://varnomen.hgvs.org/)
    Args:
      ref_align: Aligned target seuquence
      query_align: Aligned query sequence
    Returns: Variants list (Format as required by StructMAn tool)
    """
    # Construct a reference seq/alignment index for actual reference seq coordinate mapping
    ref_coor = ref_coordinates(ref_align)

    # Non-gap char's index
    ref_non_gap_chars_ind = [m.start() for m in re.finditer(r"[A-Z]", ref_align)]

    # Sequence Variant Nomenclature
    # General: Two variants separated by one or more chars should preferably
    # be described individually and not as a "delins"
    # Substitution: Compared to a reference seq, ONE char is replaced by ONE other char
    substitutions = set()

    # Insertion: One or more gaps in the reference seq
    insertions = set()

    # Deletion: One or more gaps in the query seq
    deletions = set()

    # Deletion/Insertion (delins)
    # 1. Changes involving two or more consecutive chars
    # 2. Two or more continuous substitutions
    delins = set()

    # Counters
    # Add index when encountered the respective variant and this will be used in checking for neighbours
    ins_ind = []
    del_ind = []
    subs_ind = []

    # Varposants collectposon
    for pos, _ in enumerate(ref_align):
        if ref_align[pos] != query_align[pos]:
            # Insertion
            if ref_align[pos] == "-":
                ins_ind.append(pos)
                if pos - 1 in del_ind or pos - 1 in subs_ind or pos - 1 in delins:
                    delins.update([pos - 1, pos])
                else:
                    insertions.add(pos)
            # Deletion
            elif query_align[pos] == "-":
                del_ind.append(pos)
                if pos - 1 in ins_ind or pos - 1 in subs_ind or pos - 1 in delins:
                    delins.update([pos - 1, pos])
                else:
                    deletions.add(pos)
            # Substitution
            else:
                subs_ind.append(pos)
                if pos - 1 in ins_ind or pos - 1 in del_ind or pos - 1 in subs_ind:
                    delins.update([pos - 1, pos])
                else:
                    substitutions.add(pos)

    # Group neighbors
    delins = find_neighbors(insertions, delins)
    delins = find_neighbors(deletions, delins)
    delins = find_neighbors(substitutions, delins)

    # Remove duplicates if already present in delins
    substitutions = substitutions - delins
    insertions = group_consecutives(insertions - delins)
    deletions = group_consecutives(deletions - delins)
    delins = group_consecutives(delins)

    # Mapping to ref coordinates and formatting as per StructMAn's input requirements
    variants = []

    if substitutions:
        for i in substitutions:
            variants.append(f"{ref_coor[i][1]}{ref_coor[i][0]}{query_align[i]}")

    if insertions:
        for i, j in insertions:
            # If gap(s) in the beginning
            if i == 0:
                variants.append(
                    f"{ref_coor[j+1][1]}{ref_coor[j+1][0]}_{ref_coor[j+1][1]}{ref_coor[j+1][0]}ins{query_align[i:j+1]}"
                )
            # If gap(s) in the end
            elif i == len(ref_coor) - 1 or j == len(ref_coor) - 1:
                variants.append(
                    f"{ref_coor[i-1][1]}{ref_coor[i-1][0]}_{ref_coor[i-1][1]}{ref_coor[i-1][0]}ins{query_align[i:j+1]}"
                )
            # If gap(s) anywhere else
            else:
                variants.append(
                    f"{ref_coor[i-1][1]}{ref_coor[i-1][0]}_{ref_coor[j+1][1]}{ref_coor[j+1][0]}ins{query_align[i:j+1]}"
                )

    if deletions:
        for i, j in deletions:
            variants.append(
                f"{ref_coor[i][1]}{ref_coor[i][0]}_{ref_coor[j][1]}{ref_coor[j][0]}del{ref_align[i:j+1]}"
            )

    if delins:
        for i, j in delins:
            if ref_coor[i][1] == "-" and ref_coor[j][1] == "-":
                new_i = get_neighbor(ref_non_gap_chars_ind, i, "upstream")
                new_j = get_neighbor(ref_non_gap_chars_ind, j, "downstream")
                variants.append(
                    f"{ref_coor[new_i][1]}{ref_coor[new_i][0]}_{ref_coor[new_j][1]}{ref_coor[new_j][0]}delins{query_align[i:j+1]}".replace(
                        "-", ""
                    )
                )
            elif ref_coor[i][1] != "-" and ref_coor[j][1] != "-":
                variants.append(
                    f"{ref_coor[i][1]}{ref_coor[i][0]}_{ref_coor[j][1]}{ref_coor[j][0]}delins{query_align[i:j+1]}".replace(
                        "-", ""
                    )
                )
            elif ref_coor[i][1] == "-" and ref_coor[j][1] != "-":
                new_i = get_neighbor(ref_non_gap_chars_ind, i, "upstream")
                variants.append(
                    f"{ref_coor[new_i][1]}{ref_coor[new_i][0]}_{ref_coor[j][1]}{ref_coor[j][0]}delins{query_align[i:j+1]}".replace(
                        "-", ""
                    )
                )
            elif ref_coor[i][1] != "-" and ref_coor[j][1] == "-":
                new_j = get_neighbor(ref_non_gap_chars_ind, j, "downstream")
                variants.append(
                    f"{ref_coor[i][1]}{ref_coor[i][0]}_{ref_coor[new_j][1]}{ref_coor[new_j][0]}delins{query_align[i:j+1]}".replace(
                        "-", ""
                    )
                )

    return variants


def ref_coordinates(ref_alignment):
    """Builds a reference coordinate from the reference alignment for faster lookups
    Args:
      ref_alignment: Reference sequence alignment (String)
    Returns: Reference coordinate dictionary (key=> ref alignment index, value=> list(ref coordinate, char))
    """
    ref_seq_index = {}
    coor_counter = 1

    for ind, char in enumerate(ref_alignment):
        if char != "-":
            ref_seq_index[ind] = [coor_counter, char]
            coor_counter += 1
        else:
            ref_seq_index[ind] = [None, char]
    return ref_seq_index


def get_neighbor(iterable, number, direction):
    """Get the closest value from the iterable for the number based on the said direction
    Assumes iterable is sorted
    Args:
      iterable: A list of numbers to find the closest value in
      number: A number to which the closest value in the iterable needs to be found
      direction: Upstream or Downstream [upstream or downstream]
    Returns: The closest value of the number in the iterable
    """
    pos = bisect_left(iterable, number)

    if pos == 0:
        return iterable[0]

    if pos == len(iterable):
        return iterable[-1]

    before = iterable[pos - 1]
    after = iterable[pos]

    if direction == "upstream":
        return after
    else:
        return before


def group_consecutives(iterable):
    """Groups consecutive numbers
    Args:
      iterable: A set/list object with numbers
    Returns: A list object containing tuples with start and end values of consecutive numbers
    """
    if not iterable:
        return []

    ranges = []

    for group in consecutive_groups(sorted(iterable)):
        group = list(group)
        if len(group) == 1:
            ranges.append((group[0], group[0]))
        else:
            ranges.append((group[0], group[-1]))

    return ranges


def find_neighbors(iterable1, iterable2):
    """Finds upstream and downstream neighbors of items in iterable 2 in iterable 1
    Args:
      iterable1: A list or set of numbers
      iterable2: A list or set of numbers
    Returns: A set with all neighbors
    """

    neighbors = list(iterable2)

    # Donstream neighbors
    for i in neighbors:
        if i - 1 in iterable1:
            neighbors.append(i - 1)

    # Upstream neighbors
    for i in neighbors:
        if i + 1 in iterable1:
            neighbors.append(i + 1)

    return set(neighbors)


