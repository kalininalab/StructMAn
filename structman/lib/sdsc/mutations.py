import math

from structman.lib.sdsc import position as position_package
from structman.lib.sdsc.indel import Indel
from structman.lib.sdsc.sdsc_utils import doomsday_protocol, Slotted_obj


class MultiMutation(Slotted_obj):
    __slots__ = ['wt_prot', 'mut_prot', 'snvs', 'indels', 'database_id', 'stored', 'tags']

    def __init__(self, wt_prot, mut_prot, mutation_list, tags=None):

        self.wt_prot = wt_prot
        self.mut_prot = mut_prot
        self.snvs = {}
        self.indels = []
        if isinstance(tags, str):
            tags = set(tags.split(','))
        if tags is not None:
            self.tags = tags.copy()
        else:
            self.tags = None
        for mut in mutation_list:
            if not isinstance(mut, tuple):
                self.indels.append(mut)
            else:
                self.snvs[mut[0].pos] = (mut[0].mut_aas[mut[1]])
        self.database_id = None
        self.stored = False

    def deconstruct(self):
       for snv in self.snvs:
           self.snvs[snv].deconstruct()
       del self.snvs
       for indel in self.indels:
           self.indels[indel].deconstruct()
       del self.indels
       del self.tags
       doomsday_protocol(self)

    def get_snv_db_ids(self):
        db_ids = []
        for pos in self.snvs:
            snv = self.snvs[pos]
            db_ids.append(snv.database_id)
        return db_ids

    def get_indel_db_ids(self):
        db_ids = []
        for indel in self.indels:
            db_ids.append(indel.database_id)
        return db_ids

    def print_indels(self):
        for indel in self.indels:
            print(indel.get_notation())

    def get_score(self, proteins):
        squared_pos_score_sum = 0.
        for pos in self.snvs:
            pos_score = proteins[self.wt_prot].positions[pos].get_score()
            squared_pos_score_sum += pos_score**2

        #print(f'Squared pos score sum before indels: {squared_pos_score_sum}')

        ins_score_sum = 0.
        for indel in self.indels:
            indel_type, terminal, terminal_end = indel.get_type()
            if indel_type == 'Insertion':
                ins_score = indel.get_score(proteins)
                ins_score_sum += ins_score
            elif indel_type == 'Deletion':
                pos_scores = indel.get_position_scores(proteins)
                for score in pos_scores:
                    squared_pos_score_sum += score**2
            else:
                pos_scores, insertion_score = indel.get_scores(proteins)
                for score in pos_scores:
                    squared_pos_score_sum += score**2
                ins_score_sum += insertion_score

        #print(f'Squared pos score sum: {squared_pos_score_sum}, ins score sum: {ins_score_sum}')

        multi_mutation_score = math.sqrt(squared_pos_score_sum) + ins_score_sum
        return multi_mutation_score

    def mutate(self, proteins, config):
        if self.mut_prot is None:
            return
        if self.mut_prot not in proteins.protein_map:
            return
        if proteins[self.mut_prot].sequence is None:
            if proteins[self.wt_prot].sequence is None:
                config.errorlog.add_error('Sequence is None: %s %s' % (self.wt_prot, self.mut_prot))
                return
            mut_seq = proteins[self.wt_prot].sequence
            indel_position_order = {}
            for pos in self.snvs:
                wt_aa = proteins[self.wt_prot].positions[pos].wt_aa
                aa2 = self.snvs[pos].new_aa
                mut_seq = '%s%s%s' % (mut_seq[:pos], aa2, mut_seq[(pos + 1):])

            for indel in self.indels:
                indel_position_order[indel.left_flank_pos] = indel

            pos_sorted = sorted(indel_position_order.keys(), reverse=True)
            for pos in pos_sorted:
                mut_seq = indel_position_order[pos].return_mutate_sequence(mut_seq)
            proteins[self.mut_prot].sequence = mut_seq

            for (pos, aa) in enumerate(mut_seq):
                seq_pos = pos + 1
                position = position_package.Position(pos=seq_pos, wt_aa=aa, checked=True)
                proteins[self.mut_prot].add_positions([position])

        else:
            mut_seq = proteins[self.mut_prot].sequence

    def get_affected_positions(self):
        insertions = []
        deletions = []
        subsitutions = []
        for pos in self.snvs:
            subsitutions.append((pos, self.snvs[pos].new_aa))
        for indel in self.indels:
            ins, dels, subs = indel.get_affected_positions()
            insertions += ins
            deletions += dels
            subsitutions += subs

        return insertions, deletions, subsitutions
    
    def give_count(self):
        insertions, deletions, subsitutions = self.get_affected_positions()
        l_inserted = 0

        for _, inserted_seq in insertions:
            l_inserted += len(inserted_seq)
            
        return len(deletions), l_inserted, len(subsitutions)