import numpy as np
import copy
import itertools
import multiprocessing as mp
from math import ceil
from classless_methods import generate_comparison_matrix, calculate_degeneracy, disambiguate, reverse_complement
from Primer import *
from class_methods import hamming_distance, levenshtein_distance


class PrimerIndex:
    thresholds = {
        'gc_lb': 0.4,
        'gc_ub': 0.6,
        'melting_lb': 55.,
        'melting_ub': 75.,
        'end_at_threshold': 2,
        'end_gc_threshold': 3,
        'monorun_threshold': 3,
        'duorun_threshold': 3,
        'mfe_threshold': -5.,
        'self_complementarity_threshold': 10
    }

    def __init__(self):
        self.primer2index = {'forward': {},
                             'reverse': {}}  # contains primer sequences as keys and corresponding primer index as value
        self.index2primer = {'forward': np.empty((0)), 'reverse': np.empty((0))}  # contains Primer objects in an array
        self.conflict_matrix = None
        self.comparison_matrix = generate_comparison_matrix()

        # Used to determine whether inexact matching should be used
        self.inexact = False
        self.mismatches = 0
        # It is important to note that the similarity score between the two forward primers, will be the same
        # as the similarity score between their respective reverse primers.
        self.similar_primers = {'forward': {},
                                'reverse': {}}  # Contains clusters of similar primers (Using Hamming Distance)

    # This does not work as intended currently
    def __eq__(self, other):
        try:
            for orientation in self.set:
                for p in self.set[orientation]:
                    if not self.set[orientation][p].indices == other.set[orientation][p].indices or \
                            self.set[orientation][p].feasible == other.set[orientation][p].feasible:
                        return False
            return True

            for orientation in other.set:
                for p in other.set[orientation]:
                    if not self.set[orientation][p].indices == other.set[orientation][p].indices or \
                            self.set[orientation][p].feasible == other.set[orientation][p].feasible:
                        return False
            return True
        except:
            return False

    def add_primer(self, sequence, orientation):
        '''
        Function that adds a new primer to the primer index

        Parameters
        ----------
        sequence : str
            String representation of the primer to add.
        orientation : str
            'forward' or 'reverse' indicating the primer orientation.

        Returns
        -------
        None.

        '''
        if sequence not in self.primer2index[orientation]:
            self.primer2index[orientation][sequence] = len(self.index2primer[orientation])
            self.index2primer[orientation] = np.append(self.index2primer[orientation], Primer(sequence, orientation))
            # self.index2primer[orientation].append(Primer(sequence, orientation))
            self.index2primer[orientation][-1].check_feasibility(self.comparison_matrix,
                                                                 gc_lb=self.thresholds['gc_lb'],
                                                                 gc_ub=self.thresholds['gc_ub'],
                                                                 melting_lb=self.thresholds['melting_lb'],
                                                                 melting_ub=self.thresholds['melting_ub'],
                                                                 end_at_threshold=self.thresholds['end_at_threshold'],
                                                                 end_gc_threshold=self.thresholds['end_gc_threshold'],
                                                                 monorun_threshold=self.thresholds['monorun_threshold'],
                                                                 duorun_threshold=self.thresholds['duorun_threshold'],
                                                                 mfe_threshold=self.thresholds['mfe_threshold'],
                                                                 self_complementarity_threshold=self.thresholds[
                                                                     'self_complementarity_threshold'])

    def add_sequence(self, sequence, sequence_index, primer_sequence, orientation):
        '''
        Function that adds a sequence to a primer in the index, or create a new primer and add the sequence to it

        Parameters
        ----------
        sequence : Sequence
            Sequence object that should be related to a primer.
        sequence_index : str
            Starting index of the primer that has to be linked to this sequence.
        primer_sequence : str
            String representation of the primer to add.
        orientation : str
            'forward' or 'reverse' indicating the primer orientation.

        Returns
        -------
        bool
            True if the sequence has been successfully added, False otherwise

        '''
        if primer_sequence not in self.primer2index[orientation]:
            self.add_primer(primer_sequence, orientation)
            self.index2primer[orientation][-1].add_sequence(sequence, sequence_index)
        else:
            self.index2primer[orientation][self.primer2index[orientation][primer_sequence]].add_sequence(sequence,
                                                                                                         sequence_index)
        return self.index2primer[orientation][-1].feasible

    def remove_redundant(self):
        '''
        Function that removes all of the primers in this PrimerIndex that are infeasible. Note that the function itself
        does not check feasibility, instead see Primer.check_feasibility!

        Returns
        -------
        None.

        '''
        common_kmers = set(self.primer2index['forward'].keys()).intersection(set(self.primer2index['reverse'].keys()))
        for orientation in self.primer2index:
            kmers = list(self.primer2index[orientation].keys())
            print('Initially contains %d %s primers' % (len(kmers), orientation))
            index = 0
            to_remove = []
            while index < len(kmers):
                if not self.index2primer[orientation][index].feasible or kmers[index] in common_kmers:
                    to_remove.append(index)
                    self.primer2index[orientation].pop(kmers[index])
                index += 1
            self.index2primer[orientation] = np.delete(self.index2primer[orientation], to_remove)
            for index in range(len(self.index2primer[orientation])):
                self.primer2index[orientation][self.index2primer[orientation][index].sequence] = index
            print('Finally contains %d %s primers' % (len(self.primer2index[orientation]), orientation))
            print('Removed %d primers occurring both as forward and reverse' % (len(common_kmers)))

    @staticmethod
    def set_thresholds(thresholds):
        '''
        Function that sets the primer property thresholds to the given values. Note that this does not check for
        existing primers in the index whether they satisfy the new thresholds and thus should be set beforehand.

        Parameters
        ----------
        thresholds : dict[ String ]
            Dictionary containing the properties as keys, and the values they should be set to as values.

        Returns
        -------
        None.

        '''
        for prop in thresholds:
            try:
                PrimerIndex.thresholds[prop] = thresholds[prop]
            except:
                continue

    def merge_indices(self, other_index):
        for orientation in other_index.primer2index:
            primers_to_add = []
            indices_to_add = []
            for primer in other_index.primer2index[orientation]:
                index_in_other = other_index.primer2index[orientation][primer]
                # Check if primer is also in this index
                if primer in self.primer2index[orientation]:
                    index_in_this = self.primer2index[orientation][primer]
                    # Check if primer is feasible both in this index and other index
                    if self.index2primer[orientation][index_in_this].feasible and other_index.index2primer[orientation][
                        index_in_other].feasible:
                        for sequence in other_index.index2primer[orientation][index_in_other].indices:
                            # Check if the sequence is also linked to the primer in this index
                            if sequence in self.index2primer[orientation][
                                self.primer2index[orientation][primer]].indices:
                                # Check if primer occurs at the same indices, otherwise it should be rejected
                                if self.index2primer[orientation][index_in_this].indices[sequence] == \
                                        other_index.index2primer[orientation][index_in_other].indices[sequence]:
                                    continue
                                else:
                                    self.index2primer[orientation][index_in_this].feasible = False
                            # If sequence is not linked to the primer in this, add it
                            else:
                                self.index2primer[orientation][index_in_this].indices[sequence] = copy.deepcopy(
                                    other_index.index2primer[orientation][index_in_other].indices[sequence])
                    # If primer is infeasible in either index then set to infeasible
                    else:
                        self.index2primer[orientation][index_in_this].feasible = False
                # If primer is not in this index add it and copy information from other index
                else:
                    primers_to_add.append(primer)
                    indices_to_add.append(index_in_other)
            k = 0
            for primer in primers_to_add:
                self.primer2index[orientation][primer] = len(self.index2primer[orientation]) + k
                k += 1
            self.index2primer[orientation] = np.append(self.index2primer[orientation],
                                                       other_index.index2primer[orientation][indices_to_add])

    def check_amplicon(self, sequences, amplicon, primer_width, search_width):
        '''
        Function that generates the primers (per sequence) of length $primer_width in a search window of length
        $search_width that can be used to amplify this amplicon for all the sequences in $sequences.

        Parameters
        ----------
        sequences : list[ Sequence ]
            List of sequences to generate primers for.
        amplicon : Amplicon
            Amplicon to find primers around.
        primer_width : int
            Width of primers in number of nucleotides.
        search_width : int
            Search window around the amplicon in which we want to find primers.

        Returns
        -------
        None
        
        '''
        sequence_ids = [sequence.id_num for sequence in sequences]
        amplicon.primers = {'forward': {s: set() for s in sequence_ids}, 'reverse': {s: set() for s in sequence_ids}}
        amplicon.full_primerset = {'forward': set(), 'reverse': set()}

        for sequence in sequences:
            # Check if the start of the amplicon is a misalign, in which case correct for it
            if sequence.aligned_to_trim[amplicon.start] == sequence.aligned_to_trim[amplicon.start - 1]:
                forward_end_index = sequence.aligned_to_trim[amplicon.start] + 1
            else:
                forward_end_index = sequence.aligned_to_trim[amplicon.start]
            # Check if the character after the amplicon is a misalign, in which case correct for it
            if sequence.aligned_to_trim[amplicon.end] == sequence.aligned_to_trim[amplicon.end - 1]:
                reverse_start_index = sequence.aligned_to_trim[amplicon.end] + 1
            else:
                reverse_start_index = sequence.aligned_to_trim[amplicon.end]

            # Check which primers (both forward and reverse) are found for the corresponding sequences
            # this iterates over the possible primers within the search range
            # TODO Get the index where the primer would bind to in the amplicons sequence? Figure out what `start` and
            #  `end` represents in the Amplicon class
            for offset in range(search_width - primer_width + 1):
                # Iterate over forward primers for this sequence
                current_fwd_primer = \
                    sequence.sequence_raw[forward_end_index - primer_width - offset: forward_end_index - offset]
                # only proceed if the sequence is not "too degenerate"
                if calculate_degeneracy(current_fwd_primer) <= 4 ** 5:
                    for forward_primer in disambiguate(current_fwd_primer):
                        if forward_primer in self.primer2index['forward']:
                            if self.index2primer['forward'][self.primer2index['forward'][forward_primer]].feasible:
                                amplicon.primers['forward'][sequence.id_num].add(self.primer2index['forward'][forward_primer])
                                amplicon.full_primerset['forward'].add(self.primer2index['forward'][forward_primer])
                            # !!! - SHOULD ONLY BE EXECUTED WHEN self.inexact IS TRUE
                            if self.inexact:  # If inexact matching is enabled, we also add all the similar primers
                                self.add_similar_primers(amplicon, sequence, 'forward', forward_primer)

                # Iterate over reverse primers for this sequence
                current_rev_primer = reverse_complement(
                    sequence.sequence_raw[reverse_start_index + offset: reverse_start_index + primer_width + offset])
                if calculate_degeneracy(current_rev_primer) <= 4 ** 5:
                    for reverse_primer in disambiguate(current_rev_primer):
                        if reverse_primer in self.primer2index['reverse']:
                            if self.index2primer['reverse'][self.primer2index['reverse'][reverse_primer]].feasible:
                                amplicon.primers['reverse'][sequence.id_num].add(self.primer2index['reverse'][reverse_primer])
                                amplicon.full_primerset['reverse'].add(self.primer2index['reverse'][reverse_primer])
                            if self.inexact:  # If inexact matching is enabled, we also add all the similar primers
                                self.add_similar_primers(amplicon, sequence, 'reverse', reverse_primer)

    def update_conflict_matrix(self, primers):
        '''
        Function that generates the conflicts for all the primer pairs that can be obtained by taking combinations of primers from $primers. If this PrimerIndex
        already has a conflict matrix, it will only be updated and not generated again.

        Parameters
        ----------
        primers : list[ Primer ]
            List of Primer objects to determine conflicts between.

        Returns
        -------
        None.

        '''
        if not self.conflict_matrix:
            # Matrix entry will be equal to -1 if not yet assigned, 1 if primers have a conflict, 2 if primers don't
            # have a conflict
            self.conflict_matrix = {
                ('f', 'r'): -1 * np.ones((len(self.index2primer['forward']), len(self.index2primer['reverse'])),
                                         dtype=np.int8),
                ('f', 'f'): -1 * np.ones((len(self.index2primer['forward']), len(self.index2primer['forward'])),
                                         dtype=np.int8),
                ('r', 'r'): -1 * np.ones((len(self.index2primer['reverse']), len(self.index2primer['reverse'])),
                                         dtype=np.int8)}
        # Iterate over primer pairs
        for pair in itertools.combinations(primers, 2):
            # First primer is forward, second is reverse
            if pair[0].orientation == 'forward' and pair[1].orientation == 'reverse':
                current_index_pair = (
                    self.primer2index['forward'][pair[0].sequence], self.primer2index['reverse'][pair[1].sequence])
                if self.conflict_matrix[('f', 'r')][current_index_pair] == -1:
                    if pair[0].check_compatibility(pair[1], self.comparison_matrix,
                                                   self.thresholds['self_complementarity_threshold'])[0] > \
                            self.thresholds['self_complementarity_threshold']:
                        self.conflict_matrix[('f', 'r')][current_index_pair] = 1
                    else:
                        self.conflict_matrix[('f', 'r')][current_index_pair] = 2
            # First primer is reverse, second is forward
            elif pair[0].orientation == 'reverse' and pair[1].orientation == 'forward':
                current_index_pair = (
                    self.primer2index['forward'][pair[1].sequence], self.primer2index['reverse'][pair[0].sequence])
                if self.conflict_matrix[('f', 'r')][current_index_pair] == -1:
                    if pair[0].check_compatibility(pair[1], self.comparison_matrix,
                                                   self.thresholds['self_complementarity_threshold'])[0] > \
                            self.thresholds['self_complementarity_threshold']:
                        self.conflict_matrix[('f', 'r')][current_index_pair] = 1
                    else:
                        self.conflict_matrix[('f', 'r')][current_index_pair] = 2
            # Both primers are forward
            elif pair[0].orientation == 'forward' and pair[1].orientation == 'forward':
                current_index_pair = (
                    self.primer2index['forward'][pair[0].sequence], self.primer2index['forward'][pair[1].sequence])
                if self.conflict_matrix[('f', 'f')][current_index_pair] == -1:
                    if pair[0].check_compatibility(pair[1], self.comparison_matrix,
                                                   self.thresholds['self_complementarity_threshold'])[0] > \
                            self.thresholds['self_complementarity_threshold']:
                        self.conflict_matrix[('f', 'f')][current_index_pair] = 1
                        self.conflict_matrix[('f', 'f')][current_index_pair[1], current_index_pair[0]] = 1
                    else:
                        self.conflict_matrix[('f', 'f')][current_index_pair] = 2
                        self.conflict_matrix[('f', 'f')][current_index_pair[1], current_index_pair[0]] = 2
            # Both primers are reverse
            else:
                current_index_pair = (
                    self.primer2index['reverse'][pair[0].sequence], self.primer2index['reverse'][pair[1].sequence])
                if self.conflict_matrix[('r', 'r')][current_index_pair] == -1:
                    if pair[0].check_compatibility(pair[1], self.comparison_matrix,
                                                   self.thresholds['self_complementarity_threshold'])[0] > \
                            self.thresholds['self_complementarity_threshold']:
                        self.conflict_matrix[('r', 'r')][current_index_pair] = 1
                        self.conflict_matrix[('r', 'r')][current_index_pair[1], current_index_pair[0]] = 1
                    else:
                        self.conflict_matrix[('r', 'r')][current_index_pair] = 2
                        self.conflict_matrix[('r', 'r')][current_index_pair[1], current_index_pair[0]] = 2

    def check_conflict(self, primer_pair):
        self.update_conflict_matrix(primer_pair)
        if primer_pair[0].orientation == 'forward' and primer_pair[1].orientation == 'reverse':
            orientation = ('f', 'r')
            pair = (self.primer2index['forward'][primer_pair[0].sequence],
                    self.primer2index['reverse'][primer_pair[1].sequence])
        elif primer_pair[1].orientation == 'reverse' and primer_pair[1].orientation == 'forward':
            orientation = ('f', 'r')
            pair = (self.primer2index['forward'][primer_pair[1].sequence],
                    self.primer2index['reverse'][primer_pair[0].sequence])
        else:
            orientation = (primer_pair[0].orientation[0], primer_pair[1].orientation[0])
            pair = (self.primer2index[primer_pair[0].orientation][primer_pair[0].sequence],
                    self.primer2index[primer_pair[1].orientation][primer_pair[1].sequence])
        return self.conflict_matrix[orientation][pair]

    @staticmethod
    def generate_primer(primer_index, sequence, width, max_degeneracy=4 ** 5):
        for cur_index in range(sequence.length_raw - width + 1):
            current_fwd_primer = sequence.sequence_raw[cur_index: cur_index + width]
            if calculate_degeneracy(current_fwd_primer) <= max_degeneracy:
                for forward_primer in disambiguate(current_fwd_primer):
                    primer_index.add_sequence(sequence, cur_index, forward_primer, 'forward')
            current_rev_primer = reverse_complement(current_fwd_primer)
            if calculate_degeneracy(current_rev_primer) <= max_degeneracy:
                for reverse_primer in disambiguate(current_rev_primer):
                    primer_index.add_sequence(sequence, cur_index, reverse_primer, 'reverse')

    @staticmethod
    def generate_index(sequences, width, comparison_matrix, max_degeneracy=4 ** 5, mismatches=0):
        '''
        Static function that generates a primer index for the given sequences using a primer width of $width. For the
        multiprocessing variant see generate_index_mp

        Parameters
        ----------
        mismatches : Number of mismatches allowed if inexact matching is enabled.
        sequences : list[ Sequence ]
            List of sequences to find primers in.
        width : int
            Width of the primers to include in this index.
        comparison_matrix : dict[ (char,char) ]
            Dictionary that determines which characters should be considered equal.
        max_degeneracy : int, optional
            Maximum allowed degeneracy of a k-mer for processing it.

        Returns
        -------
        primer_index : PrimerIndex
            Primer index containing all the primers of length #width that appear in the sequences in $sequences.

        '''
        primer_index = PrimerIndex()
        primer_index.mismatches = max(0, mismatches)
        primer_index.inexact = primer_index.mismatches > 0
        i = 0
        if type(sequences) == list:  # If multiple sequences supplied
            for sequence in sequences:
                PrimerIndex.generate_primer(primer_index, sequence, width, max_degeneracy)
                # for cur_index in range(sequence.length_raw - width + 1):
                #     current_fwd_primer = sequence.sequence_raw[cur_index : cur_index + width]
                #     if calculate_degeneracy(current_fwd_primer) <= max_degeneracy:
                #         for forward_primer in disambiguate(current_fwd_primer):
                #             primer_index.add_sequence(sequence, cur_index, forward_primer, 'forward')
                #     current_rev_primer = reverse_complement(current_fwd_primer)
                #     if calculate_degeneracy(current_rev_primer) <= max_degeneracy:
                #         for reverse_primer in disambiguate(current_rev_primer):
                #             primer_index.add_sequence(sequence, cur_index, reverse_primer, 'reverse')
                i += 1
        else:  # If a singular sequence is supplied
            PrimerIndex.generate_primer(primer_index, sequences, width, max_degeneracy)
            # for cur_index in range(sequences.length_raw - width + 1):
            #     current_fwd_primer = sequences.sequence_raw[cur_index : cur_index + width]
            #     if calculate_degeneracy(current_fwd_primer) <= max_degeneracy:
            #         for forward_primer in disambiguate(current_fwd_primer):
            #             primer_index.add_sequence(sequences, cur_index, forward_primer, 'forward')
            #     current_rev_primer = reverse_complement(current_fwd_primer)
            #     if calculate_degeneracy(current_rev_primer) <= max_degeneracy:
            #         for reverse_primer in disambiguate(current_rev_primer):
            #             primer_index.add_sequence(sequences, cur_index, reverse_primer, 'reverse')
            i += 1
        return primer_index

    @staticmethod
    def generate_index_mp(sequences, width, comparison_matrix, max_degeneracy=4 ** 5, processors=1, mismatches=0):
        if processors > 1:
            sequences_partitioned = [sequences[i:i + (ceil(len(sequences) / processors))] for i in
                                     range(0, len(sequences), ceil(len(sequences) / processors))]
            with mp.Pool(processors) as pool:
                indices = pool.starmap(PrimerIndex.generate_index, zip(sequences_partitioned, itertools.repeat(width),
                                                                       itertools.repeat(comparison_matrix),
                                                                       itertools.repeat(max_degeneracy)))
            master_index = indices[0]
            for index in indices[1:]:
                master_index.merge_indices(index)
            return master_index
        else:
            return PrimerIndex.generate_index(sequences, width, comparison_matrix, mismatches)

    def primer_similarity(self, similarity_metric='hd'):
        '''
        Function is used to iterate over all the primers and compute their similarity score, and update them in the
            self.similar_primers[orientation][primer] variable

        Parameters
        ----------
        e - Number of mismatches allowed between the primers

        Returns
        -------

        '''
        print('Primer similarity function accessed')
        for primer in self.index2primer['forward']:
            for similar_primer in self.index2primer['forward']:
                if primer.__eq__(similar_primer):
                    continue
                if similarity_metric == 'hd':
                    similarity = hamming_distance(primer, similar_primer)
                else:
                    similarity = levenshtein_distance(primer, similar_primer)
                # The similarity score between two forward primers will be the same as the similarity between
                # their complements

                # TODO the similarity should not be stored in the Primer instance, since the Primer instance can be
                #  removed when the `remove_redundant()` function is called, instead it can be stored in the PrimerIndex

                # If they are similar then add the primers to each other's respective similarity set
                if similarity <= self.mismatches:
                    # Forward Primers
                    self.similar_primers['forward'][primer].add(similar_primer)
                    self.similar_primers['forward'][similar_primer].add(primer)

                    # Reverse primers
                    reverse_primer = Primer(reverse_complement(primer.sequence), 'reverse')
                    reverse_similar_primer = Primer(reverse_complement(similar_primer.seqence), 'reverse')

                    self.similar_primers['reverse'][reverse_primer].add(reverse_similar_primer)
                    self.similar_primers['reverse'][reverse_similar_primer].add(reverse_primer)

    def add_similar_primers(self, amplicon, sequence, orientation, original_primer):
        '''
        This function iterates over all the primers which are considered to be similar to the `original_primer` and
        adds them to the set of primers which should be considered for this `amplicon`.

        Parameters
        ----------
        amplicon - Amplicon for which we are adding the primers
        sequence - The sequence being considered
        orientation - The orientation of the primers (forward or reverse)
        original_primer - Primer being considered, for which we should also include the similar primers

        Returns
        -------
        '''
        # TODO here we should keep track of the indices where these inexact primers will bind to the
        #  amplicon, since inexact matching is allowed, we need to ensure that the same primer is
        #  NOT used on the same amplicon within a range of 2 * search_window !!!
        for similar_primer in self.similar_primers[orientation][original_primer]:
            if similar_primer in self.primer2index[orientation]:
                if self.index2primer[orientation][self.primer2index[orientation][similar_primer]].feasible:
                    amplicon.primers[orientation][sequence.id_num].add(self.primer2index[orientation][similar_primer])
                    amplicon.full_primerset[orientation].add(self.primer2index[orientation][similar_primer])