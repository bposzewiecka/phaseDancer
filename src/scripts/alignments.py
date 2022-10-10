from collections import Counter

import numpy as np
import pysam

from src.scripts.constants import (
    BAM_CDEL,
    BAM_CDIFF,
    BAM_CEQUAL,
    BAM_CINS,
    BAM_CMATCH,
    BAM_CSOFT_CLIP,
    DEFAULT_NUMBER,
    DELETION_NUMBER,
    INSERTION_NUMBER,
    LETTERS_TO_NUMBERS,
    NON_INSERTION_NUMBER,
    get_excluded_values,
)
from src.scripts.defaults import DISTANCE_TO_END_THRESHOLD, MAX_NUMBER_OF_COORDS
from src.utils.reads_utils import get_sequence


def transform_sequence_to_letters(sequence):

    return [LETTERS_TO_NUMBERS[base] for base in sequence]


def get_alignment_on_reference(read, reference_size, masked_low_complexity):

    read_pos = 0
    ref_pos = read.reference_start

    alignment = [DEFAULT_NUMBER] * reference_size

    query_seq = transform_sequence_to_letters(read.query_sequence)

    insertions = [NON_INSERTION_NUMBER] * reference_size

    for operation, op_size in read.cigartuples:

        if operation in [BAM_CEQUAL, BAM_CDIFF, BAM_CMATCH]:

            alignment[ref_pos : ref_pos + op_size] = query_seq[
                read_pos : read_pos + op_size
            ]

            read_pos += op_size
            ref_pos += op_size

        elif operation in [BAM_CINS, BAM_CSOFT_CLIP]:

            if operation == BAM_CINS and ref_pos not in masked_low_complexity:
                insertions[ref_pos] = INSERTION_NUMBER

            read_pos += op_size

        elif operation in [BAM_CDEL]:

            alignment[ref_pos : ref_pos + op_size] = [DELETION_NUMBER] * op_size

            for pos in range(ref_pos, ref_pos + op_size):
                if pos in masked_low_complexity:
                    alignment[pos] = DEFAULT_NUMBER

            ref_pos += op_size

        else:
            raise Exception(f"Operation {operation} of size {op_size} in cigartuples.")

    return alignment, insertions


def get_second_most_common_freqs(arr, technology):

    excluded_values = get_excluded_values(technology)

    def smc_freq(bases):
        try:
            counter = Counter(bases).most_common()
            return [freq for base, freq in counter if base not in excluded_values][1]
        except IndexError:
            return 0

    return np.apply_along_axis(smc_freq, axis=0, arr=arr)


def get_most_common_bases(arr):
    def mc_base(bases):

        try:
            return Counter(bases).most_common()[0][0]
        except IndexError:
            return 0

    return np.apply_along_axis(mc_base, axis=0, arr=arr)


def get_consensus_seq(clusters, alignment_array):

    return [get_most_common_bases(alignment_array.arr[cluster]) for cluster in clusters]


def mask_low_complexity(reference, min_low_complexity_size):

    masked = set()

    last = reference[:2]
    count = 2

    for i, base in enumerate(reference):

        if i < 2:
            continue

        if base == last[0]:
            count += 1
        else:
            if count >= min_low_complexity_size:
                masked.update(set(list(range(i - count, i))))

            count = 2

        last = [last[1], base]

    print(masked)
    return masked


class AlignmentArray:
    def __init__(self, bam_fn, fasta_fn, technology):

        self.reference_name, self.reference_seq = get_sequence(fasta_fn)
        self.reference_seq = transform_sequence_to_letters(self.reference_seq)
        self.reference_size = len(self.reference_seq)
        self.technology = technology
        self.init_array(bam_fn)

    def init_array(self, bam_fn):

        self.masked_low_complexity = mask_low_complexity(self.reference_seq, 20)

        self.reads = []
        self.arr = []

        with pysam.AlignmentFile(bam_fn) as bam:  # pylint: disable=no-member

            i = 0

            for read in bam.fetch():

                if read.is_supplementary:
                    continue

                if self.reference_size - read.reference_end > DISTANCE_TO_END_THRESHOLD:
                    continue

                (
                    alignment_on_reference,
                    insertions_on_reference,
                ) = get_alignment_on_reference(
                    read, self.reference_size, self.masked_low_complexity
                )

                alignment = alignment_on_reference

                if self.technology == "hifi":
                    alignment += insertions_on_reference

                self.arr.append(alignment)

                self.reads.append(read)

                i += 1

        self.arr = np.vstack(self.arr)

    def get_read_names(self, reads):

        return [read.query_name for i, read in enumerate(self.reads) if i in reads]

    def get_reference_size(self):

        return self.reference_size

    def get_coords(self, bases_freq_threshold, partition):

        second_most_common_freqs = get_second_most_common_freqs(
            self.arr[partition.reads], self.technology
        )

        threshold = sorted(second_most_common_freqs)[-MAX_NUMBER_OF_COORDS]

        return [
            i
            for i, freq in enumerate(second_most_common_freqs)
            if freq >= threshold and freq >= bases_freq_threshold
        ]
