from collections import Counter, defaultdict

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
    LETTERS_TO_NUMBERS,
)
from src.scripts.defaults import DISTANCE_TO_END_THRESHOLD, MAX_NUMBER_OF_COORDS
from src.utils.reads_utils import get_sequence


def transform_sequence_to_letters(sequence):

    return [LETTERS_TO_NUMBERS[base] for base in sequence]


def get_alignment_on_reference(read, reference_size):

    read_pos = 0
    ref_pos = read.reference_start

    alignment = [DEFAULT_NUMBER] * reference_size

    query_seq = transform_sequence_to_letters(read.query_sequence)

    deletions = []
    insertions = []

    for operation, op_size in read.cigartuples:

        if operation in [BAM_CEQUAL, BAM_CDIFF, BAM_CMATCH]:

            alignment[ref_pos : ref_pos + op_size] = query_seq[
                read_pos : read_pos + op_size
            ]

            read_pos += op_size
            ref_pos += op_size

        elif operation in [BAM_CINS, BAM_CSOFT_CLIP]:

            if operation == BAM_CINS:
                insertions.append((ref_pos, query_seq[read_pos : read_pos + op_size]))

            read_pos += op_size

        elif operation in [BAM_CDEL]:

            alignment[ref_pos : ref_pos + op_size] = [DELETION_NUMBER] * op_size
            ref_pos += op_size

            deletions.append((ref_pos, op_size))

        else:
            raise Exception(f"Operation {operation} of size {op_size} in cigartuples.")

    return alignment, insertions, deletions


def get_second_most_common_freqs(arr):
    def smc_freq(bases):
        try:
            counter = Counter(bases).most_common()
            return [
                freq
                for base, freq in counter
                if base not in (DELETION_NUMBER, DEFAULT_NUMBER)
            ][1]
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


class AlignmentArray:
    def __init__(self, bam_fn, fasta_fn):

        self.reference_name, self.reference_seq = get_sequence(fasta_fn)
        self.reference_seq = transform_sequence_to_letters(self.reference_seq)
        self.reference_size = len(self.reference_seq)
        self.init_array(bam_fn)

    def init_array(self, bam_fn):

        self.reads = []
        self.arr = []

        self.insertions = defaultdict(dict)
        self.deletions = defaultdict(dict)

        with pysam.AlignmentFile(bam_fn) as bam:  # pylint: disable=no-member

            i = 0

            for read in bam.fetch():

                if read.is_supplementary:
                    continue

                if self.reference_size - read.reference_end > DISTANCE_TO_END_THRESHOLD:
                    continue

                (
                    alignment_on_reference,
                    insertions,
                    deletions,
                ) = get_alignment_on_reference(read, self.reference_size)

                for ref_pos, seq in insertions:
                    self.insertions[i][ref_pos] = seq

                for ref_pos, size in deletions:
                    self.deletions[i][ref_pos] = size

                self.arr.append(alignment_on_reference)
                self.reads.append(read)

                i += 1

        self.arr = np.vstack(self.arr)

    def get_read_names(self, reads):

        return [read.query_name for i, read in enumerate(self.reads) if i in reads]

    def get_coords(self, bases_freq_threshold, partition):

        second_most_common_freqs = get_second_most_common_freqs(
            self.arr[partition.reads]
        )

        threshold = sorted(second_most_common_freqs)[-MAX_NUMBER_OF_COORDS]

        return [
            i
            for i, freq in enumerate(second_most_common_freqs)
            if freq >= threshold and freq >= bases_freq_threshold
        ]
