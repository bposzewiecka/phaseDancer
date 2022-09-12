import random
from collections import Counter

import numpy as np
from scipy.spatial import distance

from src.scripts.alignments import AlignmentArray, get_most_common_bases
from src.scripts.connected_components_clustering import (
    partition_by_connected_components,
)
from src.scripts.constants import BAM_CSOFT_CLIP, BY_STARTS, SEED
from src.scripts.defaults import (
    MAX_BAM_CSOFT_CLIP_TO_CLASSIFY_AS_RIGHT_READ,
    MIN_DIATANCE_FROM_THE_BEGINNING_TO_CLASSIFY_AS_RIGHT_READ,
    MIN_READS_TO_FORM_PARTITION,
    STARTS_DISTANCE_CLUSTER_THRESHOLD,
)
from src.scripts.partition import Partition
from src.scripts.random_forest_clustering import cluster_by_random_forests


def divide_2_partitions_and_right_reads(alignment_array):

    partitions = []
    right_reads = []
    reads = []
    begin = 0

    def assign_reads(reads, partitions, right_reads, begin, end):

        if len(reads) >= MIN_READS_TO_FORM_PARTITION or begin < 10:
            partition = Partition(
                alignment=alignment_array,
                reads=reads,
                by_type=BY_STARTS,
                start_interval=(begin, end),
            )
            partitions.append(partition)
        else:
            right_reads.extend(reads)

    for i, read in enumerate(alignment_array.reads):

        firstcigar = read.cigartuples[0]

        if firstcigar[0] != BAM_CSOFT_CLIP or (
            firstcigar[1] <= MAX_BAM_CSOFT_CLIP_TO_CLASSIFY_AS_RIGHT_READ
            and read.reference_start
            > MIN_DIATANCE_FROM_THE_BEGINNING_TO_CLASSIFY_AS_RIGHT_READ
        ):
            right_reads.append(i)
            continue

        if (
            reads
            and read.reference_start
            > alignment_array.reads[i - 1].reference_start
            + STARTS_DISTANCE_CLUSTER_THRESHOLD
        ):

            assign_reads(
                reads,
                partitions,
                right_reads,
                begin,
                alignment_array.reads[i - 1].reference_start,
            )

            if read.reference_start > alignment_array.reference_size / 2:
                return partitions, right_reads

            reads = [i]
            begin = read.reference_start

        else:
            reads.append(i)

    assign_reads(
        reads, partitions, right_reads, begin, alignment_array.reads[-1].reference_start
    )

    return partitions, right_reads


def get_counters(arr):

    return np.apply_along_axis(Counter, axis=0, arr=arr)


def get_consensus_seq(clusters, alignment_array):

    return [get_most_common_bases(alignment_array.arr[cluster]) for cluster in clusters]


def get_unique_loci_by_cluster(clusters, alignment_array):

    consensus_seqs = get_consensus_seq(clusters, alignment_array)

    consensus_seqs = np.vstack(consensus_seqs)

    unique_loci = [[]] * len(clusters)

    for locus, counter in enumerate(get_counters(consensus_seqs)):

        if len(counter) != 2:
            continue

        for val, freq in counter.most_common():

            if freq == 1:
                cluster = list(consensus_seqs[:, locus]).index(val)
                unique_loci[cluster].append((locus, val))
                print((locus, val, freq))

    return unique_loci


def classify_right_reads(right_reads, alignment_array, clusters):

    if len(clusters) == 1:
        clusters[0] += right_reads
        return

    unique_loci_by_cluster = get_unique_loci_by_cluster(clusters, alignment_array)

    for i in right_reads:

        read = alignment_array.reads[i]
        distances = []

        for unique_loci in unique_loci_by_cluster:

            distances.append(
                np.average(
                    [
                        alignment_array.arr[i][locus] != val
                        for locus, val in unique_loci
                        if locus > read.reference_start
                    ]
                )
            )

        ranking = sorted(enumerate(distances), key=lambda x: x[1])
        ranking = [rank[0] for rank in ranking]

        if (
            distances[ranking[0]] < 0.25
            and distances[ranking[1]] - distances[ranking[0]] > 0.2
        ):
            clusters[ranking[0]].append(i)


def sort_by_similarity_to_ref(clusters, alignment_array):

    consensus_seqs = get_consensus_seq(clusters, alignment_array)

    hamming_distances = [
        distance.hamming(consensus_seq, alignment_array.reference_seq)
        for consensus_seq in consensus_seqs
    ]

    sorted_clusters = sorted(
        range(0, len(clusters)), key=lambda x: hamming_distances[x]
    )

    return sorted_clusters


def cluster_reads(bam_fn, fasta_fn):

    alignment_array = AlignmentArray(bam_fn, fasta_fn)

    partitions, right_reads = divide_2_partitions_and_right_reads(alignment_array)
    partitions, clusters = partition_by_connected_components(partitions)

    random.seed(SEED)

    for partition in partitions:
        clusters += cluster_by_random_forests(partition)

    classify_right_reads(right_reads, alignment_array, clusters)

    return clusters
