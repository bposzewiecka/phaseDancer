import random
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import distance

from src.scripts.alignments import AlignmentArray, get_consensus_seq
from src.scripts.connected_components_clustering import (
    partition_by_connected_components,
)
from src.scripts.constants import BY_STARTS, SEED
from src.scripts.defaults import (
    MAX_DISTANCE_TO_JOIN_PARTITIONS,
    MIN_READS_TO_FORM_PARTITION,
)
from src.scripts.partition import Partition
from src.scripts.random_forest_clustering import cluster_by_random_forests
from src.utils.clustering_utils import select_cluster


def divide_2_partitions_and_right_reads(alignment_array):

    reference_starts = defaultdict(int)
    right_reads = []
    partitions_starts = [list(range())]

    for read in alignment_array.reads:

        reference_starts[read.reference_start] += 1

    last_start = partitions_starts[-1][-1]

    for start, freq in reference_starts.items():

        if freq < MIN_READS_TO_FORM_PARTITION or start < last_start:
            continue

        if start - last_start <= MAX_DISTANCE_TO_JOIN_PARTITIONS:
            partitions_starts[-1].append(start)
        else:
            partitions_starts.append([start])
        last_start = start

    inv_starts = {
        start: i for i, starts in enumerate(partitions_starts) for start in starts
    }

    partitions = [[] for _ in partitions_starts]

    print(f"Partitions starts: {partitions_starts}")

    for i, read in enumerate(alignment_array.reads):
        if read.reference_start in inv_starts:
            partitions[inv_starts[read.reference_start]].append(i)
        else:
            right_reads.append(i)

    partitions = [
        Partition(alignment=alignment_array, reads=reads, by_type=BY_STARTS)
        for reads in partitions
    ]

    return partitions, right_reads


def get_counters(arr):

    return np.apply_along_axis(Counter, axis=0, arr=arr)


def get_unique_loci_by_cluster(clusters, alignment_array):

    consensus_seqs = get_consensus_seq(clusters, alignment_array)

    consensus_seqs = np.vstack(consensus_seqs)

    unique_loci = [[] for _ in range(len(clusters))]

    for locus, counter in enumerate(get_counters(consensus_seqs)):

        if len(counter) != 2:
            continue

        for val, freq in counter.most_common():

            if freq == 1:

                cluster = list(consensus_seqs[:, locus]).index(val)

                unique_loci[cluster].append((locus, val))

    return unique_loci


def classify_right_reads(right_reads, alignment_array, clusters, selected_cluster):

    consensus_seqs = get_consensus_seq(clusters, alignment_array)
    selected_cluster_consensus_seq = consensus_seqs[selected_cluster]

    diffs = [
        [
            i
            for i, base in enumerate(consensus_seq)
            if base != selected_cluster_consensus_seq[i]
        ]
        for consensus_seq in consensus_seqs
    ]

    for i in right_reads:

        read = alignment_array.reads[i]
        consensus_seq_scores = [
            np.average(
                [
                    alignment_array.arr[i][coord] == consensus_seqs[j][coord]
                    for coord in diff_coords
                    if coord > read.reference_start
                ]
            )
            for j, diff_coords in enumerate(diffs)
            if j != selected_cluster
        ]
        selected_cluster_consensus_seq_scores = [
            np.average(
                [
                    alignment_array.arr[i][coord]
                    == selected_cluster_consensus_seq[coord]
                    for coord in diff_coords
                    if coord > read.reference_start
                ]
            )
            for j, diff_coords in enumerate(diffs)
            if j != selected_cluster
        ]

        if all(
            c < sc and not np.isnan(c)
            for c, sc in zip(
                consensus_seq_scores, selected_cluster_consensus_seq_scores
            )
        ):
            clusters[selected_cluster].append(i)


def sort_by_similarity_to_ref(clusters, alignment_array):

    consensus_seqs = get_consensus_seq(clusters, alignment_array)

    hamming_distances = [
        distance.hamming(consensus_seq, alignment_array.reference_seq)
        for consensus_seq in consensus_seqs
    ]

    sort_order = sorted(range(0, len(clusters)), key=lambda x: hamming_distances[x])

    return [clusters[i] for i in sort_order]


def cluster_reads(bam_fn, fasta_fn, prev_clusters, prev_selected_cluster):

    alignment_array = AlignmentArray(bam_fn, fasta_fn)

    partitions, right_reads = divide_2_partitions_and_right_reads(alignment_array)

    partitions, clusters = partition_by_connected_components(partitions)

    random.seed(SEED)

    for partition in partitions:
        clusters += cluster_by_random_forests(partition)

    clusters = sort_by_similarity_to_ref(clusters, alignment_array)

    if prev_selected_cluster is None:
        # clusters are sorted by similarity to reference sequence
        selected_cluster = 0
    else:
        selected_cluster = select_cluster(
            {
                i: [alignment_array.reads[read].query_name for read in cluster]
                for i, cluster in enumerate(clusters)
            },
            prev_clusters,
            prev_selected_cluster,
        )

    classify_right_reads(right_reads, alignment_array, clusters, selected_cluster)

    clusters_with_read_names = {
        i: [alignment_array.reads[read].query_name for read in cluster]
        for i, cluster in enumerate(clusters)
    }

    print(f"Selected cluster: {selected_cluster}")

    return clusters_with_read_names, selected_cluster
