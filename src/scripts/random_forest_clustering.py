import random
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import distance

from src.scripts.alignments import get_most_common_bases, get_second_most_common_freqs
from src.scripts.constants import DEFAULT_NUMBER, DELETION_NUMBER
from src.scripts.defaults import (
    COORD_FREQ_TRESHOLD,
    MIN_READ_ACCURACY,
    NUMBER_OF_SIMULATIONS,
)


def simulate_clustering(arr, reads):
    def cluster(rows, coords):

        sub_arr = arr[rows]
        sub_arr = sub_arr[:, coords]

        smc_freqs = get_second_most_common_freqs(sub_arr)

        indices = [i for i, freq in enumerate(smc_freqs) if freq >= freq_threshold]
        coords = [
            coord for coord, freq in zip(coords, smc_freqs) if freq >= freq_threshold
        ]

        if indices:

            index = random.randrange(len(indices))
            cluster_by = indices[index]
            coords = coords[:index] + coords[index + 1 :]

            smc = Counter(
                base
                for base in sub_arr[:, cluster_by]
                if base not in (DELETION_NUMBER, DEFAULT_NUMBER)
            ).most_common()[1][0]

            rows1 = [
                row for row, base in zip(rows, sub_arr[:, cluster_by]) if base == smc
            ]
            rows2 = [
                row
                for row, base in zip(rows, sub_arr[:, cluster_by])
                if base not in (smc, DELETION_NUMBER, DEFAULT_NUMBER)
            ]

            partitions.append((rows1, coords.copy()))
            partitions.append((rows2, coords.copy()))

        else:

            clusters.append(rows)

    def realign_results(clusters):

        consensus_seqs = [get_most_common_bases(arr[cluster]) for cluster in clusters]

        consensus_seqs = np.vstack(consensus_seqs)

        row_distances = [
            [distance.hamming(seq, row) for seq in consensus_seqs] for row in arr
        ]

        seqs_number = [0] * len(consensus_seqs)

        for distances in row_distances:
            seqs_number[np.argmin(distances)] += 1

        consensus_seqs = [
            seq
            for seq, number in zip(consensus_seqs, seqs_number)
            if number >= freq_threshold
        ]

        row_distances = [
            [distance.hamming(seq, row) for seq in consensus_seqs] for row in arr
        ]

        seqs_number = [0] * len(consensus_seqs)

        clusters = defaultdict(list)
        dist_sum = 0

        for i, distances in enumerate(row_distances):
            cluster = np.argmin(distances)
            dist_sum += distances[cluster] * ncols
            clusters[cluster].append(reads[i])
            seqs_number[cluster] += 1

        return dist_sum, clusters.values()

    nrows, ncols = np.shape(arr)

    clusters = []
    partitions = [(list(range(nrows)), list(range(ncols)))]

    freq_threshold = int(COORD_FREQ_TRESHOLD * MIN_READ_ACCURACY)

    while partitions:

        rows, coords = partitions.pop()

        if coords:
            cluster(rows, coords)
        else:
            clusters.append(rows)

    return realign_results(clusters)


def cluster_by_random_forests(partition):

    coords = partition.get_coords(COORD_FREQ_TRESHOLD)

    if len(coords) == 0:
        return [partition.reads]

    sub_arr = partition.get_sub_arr(coords)

    min_rand_dist_sums = defaultdict(lambda: 1000000000)
    min_rand_clusters = {}

    random_clusters = [
        simulate_clustering(sub_arr, partition.reads)
        for i in range(NUMBER_OF_SIMULATIONS)
    ]

    for dist_sum, clusters in random_clusters:

        cluster_len = len(clusters)

        if min_rand_dist_sums[cluster_len] > dist_sum:
            min_rand_dist_sums[cluster_len] = dist_sum
            min_rand_clusters[cluster_len] = clusters

    min_rand_dist_sums_sorted = sorted(min_rand_dist_sums.items(), key=lambda x: x[1])
    print(len(coords), len(partition.reads), min_rand_dist_sums_sorted)

    return min_rand_clusters[min_rand_dist_sums_sorted[0][0]]
