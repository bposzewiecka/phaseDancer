import random
from collections import Counter, defaultdict
from multiprocessing import get_context

import numpy as np
from scipy.spatial import distance

from src.scripts.alignments import get_most_common_bases, get_second_most_common_freqs
from src.scripts.constants import get_excluded_values
from src.scripts.defaults import (  # COORD_FREQ_TRESHOLD,; MIN_READ_ACCURACY,
    NUMBER_OF_SIMULATIONS,
)
from src.utils.clustering_utils import HammingDistances


def simulate_clustering(arr, reads, technology, reference_length):

    excluded_values = get_excluded_values(technology)

    def cluster(rows, coords):

        sub_arr = arr[rows]
        sub_arr = sub_arr[:, coords]

        smc_freqs = get_second_most_common_freqs(sub_arr, technology)

        stats = [
            (abs(len(rows) // 2 - freq), i, coord < reference_length)
            for i, (freq, coord) in enumerate(zip(smc_freqs, coords))
            if freq >= freq_threshold
        ]

        indices = [i for i, freq in enumerate(smc_freqs) if freq >= freq_threshold]

        coords = [
            coord for coord, freq in zip(coords, smc_freqs) if freq >= freq_threshold
        ]

        stats.sort()

        most_freq = [i for _, i, simple in stats if simple][:10]

        if not most_freq:
            most_freq = [i for _, i, _ in stats][:10]

        if indices:

            index = random.randrange(len(most_freq))

            cluster_by = most_freq[index]
            coords = coords[:index] + coords[index + 1 :]

            smc = Counter(
                base for base in sub_arr[:, cluster_by] if base not in excluded_values
            ).most_common()[1][0]

            rows1 = [
                row for row, base in zip(rows, sub_arr[:, cluster_by]) if base == smc
            ]
            rows2 = [
                row
                for row, base in zip(rows, sub_arr[:, cluster_by])
                if base not in (smc,) + excluded_values
            ]

            partitions.append((rows1, coords.copy()))
            partitions.append((rows2, coords.copy()))

        else:

            clusters.append(rows)

    def realign_results(clusters):

        consensus_seqs = [get_most_common_bases(arr[cluster]) for cluster in clusters]

        row_distances = [
            [distance.hamming(seq, row) for seq in consensus_seqs] for row in arr
        ]

        seqs_number = [0] * len(consensus_seqs)

        for distances in row_distances:
            seqs_number[np.argmin(distances)] += 1

        # consensus_seqs = [
        #    seq
        #    for seq, number in zip(consensus_seqs, seqs_number)
        #    if number >= freq_threshold
        # ]

        # row_distances = [
        #    [distance.hamming(seq, row) for seq in consensus_seqs] for row in arr
        # ]

        seqs_number = [0] * len(consensus_seqs)

        clusters = defaultdict(list)
        dist_sum = 0

        for i, distances in enumerate(row_distances):
            cluster = np.argmin(distances)
            dist_sum += distances[cluster] * ncols
            clusters[cluster].append(reads[i])
            seqs_number[cluster] += 1

        return dist_sum, list(clusters.values()), consensus_seqs

    nrows, ncols = np.shape(arr)

    clusters = []
    partitions = [(list(range(nrows)), list(range(ncols)))]

    freq_threshold = 8  # int(COORD_FREQ_TRESHOLD * MIN_READ_ACCURACY)

    while partitions:

        rows, coords = partitions.pop()

        if coords:
            cluster(rows, coords)
        else:
            clusters.append(rows)

    return realign_results(clusters)


def cluster_by_random_forests(partition, technology):  # pylint: disable=too-many-locals

    coords = partition.get_coords(8)  # COORD_FREQ_TRESHOLD)

    if len(coords) == 0:
        return [partition.reads]

    sub_arr = partition.get_sub_arr(coords)

    min_rand_dist_sums = defaultdict(lambda: 1000000000)
    min_rand_clusters = {}
    min_rand_consensus_seqs = {}

    reference_size = partition.get_reference_size()

    if len(coords) > 1000:

        with get_context("spawn").Pool(processes=NUMBER_OF_SIMULATIONS) as pool:
            random_clusters = pool.starmap(
                simulate_clustering,
                [
                    (sub_arr, partition.reads, technology, reference_size)
                    for _ in range(NUMBER_OF_SIMULATIONS)
                ],
            )

    else:
        random_clusters = [
            simulate_clustering(sub_arr, partition.reads, technology, reference_size)
            for _ in range(NUMBER_OF_SIMULATIONS)
        ]

    for dist_sum, clusters, consensus_seqs in random_clusters:

        cluster_len = len(clusters)

        if min_rand_dist_sums[cluster_len] > dist_sum:
            min_rand_dist_sums[cluster_len] = dist_sum
            min_rand_clusters[cluster_len] = clusters
            min_rand_consensus_seqs[cluster_len] = consensus_seqs

    min_rand_dist_sums_sorted = sorted(min_rand_dist_sums.items(), key=lambda x: x[1])

    # b = np.asarray(min_rand_consensus_seqs[min_rand_dist_sums_sorted[0][0]])
    # a = b.transpose()

    # row_distances = [[distance.hamming(seq, row) for seq in b] for row in sub_arr]

    # for name, dist in zip(partition.get_read_names(), row_distances):
    #    print(name, min(dist))

    print(
        "Random forest stats: ",
        len(coords),
        len(partition.reads),
        min_rand_dist_sums_sorted,
        flush=True,
    )

    cons = min_rand_consensus_seqs[min_rand_dist_sums_sorted[0][0]]

    print(cons)

    for dist in HammingDistances(cons).get_distances():
        print([len(cons[0]) * d for d in dist])

    clusters = min_rand_clusters[min_rand_dist_sums_sorted[0][0]]

    print([len(cluster) for cluster in clusters])

    return clusters
