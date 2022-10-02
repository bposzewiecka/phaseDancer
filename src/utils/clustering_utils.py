from collections import defaultdict

import numpy as np
import pysam
from scipy.spatial.distance import pdist

from src.scripts.constants import PALETTE


class HammingDistances:
    def __init__(self, consensus_seqs, seqs=None):

        self.consensus_seqs = np.vstack(consensus_seqs)

        if seqs is not None:
            self.seqs = np.vstack(seqs)
            self.data = np.concatenate((self.consensus_seqs, self.seqs), axis=0)
        else:
            self.data = self.consensus_seqs

        self.size = np.shape(self.data)[0]
        distances = pdist(self.data, "hamming")
        self.distances = [
            [self.dist(i, j, distances) for i in range(self.size)]
            for j in range(self.size)
        ]

    def dist(self, i, j, distances):

        if i == j:
            return 0
        if i < j:
            return distances[self.size * i + j - ((i + 2) * (i + 1)) // 2]

        return self.dist(j, i, distances)  # pylint: disable=arguments-out-of-order

    def get_seq_distances(self):

        return [
            self.distances[len(self.consensus_seqs) + i][: len(self.consensus_seqs)]
            for i in range(len(self.seqs))
        ]

    def get_distances(self):
        return self.distances


def select_cluster(clusters, prev_clusters, prev_selected_cluster):

    prev_reads = set(prev_clusters[prev_selected_cluster])

    stats = [len(prev_reads.intersection(reads)) for cluster, reads in clusters.items()]

    print("Clusters to: ", stats)

    selected_cluster = stats.index(max(stats))
    selected_reads = set(clusters[selected_cluster])

    stats = [
        len(selected_reads.intersection(reads))
        for cluster, reads in prev_clusters.items()
    ]

    print("Previous to selected: ", stats)

    return selected_cluster


def read_clusters(cluster_fn):

    clusters = defaultdict(set)

    with open(cluster_fn, encoding="UTF-8") as handler:
        for line in handler:
            cluster, read_name = line.strip().split()
            cluster = int(cluster)
            clusters[cluster].add(read_name)

    return clusters


def save_clusters(clusters_fn, clusters):

    with open(clusters_fn, "w", encoding="UTF-8") as handler:

        for i, cluster in clusters.items():
            for read_no in cluster:
                handler.write(f"{i}\t{read_no}\n")


def inverse_clusters(clusters):

    inv_reads = {}

    for cluster, reads in clusters.items():
        for read in reads:
            inv_reads[read] = cluster

    return inv_reads


def get_colored_bam(in_bam_fn, out_bam_fn, clusters):

    colors = {i: PALETTE[i % len(PALETTE)] for i in range(len(clusters))}

    inverted_clusters = inverse_clusters(clusters)

    with pysam.AlignmentFile(in_bam_fn) as in_bam:  # pylint: disable=no-member

        header = in_bam.header.to_dict()

        with pysam.AlignmentFile(  # pylint: disable=no-member
            out_bam_fn, "wb", header=header
        ) as out_bam:

            read_gen = in_bam.fetch()

            for read in read_gen:

                if read.query_name in inverted_clusters:

                    cluster = inverted_clusters[read.query_name]
                    read.set_tag("YC", colors[cluster])
                    read.set_tag("HP", f"{cluster:03d}")
                    out_bam.write(read)

        pysam.index(out_bam_fn)  # pylint: disable=no-member
