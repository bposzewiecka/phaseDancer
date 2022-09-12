from collections import defaultdict

import pysam

from src.scripts.constants import PALETTE


def select_cluster(clusters, prev_clusters, prev_selected_cluster):

    prev_reads = set(prev_clusters[prev_selected_cluster])

    stats = [len(prev_reads.intersection(cluster)) for cluster in clusters]

    return stats.index(max(stats))


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
