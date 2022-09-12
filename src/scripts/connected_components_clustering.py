import networkx as nx
import numpy as np

from src.scripts.constants import BY_CONNECTED_COMPONENTS
from src.scripts.defaults import (
    CONNECTED_COMPONENTS_HAMMING_SIMILARITY_TRESHOLD,
    CONNECTED_COMPONENTS_MAX_NUMBER_OF_READS,
    COORD_FREQ_TRESHOLD,
    MIN_COORDS_TO_CLUSTER_BY_CONNECTED_COMPONENTS,
)
from src.scripts.partition import Partition


def partition_by_connected_components(partitions):

    new_partitions = []
    clusters = []

    for partition in partitions:

        coords = partition.get_coords(COORD_FREQ_TRESHOLD)

        if len(coords) == 0:
            clusters.append(partition.reads)
            continue

        if len(coords) > MIN_COORDS_TO_CLUSTER_BY_CONNECTED_COMPONENTS:
            new_partitions += cluster_by_connected_components(partition, coords)
        else:
            new_partitions.append(partition)

    return new_partitions, clusters


def cluster_by_connected_components(partition, coords):

    graph = nx.Graph()

    sub_arr = partition.get_sub_arr(coords)

    for i, read1 in enumerate(partition.reads):

        bases1 = sub_arr[i]

        weights = []

        for j, read2 in enumerate(partition.reads):

            bases2 = sub_arr[j]
            weight = np.count_nonzero(bases1 == bases2)
            weights.append((weight, read1, read2))

        most_similar_reads = sorted(weights)[-CONNECTED_COMPONENTS_MAX_NUMBER_OF_READS:]

        for weight, read1, read2 in most_similar_reads:

            if weight / len(coords) > CONNECTED_COMPONENTS_HAMMING_SIMILARITY_TRESHOLD:
                graph.add_edge(read1, read2)

    connected_components = list(nx.connected_components(graph))
    partitions = [
        Partition(
            alignment=partition.alignment,
            reads=list(reads),
            by_type=BY_CONNECTED_COMPONENTS,
        )
        for reads in connected_components
    ]

    return partitions
