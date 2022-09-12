from src.scripts.clustering import cluster_reads
from src.scripts.read_selection import get_trimmed_reads_fasta
from src.utils.clustering_utils import (
    get_colored_bam,
    read_clusters,
    save_clusters,
    select_cluster,
)
from src.utils.reads_utils import save_reads_adding_suffix_to_name
from src.utils.yaml_utils import load_yaml, save_yaml


def save_clusters_as_bams(params, clusters):

    compressed_bam_fn = params["compressed_bam"]
    compressed_colored_bam_fn = params["compressed_colored_bam"]

    uncompressed_bam_fn = params["uncompressed_bam"]
    uncompressed_colored_bam_fn = params["uncompressed_colored_bam"]

    get_colored_bam(uncompressed_bam_fn, uncompressed_colored_bam_fn, clusters)
    get_colored_bam(compressed_bam_fn, compressed_colored_bam_fn, clusters)


def cluster_and_select_cluster(params, number):

    prev_selected_cluster_fn = params["prev_selected_cluster"]
    prev_clusters_fn = params["prev_clusters"]

    selected_cluster_fn = params["selected_cluster"]
    clusters_fn = params["clusters"]

    compressed_bam_fn = params["compressed_bam"]
    compressed_fasta_fn = params["compressed_fasta"]

    prev_selected_cluster = None
    prev_clusters = None

    if number > 0:
        prev_selected_cluster = load_yaml(prev_selected_cluster_fn)["selected"]
        prev_clusters = read_clusters(prev_clusters_fn)

    clusters = cluster_reads(compressed_bam_fn, compressed_fasta_fn)

    if prev_selected_cluster is None:
        # clusters are sorted by similarity to reference sequence
        selected_cluster = 0
    else:
        selected_cluster = select_cluster(
            clusters, prev_clusters, prev_selected_cluster
        )

    save_yaml(selected_cluster_fn, {"selected": selected_cluster})
    save_clusters(clusters_fn, clusters)

    return clusters, selected_cluster


def extract_reads(params, selected_reads):

    uncompressed_bam_fn = params["uncompressed_bam"]

    extension_size_fn = params["extension_size"]
    reads_fn = params["reads"]

    trimmed_reads = set()

    reads, extension_size = get_trimmed_reads_fasta(
        uncompressed_bam_fn,
        selected_reads,
        -2500,
        10000,
        split_mode="HALF",
    )

    trimmed_reads.update(reads)

    save_yaml(extension_size_fn, {"extension_size": extension_size})

    save_reads_adding_suffix_to_name(reads_fn, trimmed_reads, "trimmed")


def cluster_and_extract_reads(number, params):

    print(f"********************** Step {number} **********************")

    clusters, selected_cluster = cluster_and_select_cluster(params, number)
    save_clusters_as_bams(params, clusters)
    extract_reads(params, clusters[selected_cluster])
