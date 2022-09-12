# pylint: disable=line-too-long

from src.utils.reads_utils import get_sequence, save_fasta_record
from src.utils.yaml_utils import load_yaml

CONTIG_PATTERN_FN = "data/{sample}/{contig}/{assembler}/{cl_type}_{cluster}/seq_{{number:03d}}/seq_{{number:03d}}.{contig}.{cl_type}_{cluster}.{sample}.{assembler}.uncompressed.fasta"
EXTENSION_INFO_PATTERN_FN = "data/{sample}/{contig}/{assembler}/{cl_type}_{cluster}/seq_{{number:03d}}/seq_{{number:03d}}.{contig}.{cl_type}_{cluster}.{sample}.{assembler}.extension_info.yaml"
MERGED_CONTIGS_SEQ_NAME_PATTERN = "merged_contig_seq-start_{start_number}_seq-end_{end_number}_{contig}_fc_all_{sample}_{assembler}"


def merge_contigs(fasta_fn, wildcards, output_dir):

    contig_pattern_fn = output_dir + CONTIG_PATTERN_FN.format(**wildcards)
    extension_info_pattern_fn = output_dir + EXTENSION_INFO_PATTERN_FN.format(
        **wildcards
    )

    merged_contig_seq = ""

    start = int(wildcards["start_number"])
    end = int(wildcards["end_number"])

    for i in range(start, end + 1):

        contig_fn = contig_pattern_fn.format(number=i)
        _, contig_seq = get_sequence(contig_fn)

        if i == start:
            merged_contig_seq = contig_seq
        else:
            extension_info_fn = extension_info_pattern_fn.format(number=i)
            extension_info = load_yaml(extension_info_fn)
            contig_extension_size = extension_info["contig-extension-size"]
            merged_contig_seq += contig_seq[-contig_extension_size:]

    merged_contigs_seq_name = MERGED_CONTIGS_SEQ_NAME_PATTERN.format(**wildcards)
    save_fasta_record(fasta_fn, merged_contigs_seq_name, merged_contig_seq)
