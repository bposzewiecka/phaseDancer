# pylint: disable=consider-using-with

import os
import sys
from collections import defaultdict, namedtuple

from pyfaidx import Fasta

END_TEXT_PATTERN = "END"

PHASEDANCER_DATA_DIR = os.environ["PHASEDANCER_DATA_DIR"]
OUTPUT_DIR = PHASEDANCER_DATA_DIR + "/"

file_txt_pattern = (
    OUTPUT_DIR
    + "data/{0.sample}/{0.contig}/{0.assembler}/{0.cl_type}_{0.cluster}/seq_{number}/seq_{number}.{0.contig}.{0.cl_type}_{0.cluster}.{0.sample}.{0.assembler}.{index_no}.from_index.fasta"  # pylint: disable=line-too-long
)

Contig = namedtuple("Contig", ["contig", "cl_type", "cluster", "sample", "assembler"])

Alignment = namedtuple(
    "Alignment",
    [
        "contig_type",
        "contig_number",
        "contig",
        "contig_length",
        "alignment_length",
        "read_name",
    ],
)


def get_contig(contig_name):
    params = contig_name.split("_")

    return (
        params[0],
        params[1],
        Contig(
            contig=params[2],
            cl_type=params[3],
            cluster=params[4],
            sample=params[5],
            assembler=params[6],
        ),
    )


def get_alignment(paf_line):

    paf_line = paf_line.split()

    contig_type, contig_number, contig = get_contig(paf_line[0])

    contig_length = int(paf_line[1])
    alignment_start = int(paf_line[2])
    alignment_end = int(paf_line[3])
    read_name = paf_line[5]

    return Alignment(
        contig_type=contig_type,
        contig_number=contig_number,
        contig=contig,
        contig_length=contig_length,
        alignment_length=alignment_end - alignment_start,
        read_name=read_name,
    )


def write_reads(fasta, desc, alignments):

    if not alignments:
        return

    alignments_per_read = defaultdict(int)
    contig_length = alignments[0].contig_length

    for alignment in alignments:
        alignments_per_read[alignment.read_name] += alignment.alignment_length

    for read_name, alignment_size in alignments_per_read.items():
        if contig_length * 0.5 < alignment_size:
            desc.write(">" + read_name + "\n")
            desc.write(fasta[read_name][:].seq + "\n")


def multiplex_minimap_output():

    fasta = Fasta(sys.argv[1])
    index_no = sys.argv[2]

    f_log = open(f"log_paf_{index_no}.paf", "w", encoding="UTF-8")

    descriptors = {}
    alignments_by_contig = defaultdict(list)
    contig_numbers = {}

    for line in sys.stdin:

        f_log.write(line)
        f_log.flush()

        if line[0] == "[":
            continue

        alignment = get_alignment(line)
        contig = alignment.contig

        if alignment.contig_type == "seq":

            if contig not in descriptors:
                descriptor_name = file_txt_pattern.format(
                    contig, number=alignment.contig_number, index_no=index_no
                )
                descriptors[contig] = open(descriptor_name, "w", encoding="UTF-8")
                contig_numbers[contig] = alignment.contig_number

            alignments_by_contig[contig].append(alignment)

        if alignment.contig_type == "dummy" and contig in descriptors:

            desc = descriptors[contig]

            write_reads(fasta, desc, alignments_by_contig[contig])

            end_text = END_TEXT_PATTERN
            desc.write(">" + end_text + "\n")
            desc.write(("A" * 10) + "\n")
            desc.close()

            del descriptors[contig]
            del alignments_by_contig[contig]
            del contig_numbers[contig]


multiplex_minimap_output()
