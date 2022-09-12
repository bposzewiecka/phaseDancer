import pysam
from Bio import SeqIO


def get_reads_lengths(fasta_fn):

    with open(fasta_fn, encoding="UTF-8") as handler:
        return [len(record) for record in SeqIO.parse(handler, "fasta")]


def get_sequence(fasta_in_fn):

    for record in SeqIO.parse(fasta_in_fn, "fasta"):
        return record.id, str(record.seq)


def get_sequence_reverse_complement(fasta_in_fn):

    for record in SeqIO.parse(fasta_in_fn, "fasta"):
        return record.id, str(record.reverse_complement().seq)


def save_fasta_record(fasta_fn, read_name, seq):

    with open(fasta_fn, "w", encoding="UTF-8") as fasta_out:
        fasta_out.write(">")
        fasta_out.write(read_name)
        fasta_out.write("\n")
        fasta_out.write(seq)
        fasta_out.write("\n")


def save_reads_adding_suffix_to_name(fasta_fn, reads, suffix):

    with open(fasta_fn, "w", encoding="UTF-8") as fasta_out:

        for read_name, seq in reads.items():
            fasta_out.write(f">{read_name}_{suffix}\n")
            fasta_out.write(seq + "\n")


def get_reference_length(bam_fn):

    with pysam.AlignmentFile(bam_fn) as bam:  # pylint: disable=no-member

        if len(bam.lengths) > 1:
            raise Exception(f"More than one reference in file {bam_fn}.")

        return bam.lengths[0]
