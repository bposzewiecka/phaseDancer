import pysam

from src.scripts.constants import (
    BAM_CDEL,
    BAM_CEQUAL,
    BAM_CINS,
    BAM_CMATCH,
    BAM_CSOFT_CLIP,
)
from src.scripts.defaults import MIN_COVERAGE_FOR_EXTENSION, MIN_EXTENSION_SIZE
from src.utils.logging_utils import pLogger
from src.utils.reads_utils import get_reference_length

PHASE_NAME = "READ_TRIMMING"


def split_by_reference_coordinate(read, reference_split):

    if read.reference_start > reference_split:
        raise Exception("Reference start greater then reference split.")

    if reference_split > read.reference_end:
        raise Exception("Reference split greater then reference end.")

    read_pos = 0
    ref_pos = read.reference_start

    for operation, op_size in read.cigartuples:

        if operation in [BAM_CMATCH, BAM_CEQUAL]:
            read_pos += op_size
            ref_pos += op_size
        elif operation in [BAM_CINS, BAM_CSOFT_CLIP]:
            read_pos += op_size
        elif operation in [BAM_CDEL]:
            ref_pos += op_size
        else:
            raise Exception(f"Operation {operation} of size {op_size} in cigartuples.")

        if ref_pos > reference_split:
            # only when op is in [ BAM_CMATCH, BAM_CEQUAL] or [BAM_CDEL]

            read_split = read_pos

            if operation in [BAM_CMATCH, BAM_CEQUAL]:
                read_split -= ref_pos - reference_split + 1

            return (
                read_split,
                read.query_sequence[:read_split],
                read.query_sequence[read_split:],
            )

    raise Exception("Split reference coordinate not found.")


def get_extension_size(read_ends):

    read_ends.sort()

    print("Read ends length", read_ends)

    if len(read_ends) >= MIN_COVERAGE_FOR_EXTENSION:
        return max(
            read_ends[-MIN_COVERAGE_FOR_EXTENSION] // 1000 * 1000 - 1000,
            MIN_EXTENSION_SIZE,
        )

    return MIN_EXTENSION_SIZE


def get_trimmed_reads_fasta(bam_fn, selected_reads, start, end, split_mode):

    reads = {}

    reference_length = get_reference_length(bam_fn)

    if split_mode == "HALF":
        split_coordinate = reference_length // 2
    elif split_mode == "END":
        split_coordinate = reference_length - 100
    else:
        raise Exception

    with pysam.AlignmentFile(bam_fn, "rb") as bam_in:  # pylint: disable=no-member

        read_ends = []

        for read in bam_in.fetch():

            if read.query_name not in selected_reads:
                continue

            # print(
            #    read.reference_start,
            #    read.reference_end,
            #    read.cigartuples[0],
            #    read.cigartuples[-1],
            # )

            if read.is_secondary or read.is_supplementary:
                continue

            if (
                read.reference_start > 100
                and read.cigartuples[0][0] == BAM_CSOFT_CLIP
                and read.cigartuples[0][1] > 50
            ):
                # print("Begining", read.reference_start, read.cigartuples[0][1])

                continue

            if (
                reference_length - read.reference_end > 100
                and read.cigartuples[-1][0] == BAM_CSOFT_CLIP
                and read.cigartuples[-1][1] > 50
            ):
                # print("End", read.reference_end, read.cigartuples[-1][1])
                continue

            if (
                read.reference_start > 100
                and reference_length - read.reference_end > 100
            ):
                # print("In the middle")
                continue

            try:
                _, seq1, seq2 = split_by_reference_coordinate(read, split_coordinate)
            except Exception:  # pylint: disable=broad-except
                # print("Exception")
                continue

            seq = seq1[start - split_coordinate :] + seq2[: end + split_coordinate]

            reads[read.query_name] = seq
            read_ends.append(read.cigartuples[-1][1])

            pLogger.debug(  # pylint: disable=logging-fstring-interpolation
                f"{PHASE_NAME}: Read {read.query_name} was trimmed to {len(seq)}."
            )

    pLogger.info(  # pylint: disable=logging-fstring-interpolation
        f"{PHASE_NAME}: {len(reads)} reads was trimmed."
    )

    return reads, get_extension_size(read_ends)
