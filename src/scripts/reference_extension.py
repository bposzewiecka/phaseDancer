# pylint: disable=c-extension-no-member

import edlib

from src.utils.reads_utils import (
    get_sequence,
    get_sequence_reverse_complement,
    save_fasta_record,
)
from src.utils.text_utils import get_next_contig_name
from src.utils.yaml_utils import save_yaml


def align_assembled_contig(contig_fn, assembled_contig_fn):

    print(contig_fn)
    print(assembled_contig_fn)

    contig_name, contig = get_sequence(contig_fn)
    _, assembled_contig = get_sequence(assembled_contig_fn)
    _, assembled_contig_reverse_complement = get_sequence_reverse_complement(
        assembled_contig_fn
    )

    alignment = edlib.align(
        contig, assembled_contig, mode="HW", task="path"
    )  # pylint: disable=c-extension-no-member
    alignment_rev_complement = edlib.align(
        contig, assembled_contig_reverse_complement, mode="HW", task="path"
    )  # pylint: disable=c-extension-no-member

    if alignment["editDistance"] < alignment_rev_complement["editDistance"]:
        if alignment["editDistance"] > 600:
            raise Exception
        return alignment, assembled_contig, contig_name, contig

    if alignment_rev_complement["editDistance"] > 600:
        raise Exception
    return (
        alignment_rev_complement,
        assembled_contig_reverse_complement,
        contig_name,
        contig,
    )


def extend_contig(
    contig_fn,
    assembled_contig_fn,
    extended_contig_fn,
    extension_info_fn,
    contig_extension_size,
):

    result, assembled_contig, contig_name, contig = align_assembled_contig(
        contig_fn, assembled_contig_fn
    )

    contig_size = len(contig)

    _, end_coord_contig = result["locations"][0]
    extended_left = assembled_contig[end_coord_contig + 1 :][:contig_extension_size]
    extended_left_size = len(extended_left)

    print(f"Contig extended by {len(extended_left)}", result)

    extended_contig = contig[extended_left_size - contig_size :] + extended_left
    extended_contig_name = get_next_contig_name(contig_name)

    save_fasta_record(extended_contig_fn, extended_contig_name, extended_contig)

    result["contig-extension-size"] = extended_left_size
    result["locations"] = [
        {"start": start, "end": end} for start, end in result["locations"]
    ]

    save_yaml(extension_info_fn, result)
