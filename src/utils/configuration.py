# pylint: disable=line-too-long

import os.path
import re
from collections import namedtuple

from src.utils.reads_utils import get_reads_lengths
from src.utils.yaml_utils import load_yaml

Property = namedtuple(
    "Property",
    [
        "name",
        "level",
        "prop_types",
        "sample_required",
        "contig_required",
        "default_value",
        "validation",
    ],
)

TECHNOLOGIES = ("ont", "pb", "hifi")

TECHNOLOGIES_CONFIG_FN = "config_technologies.yaml"

CONFIG_FN = "config.yaml"


def get_number_of_indices(wildcards, config):

    number_of_indices = 1

    if "number-of-indices" in config["samples"][wildcards.sample]:
        number_of_indices = int(
            config["samples"][wildcards.sample]["number-of-indices"]
        )

    return number_of_indices


def technology_validation(
    key, value, message_prefix
):  # pylint: disable=unused-argument

    if value not in TECHNOLOGIES:
        technologies = ", ".join(TECHNOLOGIES)
        return [
            f'{message_prefix} Technology "{value}" is not supported. Possible values are: {technologies}.'
        ]

    return []


def mappings_validation(key, value, message_prefix):  # pylint: disable=unused-argument

    errors = []

    for i, ref in enumerate(value):
        if isinstance(ref, str):
            errors.append(
                f'{message_prefix} {i}-th element of "mapping" should (0-base) be a string.'
            )

    for ref in value:
        ref_fn = f"data/refs/{ref}.fasta"
        if not os.path.isfile(ref_fn):
            errors.append(
                f'{message_prefix} File "{ref_fn}" required for mapping not found.'
            )

    return errors


def non_negative_integer(key, value, message_prefix):

    errors = []

    if value < 0:
        errors.append(
            f'{message_prefix} Value of "{key}" should be a non-negative integer.'
        )

    return errors


def positive_integer(key, value, message_prefix):

    errors = []

    if value <= 0:
        errors.append(
            f'{message_prefix} Value of "{key}" should be a positive integer.'
        )

    return errors


properties = [
    Property(
        name="technology",
        level="sample",
        prop_types=[str],
        sample_required=True,
        contig_required=False,
        default_value=None,
        validation=technology_validation,
    ),
    Property(
        name="contigs",
        level="sample",
        prop_types=[dict, list],
        sample_required=True,
        contig_required=False,
        default_value=None,
        validation=None,
    ),
    Property(
        name="contig-size",
        level="both",
        prop_types=[int],
        sample_required=False,
        contig_required=True,
        default_value=None,
        validation=positive_integer,
    ),
    Property(
        name="contig-extension-size",
        level="both",
        prop_types=[int],
        sample_required=False,
        contig_required=True,
        default_value=None,
        validation=positive_integer,
    ),
    Property(
        name="browser",
        level="both",
        prop_types=[bool],
        sample_required=False,
        contig_required=False,
        default_value=True,
        validation=None,
    ),
    Property(
        name="iterations",
        level="both",
        prop_types=[int],
        sample_required=False,
        contig_required=True,
        default_value=0,
        validation=non_negative_integer,
    ),
    Property(
        name="mappings",
        level="both",
        prop_types=[list],
        sample_required=False,
        contig_required=False,
        default_value=[],
        validation=mappings_validation,
    ),
]

properties_dict = {prop.name: prop for prop in properties}

sample_keys = [prop.name for prop in properties if prop.level in ("sample", "both")]
contig_keys = [prop.name for prop in properties if prop.level in ("contig", "both")]

sample_required_keys = [
    prop.name
    for prop in properties
    if prop.sample_required and prop.default_value is None
]
contig_required_keys = [
    prop.name
    for prop in properties
    if prop.contig_required and prop.default_value is None
]

TECHNOLOGIES = ("ont", "pb", "hifi")

property_types = {prop.name: prop.prop_types for prop in properties}
default_values = {prop.name: prop.default_value for prop in properties}

type_info = {
    str: "string",
    int: "number",
    bool: "boolean",
    dict: "dictionary",
    list: "list",
}


def load_config(config_fn):
    try:
        return load_yaml(config_fn)
    except BaseException as exp:
        raise Exception(
            f'Problem with loading the configuration file "{config_fn}".'
        ) from exp


def get_value(key, contig, sample):

    if contig is not None and key in contig:
        return contig[key]

    return sample.get(key, default_values[key])


def validate_name(name, message_prefix, level):

    if re.match("^[a-zA-Z][a-zA-Z0-9-]*$", name, re.ASCII) is None:
        return [
            f'{message_prefix} Invalid {level} name "{name}". {level.capitalize()} name should begin with a letter (a-zA-Z) and contain only letters (a-zA-Z), digits (0-9) and hyphens.'
        ]

    return []


def validate_key_names(obj, obj_keys, message_prefix, level):
    errors = []

    for key in obj:
        if key not in obj_keys:
            errors.append(
                f'{message_prefix} Key "{key}" is not allowed in the {level} configuration.'
            )

    return errors


def validate_value_types(obj, message_prefix):
    errors = []

    for key, value in obj.items():
        if type(value) not in property_types[key]:
            possible_types = " or ".join(
                [type_info[property_type] for property_type in property_types[key]]
            )
            errors.append(
                f'{message_prefix} Value of "{key}" should be a {possible_types}.'
            )

    return errors


def validate_values(obj, message_prefix):
    errors = []

    for key, value in obj.items():

        if properties_dict[key].validation:

            errors += properties_dict[key].validation(key, value, message_prefix)

    return errors


def validate_props(obj_name, obj, obj_keys, message_prefix, level):

    if isinstance(obj, dict):
        return [f"{message_prefix} Configuration should be a dictionary."]

    errors = validate_name(obj_name, message_prefix, level)

    if errors:
        return errors

    errors = validate_key_names(obj, obj_keys, message_prefix, level)

    if errors:
        return errors

    errors = validate_value_types(obj, message_prefix)

    if errors:
        return errors

    errors = validate_values(obj, message_prefix)

    return errors


def validate_sample_files(sample_name, message_prefix):

    file_names = [
        f"data/{sample_name}/{sample_name}." + ext for ext in ["fasta", "fai", "mmi"]
    ]

    errors = []

    for file_name in file_names:
        if not os.path.isfile(file_name):
            errors.append(
                f'{message_prefix} Required  file "{file_name}" does not exist.'
            )

    return errors


def validate_required(required_keys, objs, message_prefix):

    errors = []

    obj_keys = [key for obj in objs for key in obj.keys()]

    for required_key in required_keys:
        if required_key not in obj_keys:
            errors.append(f'{message_prefix} Required key "{required_key}" not found.')

    return errors


def validate_fasta(contig_name, contig, sample_name, sample, message_prefix):

    errors = []

    fasta_fn = f"data/{sample_name}/{contig_name}.fasta"

    if not os.path.isfile(fasta_fn):
        return [
            f'{message_prefix} Fasta file "{fasta_fn}" with start contig does not exist.'
        ]

    contig_length = get_value("contig-size", contig, sample)

    fasta_fn = f"data/{sample_name}/{contig_name}.fasta"

    try:
        read_lengths = get_reads_lengths(fasta_fn)

        if len(read_lengths) == 1 and read_lengths[0] != contig_length:
            errors.append(
                f'{message_prefix} Size of contig from "{fasta_fn}" file ({read_lengths[0]}) is not equal to the  size specified in the configuration file ({contig_length}).'
            )
        if len(read_lengths) != 1:
            errors.append(
                f'{message_prefix} Fasta file "{fasta_fn}" should contain one read. Now it contains {len(read_lengths)} reads.'
            )

    except BaseException:  # pylint: disable=broad-except
        errors.append(f'{message_prefix} Fasta file "{fasta_fn}" is malformed.')

    return errors


def validate_contig_level(contig, sample, message_prefix):

    errors = validate_required(contig_required_keys, [contig, sample], message_prefix)

    if errors:
        return errors

    if (
        get_value("iterations-left", contig, sample)
        + get_value("iterations-right", contig, sample)
        == 0
    ):
        errors.append(
            f'{message_prefix} "iterations-right" or "iterations-left" should not be 0.'
        )

    return errors


def validate_sample(sample_name, sample, check_files):

    errors = []

    message_prefix = f'Sample "{sample_name}":'

    errors = validate_props(sample_name, sample, sample_keys, message_prefix, "sample")

    if errors:
        return errors

    errors = validate_required(sample_required_keys, [sample], message_prefix)

    if errors:
        return errors

    if check_files:
        errors = validate_sample_files(sample_name, message_prefix)
        if errors:
            return errors

    contigs = sample["contigs"]

    if isinstance(contigs, dict):

        for contig_name, contig in contigs.items():
            errors += validate_contig(
                contig_name,
                contig,
                sample_name,
                sample,
            )

    else:
        errors = validate_contig_level(sample, {}, message_prefix)
        contigs = {contig_name: {} for contig_name in contigs}

    if errors:
        return errors

    if check_files:
        for contig_name, contig in contigs.items():
            errors += validate_fasta(
                contig_name, contig, sample_name, sample, message_prefix
            )

    return errors


def validate_contig(contig_name, contig, sample_name, sample):

    errors = []

    message_prefix = f'Sample "{sample_name}", contig "{contig_name}":'

    errors = validate_props(contig_name, contig, contig_keys, message_prefix, "contig")

    if errors:
        return errors

    errors = validate_contig_level(contig, sample, message_prefix)

    return errors


def validate_config(config_fn, check_files=True):

    config = load_config(config_fn)

    if isinstance(config, dict):
        return [f"Config file {config_fn} should be a dictionary."]

    if "samples" not in config or len(config) != 1:
        return ['Configuration should have one root key "samples".']

    if isinstance(config["samples"], dict):
        return ['Value of root key "samples" should be a dictionary.']

    errors = []

    for sample_name, sample in config["samples"].items():
        errors += validate_sample(sample_name, sample, check_files)

    return errors


def snakemake_validate_config(config_fn, check_files):

    errors = validate_config(config_fn, check_files=check_files)

    if errors:
        print("Errors:")

        for error in errors:
            print(error)

    return errors


def get_contigs_to_extend(config):

    if "contigs" not in config:
        raise Exception(
            "Specify config parameter (see Documentation - Step 8: Starting the main algorithm)"
        )

    config_text = config["contigs"]

    def convert_to_dict(objs):
        if isinstance(objs, dict):
            return objs
        return {obj: {} for obj in objs}

    def get_params(contig, contig_name, sample, sample_name):
        return {
            "sample": sample_name,
            "contig": contig_name,
            "iterations": get_value("iterations", contig, sample),
        }

    if config_text == "all":
        return [
            get_params(contig, contig_name, sample, sample_name)
            for sample_name, sample in config["samples"].items()
            for contig_name, contig in convert_to_dict(sample["contigs"]).items()
        ]

    sample_names = []

    params_list = []

    for sample_config in config_text.split(";"):
        sample_config = sample_config.split(":")

        sample_name = sample_config[0]

        if sample_name not in config["samples"]:
            raise Exception(
                f'Sample "{sample_name}" not found in the configuration file "{CONFIG_FN}".'
            )

        if sample_name in sample_names:
            raise Exception(
                f'Sample "{sample_name}" listed more than once in the "contigs" parameter.'
            )

        sample_names.append(sample_name)

        sample = config["samples"][sample_name]

        if len(sample_config) == 1:

            params_list += [
                get_params(contig, contig_name, sample, sample_name)
                for contig_name, contig in convert_to_dict(sample["contigs"]).items()
            ]

        elif len(sample_config) == 2:

            contigs_text = sample_config[1]

            contigs = convert_to_dict(sample["contigs"])

            contig_names = []

            for contig_name in contigs_text.split(","):

                if contig_name not in contigs.keys():
                    raise Exception(
                        f'Sample "{sample_name}": Contig "{contig_name}" not found in the configuration file "{CONFIG_FN}".'
                    )

                if contig_name in contig_names:
                    raise Exception(
                        f'Sample "{sample_name}": Contig "{contig_name}" listed more than once in the "contigs" parameter"'
                    )

                contig_names.append(contig_name)

                contig = contigs[contig_name]

                params_list.append(get_params(contig, contig_name, sample, sample_name))

        else:
            raise Exception(f'Sample "{sample_name}": Contig list is not valid')

        return params_list
