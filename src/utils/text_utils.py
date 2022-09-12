def get_prev_number(number):
    prev_number = int(number) - 1
    return f"{prev_number:03d}"


def get_next_number(number):
    next_number = int(number) + 1
    return f"{next_number:03d}"


def get_range(start=0, end=None):
    return [f"{i:03d}" for i in range(int(start), int(end) + 1)]


def get_next_contig_name(contig_name):
    params = contig_name.split("_")
    params[1] = get_next_number(params[1])
    return "_".join(params)
