from collections import defaultdict

BY_ROOT = "BY_ROOT"
BY_STARTS = "BY_STARTS"
BY_CONNECTED_COMPONENTS = "BY_CONNECTED_COMPONENTS"
BY_RANDOM_CLUSTERING = "BY_RANDOM_CLUSTERING"
SEED = 666


BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9


LETTERS_TO_NUMBERS = defaultdict(lambda: DEFAULT_NUMBER)
LETTERS_TO_NUMBERS["A"] = 1
LETTERS_TO_NUMBERS["C"] = 2
LETTERS_TO_NUMBERS["G"] = 3
LETTERS_TO_NUMBERS["T"] = 4

DELETION_NUMBER = 5
DEFAULT_NUMBER = 9

NUMBERS_TO_LETTERS = {value: key for key, value in LETTERS_TO_NUMBERS.items()}

PALETTE_HEX = (
    "d61b60",
    "e2188e",
    "ab439a",
    "472e8c",
    "2a3390",
    "27337b",
    "0183c1",
    "00a885",
    "f3f0a1",
    "f7adbe",
    "f0b3d2",
    "bf9bcb",
    "70ceea",
    "a0d7d1",
    "9f968d",
    "0198cf",
    "57b948",
    "ffe709",
    "f9a05c",
    "f46d7b",
    "f04f9f",
    "ce4b9b",
    "85764f",
    "8c734b",
)


def hex2rgb(color):
    def hex2int(val):
        return int("0x" + val, 16)

    return f"{hex2int(color[0:2])},{hex2int(color[2:4])},{hex2int(color[4:6])}"


PALETTE = [hex2rgb(color) for color in PALETTE_HEX]
