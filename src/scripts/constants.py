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
    "ffedbb",
    "ffe5ca",
    "fcdcc5",
    "fdd0af",
    "fcd3c1",
    "fbd3d3",
    "e7bad7",
    "c7b1d5",
    "adb9df",
    "adc5e7",
    "bae4f0",
    "bce4e6",
    "aedfdc",
    "bee3d2",
    "e0eed4",
    "e2ecaf",
    "ffe292",
    "ffdaa3",
    "fcc79b",
    "f8a980",
    "f8aa96",
    "f8a19a",
    "d991bf",
    "9d86be",
    "6f85c1",
    "7da7d9",
    "80d2e3",
    "87d1d0",
    "7bcbbe",
    "8dceb6",
    "c2e0ae",
    "dde89a",
    "ffdc7e",
    "fec47a",
    "f9a870",
    "f9966c",
    "f38c76",
    "f37a6f",
    "d881b6",
    "816ab0",
    "5564ab",
    "5f89c9",
    "4ec6de",
    "57c4c7",
    "5ac5b1",
    "65c295",
    "add68a",
    "d7e37d",
)


def hex2rgb(color):
    def hex2int(val):
        return int("0x" + val, 16)

    return f"{hex2int(color[0:2])},{hex2int(color[2:4])},{hex2int(color[4:6])}"


PALETTE = [hex2rgb(color) for color in PALETTE_HEX]
