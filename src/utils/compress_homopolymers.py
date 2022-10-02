import sys
from itertools import groupby


def compress_homopolymers(seq):

    return "".join([k for k, g in groupby(seq)])


COMPRESS = sys.argv[1] == "compressed"

while True:
    name = sys.stdin.readline()

    if name == "":
        break

    sequence = sys.stdin.readline()
    print(name, end="")

    if COMPRESS:
        print(compress_homopolymers(sequence), end="")
    else:
        print(sequence, end="")
