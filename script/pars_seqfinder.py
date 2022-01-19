#!/usr/bin/env python3
# coding: utf-8

import regex
import pyfastx
import sys

# Parse arguments
fasta_file = sys.argv[1]
output_file = sys.argv[2]
# pattern = sys.argv[2]
# mismatch = sys.argv[3]
# insertion = sys.argv[4]
# deletion = sys.argv[5]
# errors = sys.argc[6]

# parS sites pattern values
pattern = "TGTTTCACGTGAAACA"

mismatch = 3
insertion = 0
deletion = 0
errors = max([mismatch, insertion, deletion])
all = False

# Read fasta file
def main(fasta_file, pattern, mismatch=0, insertion=0, deletion=0, errors=None):

    # Defined errors value.
    if errors == None:
        errors = max([mismatch, insertion, deletion])
    else:
        errors = max([mismatch, insertion, deletion, errors])

    # Build regex for the search
    pattern = regex.compile(
        "("
        + pattern
        + "){s<="
        + str(mismatch)
        + ",i<="
        + str(insertion)
        + ",d<="
        + str(deletion)
        + ",e<="
        + str(errors)
        + "}"
    )

    with open(output_file, "w") as out:
        for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
            start = 1
            res = regex.search(pattern, seq)
            while res != None:
                if (
                    res[0][1] == "G"
                    and res[0][2] == "T"
                    and res[0][3] == "T"
                    and res[0][5] == "C"
                    and res[0][6] in ["A", "T"]
                    and res[0][7] in ["C", "T"]
                    and res[0][8] in ["A", "G"]
                    and res[0][9] in ["G", "T"]
                    and res[0][10] in ["A", "G"]
                    and res[0][11] in ["A", "G"]
                    and res[0][12] == "A"
                    and res[0][13] in ["A", "T"]
                    and res[0][14] == "C"
                ):
                    out.write(
                        "\t".join(
                            [
                                name,
                                str(res.span()[0] + start),
                                str(res.span()[1] + start),
                                "parS",
                                res[0],
                                "\n",
                            ]
                        )
                    )
                elif (
                    res[0][1] == "G"
                    and res[0][2] == "T"
                    and res[0][3] == "T"
                    and res[0][5] == "C"
                    and res[0][12] == "A"
                    and res[0][14] == "C"
                    and all
                ):
                    out.write(
                        "\t".join(
                            [
                                name,
                                str(res.span()[0] + start),
                                str(res.span()[1] + start),
                                "parS-",
                                res[0],
                                "\n",
                            ]
                        )
                    )
                seq = seq[res.span()[1] :]
                start += res.span()[1]
                res = regex.search(pattern, seq)


main(fasta_file, pattern, mismatch, insertion, deletion)
