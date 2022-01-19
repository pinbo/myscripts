#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys

settings_common = "PRIMER_TASK=check_primers" + "\n" + \
        "PRIMER_FIRST_BASE_INDEX=1" + "\n" + \
        "PRIMER_PICK_ANYWAY=1" + "\n" + \
        "PRIMER_SALT_DIVALENT=1.5" + "\n" # for PACE final MgCl2=2.2, default is 1.5

# input: 3 columns sep by tab: primer name, left primer and right primer
primerFile = sys.argv[1]
outFile = sys.argv[2]
p3input = open(outFile, "w")
with open(primerFile) as infile:
    for line in infile:
        line = line.strip()
        if line:
            primerID, leftPrimer, rightPrimer = line.split("\t")
            settings = settings_common + \
            "SEQUENCE_ID=" + primerID + "\n" + \
            "SEQUENCE_PRIMER=" + leftPrimer + "\n" + \
            "SEQUENCE_PRIMER_REVCOMP=" + rightPrimer + "\n" + \
            "=\n"
            p3input.write(settings)

p3input.close()