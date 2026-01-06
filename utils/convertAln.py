#! /usr/bin/python

# SPDX-FileCopyrightText: 2026 Julien Y. Dutheil <jy.dutheil@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Convert an alignment into a single-block MAF file.
import sys
from Bio import AlignIO
alignment = AlignIO.read(open(sys.argv[1] + ".fasta"), "fasta")
print("Alignment length %i" % alignment.get_alignment_length())
# Convert to maf:
maf = open(sys.argv[1] + ".maf", 'w')
maf.write("##maf version=1\n")
maf.write("#Generated from a single alignment\n\n")
maf.write("a\n")
for record in alignment :
    n = sum([1 for x in record.seq if x != '-'])
    maf.write("s %s.1\t0\t%d\t+\t%d\t%s\n" % (record.id, n, n, record.seq))
maf.close()
