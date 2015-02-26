#!/usr/bin/env python

import optparse
from Bio import SeqIO

usage_line = """ test """

usage = usage_line

parser = optparse.OptionParser(usage=usage)
parser.add_option("-i", action = "store", type = "string", dest = "input", help = "input file name")
parser.add_option("-o", action = "store", type = "string", dest = "output", help = "output file name")
parser.add_option("--min", action = "store", type = "string", dest = "min", help = "minimum fragment length")
parser.add_option("--max", action = "store", type = "string", dest = "max", help = "maximum fragment length")

options, args = parser.parse_args()

outfile = open(options.output, "w")

for sequence in SeqIO.parse(options.input, "fasta"):
	if len(sequence.seq) >= options.min and len(sequence.seq) <= options.min:
		out = ">" + sequence.id + "\n" + sequence.seq + "\n"
		outfile.write(out)

outfile.close()