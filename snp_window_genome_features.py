#!/usr/bin/env python

import optparse

usage_line = """
snp_window_genome_features.py

Version 1.0 (28 July, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

Script that takes a list of SNPs (scaffold ID and SNP position, tab delimited) and a genome GFF file \
and extractions genomic features within a designed window around the SNP. Useful for pulling out genes \
located in the region surrounding an interesting SNP (inferred from genome scans). Input is a \
tab-delimited file with scaffold ID in the first column and SNP position in the second column, \
a matching GFF file from the genome of interest, and a window size (# bp on each side of SNP. \
User can also designate an output file, where GFF lines from features within the window will be \
written. Note that this script matches chromosome/scaffold IDs using a "startswith" function that \
matches to the beginning of the element in the GFF file, so appropriate changes need to be made to \
the SNP input or to this script if other circumstances exist.

python snp_window_genome_features.py --snps <snps.txt> --gff <genome.gff> [--out <output.gff> --window \
<window_size>]
"""

#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)

parser.add_option("--snps", action="store", type = "string", dest = "snps", help = "tab-delimited file of chromosome and SNP position")
parser.add_option("--gff", action="store", type = "string", dest = "gff", help = "full GFF annotation file")
parser.add_option("--out", action="store", type = "string", dest = "out", help = "output file prefix to write GFF and summary output [filtered.out]", default = "filtered.out")
parser.add_option("--window", action="store", type = "string", dest = "window", help = "window size [10000]", default = "10000")

options, args = parser.parse_args()

#################################################
###         Find Matches within Window        ###
#################################################

def find_matches(gff, scaff, pos):
	foo = gff.rstrip().split("\t")
	window_low = int(pos) - int(options.window)
	window_high = int(pos) + int(options.window)
	if foo[0].startswith(scaff):
		if window_low <= int(foo[3]) <= window_high or window_low <= int(foo[4]) <= window_high:
			summary_line = str(scaff)+"\t"+str(pos)+"\t"+str(window_low)+"\t"+str(window_high)+"\t"+str(gff)
			gff_line = str(gff)
			summary_out.write(summary_line)
			gff_out.write(gff_line)
			print summary_line


#################################################
###              Parse SNPs file              ###
#################################################

snp_list = []				
for line in open(options.snps, "r"):
	if not line.strip().startswith("#"):
		bar = line.rstrip().split("\t")
		snp_list.append([bar[0], bar[1]])


#################################################
###    Parse GFF and call Window Matching     ###
#################################################

summary_out = open(options.out+".summary.txt", "w")
gff_out = open(options.out+".gff", "w")
for row in open(options.gff, "r"):
	if not line.strip().startswith("#"):
		for item in snp_list:
			find_matches(row, item[0], item[1])
summary_out.close()
gff_out.close()