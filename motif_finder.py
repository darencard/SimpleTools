#!/usr/local/env python

#print __name__

import optparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

usage_line = """
motif_finder.py

Version 1.0 (24 March, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

This script searches a given reference DNA sequence for a specified DNA motif and reports the coordinates where a given motif is found in the reference sequence (coordinate based upon the first position of the motif). This works on both DNA strands and indicates the (+) or (-) strand in the output.

User must specify a genome or reference sequence (in fasta format), the sequence of the motif they are interested in (all caps), and an output file where the results will be written in tab-delimited format.

Note that there is likely a faster way of doing this (optimized function), but this should work.

python motif_finder.py --input <reference_genome> --output <output_file_name> --seq <query_sequence>"""


#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("--input", action = "store", type = "string", dest = "input", help = "input genome/sequence (in fasta format)")
parser.add_option("--output", action = "store", type = "string", dest = "output", help = "output file name for tab-delimited coordinate results")
parser.add_option("--seq", action = "store", type = "string", dest = "seq", help = "the query DNA motif you want to find")

options, args = parser.parse_args()

###################################################
### 	Find Motif and Report Coordinates       ###
###################################################

def motif_finder(reference, ref_ids, query_seq, output):
	# Read contig
	DNA = reference
	# Read query
	query = Seq(query_seq, IUPAC.unambiguous_dna)
	# Make query a str and get str of complement (for minus strand)
	query_reg = str(query)
	query_comp = str(query.complement()) # must also search for complement
		
	# For bp in range from 0 to length of sequence - RE motif length, iterating by 1bp
	for i in range(0, len(DNA)-len(query), 1):
		# Rare test sequence is i + length of rare RE motif
		testseq = str(DNA[i:i+len(query)])
		pos = i+1
		# If test sequence equals query sequence (plus strand), print line to terminal and output file
		if testseq == query_reg:
			# out format = tab separated columns of sequence/contig id, position/coordinate, and + for strand
			line = str(ref_ids)+"\t"+str(pos)+"\t+\n"
			output.write(line)
			print line
		
		# If test sequence equal complement of query sequence (minus strand), print line to terminal and output file	
		elif testseq == query_comp:
			# out format = tab separated columns of sequence/contig id, position/coordinate, and - for strand
			line = str(ref_ids)+"\t"+str(pos)+"\t-\n"
			output.write(line)
			print line

outfile = open(options.output, "w") # open output file

# open fasta input file and loop through sequences (contigs/scaffolds/etc.)
for sequence in SeqIO.parse(options.input, 'fasta'):
	motif_finder(sequence.seq, sequence.id, options.seq, outfile)

outfile.close()				