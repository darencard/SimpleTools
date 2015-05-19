#!/usr/local/env python

#print __name__

import optparse
from Bio import SeqIO
import collections

usage_line = """
ddRAD_digest_coordinates.py

Version 1.0 (24 February, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

This script "re-coordinates" a genome, meaning that it takes coordinates set individually on different chromosomes/contigs and adjust the coordinate values to a new, genome-wide scale. It therefore puts coordinates on an absolute genome scale rather than a relative chromosome scale. This is useful if one want to plot features across a whole genome on the same plot.

User must provide the reference genome/sequence so that contig lengths can be extracted and a tab-delimited text file that includes chromosome/contig/scaffold IDs in one column (which correspond to the reference sequence) and relative coordinates in another column. User must specify which columns those values are found.

Output is a text file with the new coordinates in the first column, the chromosome ID in the second column, the chromosome coordinate in the third column, the total chromosome length in the fourth column, and the running chromosome length in the fifth column. Need to expand this script to also take an input value and output it (e.g., a p-value from a GWAS).

python recoord.py --ref <reference_genome> --track <input_coords_file> --chr <chr_column> --coord <coord_column> --output <output_file>"""


#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("--ref", action = "store", type = "string", dest = "ref", help = "input genome/sequence (in fasta format)")
parser.add_option("--track", action = "store", type = "string", dest = "track", help = "input coordinate track with chromosome and coordinate columns")
parser.add_option("--chr", action = "store", type = "string", dest = "chr", help = "column containing chromosome IDs (1, 2, ..., N)")
parser.add_option("--coord", action = "store", type = "string", dest = "coord", help = "column containing coordinates (1, 2, ..., N)")
parser.add_option("--output", action = "store", type = "string", dest = "outfile", help = "the output file to write new coordinates")

options, args = parser.parse_args()



#################################################
###    	Parse Reference Contig Lengths        ###
#################################################

def make_chr_list(reference, dictionary, running_length, counter):
	chr_dict = dictionary
	run_len = int(running_length)
	counter = counter
	for seq in SeqIO.parse(reference, 'fasta'):
		chrom_len = len(seq.seq)
		run_len += int(chrom_len)
		counter += 1
#		print counter
#		print chrom_len
#		print run_len
		chr_dict[seq.id] = run_len
	return chr_dict
	return running_length



#################################################
###  Adjust Coordinates to New Genome Scale   ###
#################################################
	
def re_coord(dict, infile, chromosome, coordinate, outfile):
	chrom_dict = dict
	input = infile
	chr = int(chromosome)-1
	coord = int(coordinate)-1
	output = outfile
	output.write("new_coord\told_chr\told_coord\tlen_chr\tlen_past_chr\n")
	
	for row in open(input).read().splitlines():
		foo = row.split()
#		print foo
		chrom = foo[chr]
		coor = foo[coord]
#		print "Chromosome: "+chrom
#		print "Chromosome index: "+str(chrom_dict.keys().index(str(chrom)))
		if int(chrom_dict.keys().index(str(chrom))) == 0:
			new_coord = coor
#			print "Adjusted coordinate: "+str(new_coord)
			out = str(new_coord)+"\t"+str(chrom)+"\t"+str(coor)+"\t"+str(chrom_dict[chrom])+"\t0\n"
			output.write(out)
			print out
		else:
			idx_last_chr = chrom_dict.keys().index(str(chrom))
#			print "Index of last chromosome: "+str(idx_last_chr)
			len_past_chr = chrom_dict[chrom_dict.keys()[idx_last_chr - 1]]
#			print "Running length of all previous chromosomes: "+str(len_past_chr)
			new_coord = int(len_past_chr) + int(coor)
#			print "Coordinate on this chromosome: "+str(int(coor))
#			print "New, adjusted coordinate: "+str(new_coord)
			out = str(new_coord)+"\t"+str(chrom)+"\t"+str(coor)+"\t"+str(chrom_dict[chrom] - len_past_chr)+"\t"+str(len_past_chr)+"\n"
			print out
			output.write(out)
	
	

#################################################
###           	Main Program 		          ###
#################################################

chr_dict = collections.OrderedDict()

total_chr_len = 0
counter = 0

make_chr_list(options.ref, chr_dict, total_chr_len, counter)
#print chr_dict

new_coords = open(options.outfile, "w")

re_coord(chr_dict, options.track, options.chr, options.coord, new_coords)

new_coords.close()