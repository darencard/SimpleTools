#!/usr/local/env python

#print __name__

import optparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

usage_line = """
ddRAD_digest_coordinates.py

Version 1.0 (24 February, 2015)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

This script performs in in-silico double restriction enzyme digest on a genome sequence and reports the coordinates of digest points that result in fragments within a certain size range window. The script does not provide actual nucleotide fragments as output. Moreover, in cases where restriction sites are nearby (which occurs more with a common cutter), it reports coordinates for overlapping fragments, which would not normally occur in a real double digest RADseq library. Since this script is just meant to report the genomic location of fragments, this problem has been ignored, as a fragment that would normally be unobserved overlaps significantly with one that would have been observed.

User must provide a genome or reference sequence (in fasta format), the rare and common cutters (r1 and r2, respectively), the point in the restriction site after which the respective enzyme cut (c1 and c2, respectively), a lower and upper limit for fragment sizes, and an output file name. Output is a tab-delimited file with three columns: (1) the reference contig/scaffold ID, (2) the coordinate where the rare-cutter digests the DNA, and (3) the coordinate where the common-cutter digests the DNA that is

python ddRAD_digest_coordinates --input <reference_genome> --output <output_file_name> --r1 <rare_cutter> --c1 <rare_cut_location> --r2 <common_cutter> --c2 <common_cut_location --low <lower_fragment_limit> --up <upper_fragment_limit>"""


#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("--input", action = "store", type = "string", dest = "input", help = "input genome/sequence (in fasta format)")
parser.add_option("--output", action = "store", type = "string", dest = "output", help = "output file name for tab-delimited coordinate results")
parser.add_option("--r1", action = "store", type = "string", dest = "renz1", help = "first restriction enzyme (rare cutter)")
parser.add_option("--c1", action = "store", type = "string", dest = "cut1", help = "the base in the first restriction enzyme (rare cutter) after which the enzyme cuts (e.g., 6 in SbfI)")
parser.add_option("--r2", action = "store", type = "string", dest = "renz2", help = "second restriction enzyme (common cutter)")
parser.add_option("--c2", action = "store", type = "string", dest = "cut2", help = "the base in the second restriction enzyme (common cutter) after which the enzyme cuts (e.g., 1 in MspI)")
parser.add_option("--low", action = "store", type = "string", dest = "low", help = "lower limit for fragment size selection")
parser.add_option("--up", action = "store", type = "string", dest = "up", help = "upper limit for fragement size selection")

options, args = parser.parse_args()

##############################################################
### 	Perform Double Digest and Report Coordinates       ###
##############################################################

def double_digest(sequence, id, r_renz, r_cut_pos, c_renz, c_cut_pos, low, up, output):
	# Read contig
	DNA = sequence
	
	# Handles for various restriction enzymes, cut placements, and their complements for the rare-cutting enzyme
	rare_renz = Seq(r_renz, IUPAC.unambiguous_dna)
	rare_compmotif = str(rare_renz.complement()) # must also search for complement
	rare_motif = str(rare_renz)
	rare_cut_motif = int(r_cut_pos)
	rare_cut_compmotif = int(len(rare_motif))-int(r_cut_pos)
	
	# Handles for various restriction enzymes, cut placements, and their complements for the common-cutting enzyme
	common_renz = Seq(c_renz, IUPAC.unambiguous_dna)
	common_compmotif = str(common_renz.complement()) # must also search for complement
	common_motif = str(common_renz)
	common_cut_motif = int(c_cut_pos)
	common_cut_compmotif = int(len(common_motif))-int(c_cut_pos)
	
	# For bp in range from 0 to length of sequence - RE motif length, iterating by 1bp
	for i in range(0, len(DNA)-len(rare_motif), 1):
		# Rare test sequence is i + length of rare RE motif
		rare_testseq = str(DNA[i:i+len(rare_motif)])
		rare_pos = i+1
		if rare_testseq == rare_motif:
		# if rare enzyme test sequence equals the rare restriction enzyme motif, report position as position in loop + cut location in enzyme
			rare_digest = rare_pos + rare_cut_motif
			# whenever there is a rare enzyme cut, scan a window of basepairs upstream (based on lower/upper limits designed) for a common enzyme cut
			for j in range(rare_digest+int(low), (rare_digest+int(up))-len(rare_motif),1):
				common_testseq = str(DNA[j:j+len(common_motif)])
				common_pos = j+1
				# if common enzyme test sequence equals the common restriction enzyme motif, report position
				if common_testseq == common_motif:
					common_digest = common_pos + common_cut_motif
					if common_digest < len(DNA):
						rare_j_line = id+"\t"+str(rare_digest)+"\t"+str(common_digest)+\t+"+\n"
						output.write(rare_j_line)
						print rare_j_line
			# whenever there is a rare enzyme cut, scan a window of basepairs downstream (based on lower/upper limits designed) for a common enzyme cut
#			for k in range(rare_digest-int(up), (rare_digest-int(low))-len(rare_motif),1):
#				common_testseq = str(DNA[k:k+len(common_motif)])
#				common_pos = k+1
#				if common_testseq == common_motif:
#					common_digest = common_pos + common_cut_motif
#					if common_digest > 0:
#						rare_k_line = id+"\t"+str(rare_digest)+"\t"+str(common_digest)+"\n"
#						output.write(rare_k_line)
#						print rare_k_line
		elif rare_testseq == rare_compmotif:
		# must do the same as above but with complement sequences (for opposite strand)
			rare_digest = rare_pos + rare_cut_compmotif
			# whenever there is a rare enzyme cut, scan a window of basepairs upstream (based on lower/upper limits designed) for a common enzyme cut. This actually ends up being downstream on the strand we care about.
#			for j in range(rare_digest+int(low), (rare_digest+int(up))-len(rare_motif),1):
#				common_testseq = str(DNA[j:j+len(common_compmotif)])
#				common_pos = j+1
#				if common_testseq == common_compmotif:
#					common_digest = common_pos + common_cut_compmotif
#					if common_digest < len(DNA):
#						rare_j_line = id+"\t"+str(rare_digest)+"\t"+str(common_digest)+"\n"
#						output.write(rare_j_line)
#						print rare_j_line
			# whenever there is a rare enzyme cut, scan a window of basepairs downstream (based on lower/upper limits designed) for a common enzyme cut. This actually ends up being upstream on the strand we care about.
			for k in range(rare_digest-int(up), (rare_digest-int(low))-len(rare_compmotif),1):
				common_testseq = str(DNA[k:k+len(common_compmotif)])
				common_pos = k+1
				if common_testseq == common_compmotif:
					common_digest = common_pos + common_cut_compmotif
					if common_digest > 0:
						rare_k_line = id+"\t"+str(rare_digest)+"\t"+str(common_digest)+\t+"-\n"
						output.write(rare_k_line)
						print rare_k_line

output = open(options.output, "w") # open output file

# open fasta input file and loop through sequences (contigs/scaffolds/etc.)
for sequence in SeqIO.parse(options.input, 'fasta'):
	double_digest(sequence.seq, sequence.id, options.renz1, options.cut1, options.renz2, options.cut2, options.low, options.up, output)

output.close()				