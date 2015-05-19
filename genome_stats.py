#!/usr/bin/env python

import gzip
from Bio import SeqIO

#################################################
###           Parse command options           ###
#################################################

usage = """Usage: filter_parse.py [options]"""
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("-i", action="store", type = "string", dest = "input", help = "input file")
parser.add_option("-o", action="store", type = "string", dest = "output", help = "output report file")
parser.add_option("-n", action="store", type = "int", dest = "nlimit", help = "threshold number of Ns used to split scaffolds")

options, args = parser.parse_args()

def parser():
	file = open(options.input, "rU")
	for record in SeqIO.parse(file, "fasta"):
		return record
	

def gzparser():
	file = gzip.open(options.input, "r")
	for record in SeqIO.parse(file, "fasta"):
		return record

	
def bz2parser():
	file = bz2.open(options.input, "r")
	for record in SeqIO.parse(file, "fasta"):
		return record

def nsplit(record):
	window = options.nlimit
	sequence = record.seq
	for i in range(0, len(sequence)-window, window):
		testseq = sequence[i:i+len(window)]
		pos = i + len(window)
		if "A" or "C" or "G" or "T" not in testseq:
			
			

	
def main():
	infile = options.input
	if infile.endswith(".gz"):
		gzparser()
	elif infile.endswith(".bz2"):
		gz2parser()
	else:
		parser()
	nsplit(record)