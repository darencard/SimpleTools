#!/usr/bin/env python

#print __name__

import os
import optparse

usage_line = """
delimiter_swap.py

Version 1.0 (28 August, 2014)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

A simple script that converts comma-separated value (csv) files into tab-separated value (tsv) files. \
Uses simple sed find-and-replace unix command.

python delimiter_swap.py -i <input> -o <output> -t <csv/tsv>"""

#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("-i", action="store", type = "string", dest = "input", help = "input file")
parser.add_option("-o", action="store", type = "string", dest = "output", help = "output file")
parser.add_option("-t", action="store", type = "string", dest = "type", help = "input file type (csv or tsv)")

options, args = parser.parse_args()


#################################################
###           Convert CSV to TSV	          ###
#################################################

def csv2tsv():
	if options.input is None:
		print "\n***Error: specify input file!***\n"
	elif options.output is None:
		print "\n***Error: specify output file!***\n"
	elif options.type is None:
		print "\n***Error: specify input file type (csv or tsv)!***\n"
	else:
		if options.type == "csv":
			print "\n***Converting from comma-separated to tab-separated***\n"
			os.system("sed 's/,/\t/g' "+options.input+" > "+options.output)
		elif options.type == "tsv":
			print "\n***Converting from tab-separated to comma-separated***\n"
			os.system("sed 's/\t/,/g' "+options.input+" > "+options.output)
		else:
			print """\n***Error: file type must be "csv" or "tsv"!***\n"""
	
		
csv2tsv()
	