#!/usr/bin/env python

#print __name__

import os
import optparse

usage_line = """
remove_duplicate_rows.py

Version 1.0 (28 August, 2014)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or desirable functioning.

A simple script that removes duplicate rows based upon a designated column number (1, 2, ..., n). \
Passing the designated column option of 0 will filter rows that are completely identical (i.e., \
not by column). Uses a Unix awk command and both tab separated value (tsv) and comma separated \
value (csv) files are accepted.

python remove_duplicate_rows.py -i <input_csv> -o <output_tsv> -c <0,1,2...>"""

#################################################
###           Parse command options           ###
#################################################

usage = usage_line
                        
parser = optparse.OptionParser(usage=usage)
parser.add_option("-i", action="store", type = "string", dest = "input", help = "input file")
parser.add_option("-o", action="store", type = "string", dest = "output", help = "output file")
parser.add_option("-c", action="store", type = "string", dest = "column", help = "column to use when filtering duplicates")

options, args = parser.parse_args()


#################################################
###       Filter by specified column	      ###
#################################################

def rmdups():
	if options.input is None: 
		print "\n***Error: specify input file!***\n"
	elif options.output is None:
		print "\n***Error: specify output files!***\n"
	elif options.column is None:
		print "\n***Error: specify column to use for filtering duplicates!***\n"
	else:
		if options.column == "0":
			column = "entire row (or all columns)"
		else:
			column = "column "+options.column
		print "\n***Removing duplicate rows based upon "+str(column)+"***\n"
		line = str("awk '{ if ($"+options.column+" in stored_lines) x="+options.column+"; else print; stored_lines[$"+options.column+"]="+options.column+" }' "+options.input+" > "+options.output)
		os.system(line)
		
rmdups()
	