#!/usr/bin/python

########################################################################
# 15 May 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, os, re
import pybedtools
import pysam
import argparse
import operator

def combine_sam_files(list_of_sams, outname):
	count = 0
	outsam = outname + "/" + outname + ".sam"
	print "==> Combining sam files...\n"
	for sam in list_of_sams:
		original_file = pysam.Samfile(sam)
		if count == 0:
			new_file = pysam.Samfile(outsam, mode='wh', template=original_file)
			for read in original_file: 
				new_file.write(read)
		else:
			for read in original_file: 
				new_file.write(read)
		count += 1

if __main__ == "__main__":
	parser = argparse.ArgumentParser(description='Combines 2 or more ChIP-seq sam files.\n')
	parser.add_argument('-r','--COMBINE', help='This takes 2 or more sam files and combines them and processes them to bigwigs. Output name is specified by -o', nargs="+", required=False)
	parser.add_argument('-o','--OUTNAME', help='Only use if using the option -r! This will create a directory and files with this name.', required=False)
	args = vars(parser.parse_args())
	combine_sam_files(args["COMBINE"], args["OUTNAME"])