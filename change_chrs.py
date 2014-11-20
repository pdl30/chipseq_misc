#!/usr/bin/python

########################################################################
# 09 Oct 2014
# Patrick Lombard, Centre for Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, re, sys
import argparse
import subprocess
 
def convert_ucsc_ens(ifile, ofile):
	output = open(ofile, "w")
	with open(ifile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			chrom = word[0].strip("chr")
			del word[0]
			output.write("{}\t".format(chrom)),
			output.write('\t'.join(word)),
			output.write('\n'),
	output.close()

def convert_ens_ucsc(ifile, ofile):
	output = open(ofile, "w")
	with open(ifile) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			chrom = "chr"+word[0]
			del word[0]
			output.write("{}\t".format(chrom)),
			output.write('\t'.join(word)),
			output.write('\n'),
	output.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Changes UCSC to Ensembl format\n ')
	parser.add_argument('-i', '--input', help='Input file name', required=True)
	parser.add_argument('-e', action="store_true", help='Reverses default behaviour', required=False)
	parser.add_argument('-o', '--output', help='output file name', required=False)
	args = vars(parser.parse_args())
	if args["e"]:
		convert_ens_ucsc(args["input"], args["output"])
	else:
		convert_ucsc_ens(args["input"], args["output"])