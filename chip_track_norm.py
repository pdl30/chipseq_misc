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
import pkg_resources
import pychiptools

def genomeCoverage(name, house=None, cov=None):
	print "==> Converting bed to bedGraph...\n"
	inbed = pybedtools.BedTool(name+"_ucsc.BED")
	if house:
		outcov = inbed.genome_coverage(bg=True, genome='mm10', scale=house)
		output = name+"_house.bedGraph"
		outcov.saveas(output)
	elif cov:
		outcov = inbed.genome_coverage(bg=True, genome='mm10', scale=cov)
		output = name+"_cov.bedGraph"
		outcov.saveas(output)
	return output

def bedgraphtobigwig(bedgraph, chrom):
	bw = re.sub(".bedGraph$", ".bw", bedgraph)
	print "==> Converting bedGraph to bigWig...\n"
	command = ["bedGraphToBigWig", bedgraph, chrom, bw]
	subprocess.call(command)

def normalise_to_housekeeper(count_file):
	#print "==> Normalising to Housekeeper...\n"
	with open(count_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			if word[0] == "ENSMUSG00000057666": #Gapdh, substitute with what you want to use. REmove from production?
				housekeeper = int(word[1])	
	return housekeeper

def normalise_to_coverage(cov_file):
	#print "==> Normalising to Coverage file...\n"
	total_coverage = 0
	with open(cov_file) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			coverage = int(word[-4])
			total_coverage += coverage
	return total_coverage

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Processes ChIP-seq samples to bigWig tracks.\n')
	parser.add_argument('-i','--input', help='Input BED file in UCSC format', required=False)
	parser.add_argument('-g','--genome', help='Genome the samples are aligned to, options include mm10/mm9/hg19', required=True)
	parser.add_argument('-a','--house', help='Housekeeper normalisation. Input file is HTSEQ-count file containing gene for normalisation on first line', required=False)
	parser.add_argument('-c', '--cov', help='Nomralisation using bedtools coverage file', required=False)
	args = vars(parser.parse_args())
	chrom = pkg_resources.resource_filename('pychiptools', 'data/{}.chrom.sizes'.format(args["genome"]))
	
	path0 = os.getcwd()
	name = re.sub("_ucsc.BED$", "", args["input"])

	if args["house"]:
		house = normalise_to_housekeeper(args["house"])
		scale = float(1000)/int(house) #Works and checked
		bedgraph = genomeCoverage(name, house=scale)
		bedgraphtobigwig(bedgraph, chrom)
	elif args["cov"]:
		cov = normalise_to_coverage(args["cov"])
		scale = float(1000000)/int(cov) #Works and checked
		bedgraph = genomeCoverage(name, cov=scale)
		bedgraphtobigwig(bedgraph, chrom)