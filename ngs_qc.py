#!/usr/bin/python

########################################################################
# 08 Oct 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import argparse
import subprocess
import os, re, sys
from multiprocessing import Pool
import ConfigParser
import itertools

def change_bed(idir, threads):
	output=  open("{0}/{0}_htseq.BED".format(idir), "w")
	with open("{0}/{0}_ucsc.BED".format(idir)) as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			width = int(word[2]) - int(word[1])
			output.write("{}\t{}\t{}\t{}\t{}\n".format(word[0], word[1], word[2], width, word[5]))
	output.close()

def htseqtools(conditions):
	rscript = "library(htSeqTools)\n"
	for key in conditions:
		rscript += "{1} <- read.table('{0}/{0}_htseq.BED')\n".format(key, conditions[key])
		rscript += "colnames({}) <- c('space', 'start', 'end', 'width', 'strand')\n".format(conditions[key])
		rscript += "{0} <- RangedData({0})\n".format(conditions[key])

	rscript += "htSample2 <- RangedDataList("
	c = 0
	for key in conditions:
		if c == 0:
			rscript += "{0}={0}".format(key)
		else:
			rscript += ",{0}={0}".format(key)
		c += 1
	rscript += ")\n"
	rscript += "cmds1 <- cmds(htSample2,k=2, mc.cores={})\n".format(threads)

	return rscript

def function1(args):
	return change_bed(*args)

def ConfigSectionMap(section):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='QC results from htSeqTools package. Not even close to finished!\n')
	parser.add_argument('-config', help='Config file containing parameters, please see documentation for usage!', required=True)
	parser.add_argument('-threads', help='Number of threads to use, default=4', default=4, required=True) #Careful here!! v.v. high memory usage
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])

	conditions = ConfigSectionMap("Conditions")
	pool = Pool(8)
	keys = list(conditions.keys())
	#pool.map(function1, itertools.izip(keys, itertools.repeat(int(args["threads"]))))
	#pool.close()
	#pool.join()
	rscript = htseqtools(conditions, int(args["threads"]))
	rcode = open("htseqtools.R", "w")
	rcode.write(rscript)
	rcode.close()