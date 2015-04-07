#!/usr/bin/python

########################################################################
# 12 Jan 2015
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import re, os, sys
import argparse
import tempfile
import pysam
import ConfigParser

def find_peak_coverage(combined, samples, pvals, sizes):
	data = {}
	for pval in pvals:
		data[pval] = {}
		for sample in samples:
			data[pval][sample] = {}
			b = tempfile.NamedTemporaryFile(delete=False)
			b.close()
			command = "bedtools coverage -abam {} -b {} > {}\n".format(samples[sample], combined[pval], b.name)
			subprocess.call(command, shell=True)
			with open(b.name) as f:
				for line in f:
					line = line.rstrip()
					word = line.split("\t")
					norm = int(word[3])/float(sizes[sample])
					norm *= 1000000
					if (word[0], word[1], word[2]) in data[pval][sample]:
						data[pval][sample][(word[0], word[1], word[2])] += norm
					else:
						data[pval][sample][(word[0], word[1], word[2])] = norm
	return data

def sam_size(conditions):
	results=  {}
	for peak in conditions:
		if conditions[peak] not in results:
			size = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(conditions[peak]) ])
			results[conditions[peak]] = size
	return results

def combine_all_peaks(conditions):
	data = {}
	for pval in conditions:
		f = tempfile.NamedTemporaryFile(delete=False)
		for sample in samples:
			name = sample + "_" + pval + "_400bp.bed"
			with open(name) as k:
				for line in k:
					line = line.rstrip()
					word = line.split("\t")
					if word[0] == "chrM":
						new_chr = "MT"
					else:
						new_chr = word[0].strip("chr")
			
					f.write("{}\t{}\t{}\n".format(new_chr, word[1], word[2])),
		f.close()
		g = tempfile.NamedTemporaryFile(delete=False)
	#	h = tempfile.NamedTemporaryFile(delete=False)
		g.close()
	#	h.close()
		command1 = "sortBed -i {} > {}".format(f.name, g.name)
	#command2 = "mergeBed -i {} > {}".format(g.name, h.name)
		subprocess.call(command1, shell=True)
	#	subprocess.call(command2, shell=True)
		data[pval] = g.name
	return data

def write_results(coverage):
	for pval in coverage:
		output = open("{}_result.tsv".format(pval), "w")
		for key in sorted(coverage[pval]["Etv5_Flag_KI20_N16h_ucsc"]):
			output.write("{}\t{}\t{}".format(key[0], key[1], key[2])),
			for sample in sorted(coverage[pval]):
				output.write("\t{}".format(coverage[pval][sample][key])),
			output.write("\n"),

def ConfigSectionMap(section, Config):
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

def main():
	parser = argparse.ArgumentParser(description='Creates a table for coverage of a set of peaks. Still needs work\n')
	parser.add_argument('-c', 'config.ini', help='[Conditions] contains the peaks as keys and bam files as values')
	parser.add_argument('-o', '--output', help='Output file')
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])

	conditions = ConfigSectionMap(Config, "Conditions")

	sizes = sam_size(conditions)

	combined = combine_all_peaks(conditions)

	coverage = find_peak_coverage(combined, samples, pvals, sizes)
	write_results(coverage)

main()