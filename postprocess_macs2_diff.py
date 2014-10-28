#!/usr/bin/python

########################################################################
# 1 August 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import os, sys, re
import subprocess
import argparse
import ConfigParser

def remove_chr(ibed1, cond1, ibed2, cond2, obed):
	output = open(obed, "w")
	with open(ibed1) as f:
		header = next(f)
		for line in f:
			line = line.rstrip()
			word= line.split("\t")
			name = re.sub("chr", "", word[0])
			output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(name, word[1], word[2], word[3], word[4], cond1)),
	with open(ibed2) as f:
		header = next(f)
		for line in f:
			line = line.rstrip()
			word= line.split("\t")
			name = re.sub("chr", "", word[0])
			output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(name, word[1], word[2], word[3], word[4], cond2)),
	output.close()
#Process closestBed output for macs2 diff
def process_close_out(ifile, outfile):
	output = open(outfile, "w")
	output.write("Chromosome\tStart\tEnd\tDiff\tHigher Sample\tEnsemblID\tGene Name\tDistance to Gene\tGene Start\tGene End\n"),
	with open(ifile) as f:
		header = next(f)
		for line in f:
			line = line.rstrip()
			word= line.split("\t")
			info = word[14].split(";")
			gene_id = info[0].split(" ")
			gene = gene_id[1].strip("\"")
			gene_name = info[3].split(" ")
			gene_name = gene_name[2].strip("\"")
			output.write("chr{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(word[0], word[1], word[2], word[4], word[5], gene, gene_name, word[15], word[9], word[10])),
	output.close()

def postprocess_bgdiff(cond1, cond2, gtf):
	bed1 = "{}_vs_{}_c3.0_cond1.bed".format(cond1, cond2)
	bed2 = "{}_vs_{}_c3.0_cond2.bed".format(cond1, cond2)
	remove_chr(bed1, cond1, bed2, cond2,  "tmp1.bed")
	command1 = "closestBed -D a -a {} -b {} > {}".format("tmp1.bed", gtf, "tmp2.bed")
	subprocess.call(command1, shell=True)
	process_close_out("tmp2.bed", "{}_vs_{}_closest.tsv".format(cond1, cond2))

def ConfigSectionMap(Config, section):
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

def get_config_args(args, arglist):
	Config = ConfigParser.ConfigParser()
	Config.optionxform = str
	Config.read(args["config"])
	dict_list = []
	for arg in arglist:
		dict_list.append(ConfigSectionMap(Config, arg))
	return dict_list

def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='postprocess_diffreps\n')
	parser.add_argument('-c' '--config', help='ConfigParser input', required=True)
	parser.add_argument('-g' '--gtf', help='Input GTF', required=True)
	args = vars(parser.parse_args())
	conditions, controls, comparisons = get_config_args(args, ["Conditions", "Controls", "Comparisons"])
	inv_conds = reverse_dict(conditions)
	for comp in comparisons:
		c = comparisons[comp].split(",")
		comps = [x.strip(' ') for x in c]
		postprocess_bgdiff(comps[0], comps[1], args["gtf"])