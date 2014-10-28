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
from collections import defaultdict
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector, FloatVector, StrVector

def annotate_ensembl(dict_obj):
	ens = importr("biomaRt")
	ensembl = ro.r.useMart("ensembl")
	genome="mmusculus_gene_ensembl"
	ensembl = ro.r.useDataset(genome, mart=ensembl)
	values = []
	for key1 in dict_obj.keys():
		values.append(key1)
	C1BM = ro.r.getBM(attributes=StrVector(["ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "description", "gene_biotype"]), 
		filters="ensembl_gene_id", values=values, mart=ensembl)
	gene = list(C1BM.rx(True,1))
	chr1 = list(C1BM.rx(True,2))
	tss = list(C1BM.rx(True,3))
	end = list(C1BM.rx(True,4))
	st = list(C1BM.rx(True,5))
	name = list(C1BM.rx(True,6))
	des = list(C1BM.rx(True,7))
	bio = list(C1BM.rx(True,8))
	data = {}
	for index, g in enumerate(gene):
		data[g] = (chr1[index], tss[index], end[index], st[index], name[index], des[index], bio[index])
	return data
def reverse_dict(idict):
	inv_map = {}
	for k, v in idict.iteritems():
		inv_map[v] = inv_map.get(v, [])
		inv_map[v].append(k)
	return inv_map

def postprocess_diffreps(cond1, cond2, gtf):
	resbed = "{}_vs_{}_diffReps.txt".format(cond1, cond2)
	output = open("tmp1.bed", "w")
	with open(resbed) as f:
		for line in f:
			if line.startswith("#") or line.startswith("Chrom"):
				pass
			else:
				line = line.rstrip()
				word = line.split("\t")
				if float(word[12]) < 1e-7:
					name = word[0].strip("chr")
					output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name, word[1], word[2],word[6], word[7], word[10], word[11], word[12], word[13])),
	output.close()
	command = "windowBed -w 50000 -a {} -b {} > {}".format("tmp1.bed", gtf, "tmp2.bed")
	subprocess.call(command, shell=True)
	data = defaultdict(list)
	ens = {}
	with open("tmp2.bed") as f:
		for line in f:
			line = line.rstrip()
			word = line.split("\t")
			data[(word[0], word[1], word[2])].append(line)
			info = word[16].split(";")
			ensgene = info[0].strip("gene_id ")
			ensgene = ensgene.strip("\"")
			ens[ensgene] = 1
	anno = annotate_ensembl(ens)
	output1 = open("{}_vs_{}_diffReps_result.txt".format(cond1, cond2), "w")
	output1.write("Chrom\tStart\tEnd\tPvalue\tLFC\tNearest Genes\n")
	for key in data:
		gene_list = []
		word2 = data[key][0].split("\t")
		output1.write("{}\t{}\t{}\t{}\t{}\t".format(key[0], key[1], key[2], word2[7], word2[6])),
		for v in data[key]:
			words = v.split("\t")
			info = words[16].split(";")
			gene = info[0].strip("gene_id ")
			gene = gene.strip("\"")
			if gene in anno:
				gene_list.append(anno[gene][4])
		gene_list = set(gene_list)
		gene_list = list(gene_list)
		for gene in gene_list:
			output1.write("{},".format(gene)),
		output1.write("\n"),
	output1.close()

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

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='postprocess_diffreps\n')
	parser.add_argument('-c', '--config', help='ConfigParser input', required=True)
	parser.add_argument('-g', '--gtf', help='Input GTF', required=True)
	args = vars(parser.parse_args())
	conditions, controls, comparisons = get_config_args(args, ["Conditions", "Controls", "Comparisons"])
	inv_conds = reverse_dict(conditions)
	for comp in comparisons:
		c = comparisons[comp].split(",")
		comps = [x.strip(' ') for x in c]
		postprocess_diffreps(comps[0], comps[1], args["gtf"])