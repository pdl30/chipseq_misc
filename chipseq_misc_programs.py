#!/usr/bin/python

########################################################################
# 1 August 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import subprocess
import sys, os, re
import argparse
import pkg_resources

def convert_sam_bam(sam):
	name = re.sub(".sam", "", sam)
	command1 = "samtools view -bS {} > {}".format(sam, name+".bam")
	command2 = "samtools sort {} {}".format(name+".bam", name+"_sort")
	command3 = "samtools index {}".format(name+"_sort.bam")
	subprocess.call(command1, shell=True)
	subprocess.call(command2, shell=True)
	subprocess.call(command3, shell=True)

def ngs_plot(bam_file, size, region=None, bed=None):
	outname = re.sub("_sort", "", bam_file)
	outname = re.sub(".bam", "_ngsplot", outname)
	if not bed:
		command = "ngs.plot.r -G mm10 -R {3} -C {0} -O {1} -L {2} -FL 150".format(bam_file, outname, size, region)
	else:
		command = "ngs.plot.r -G mm10 -R bed -E {3} -C {0} -O {1} -L {2} -FL 150".format(bam_file, outname, size, bed)
	subprocess.call(command.split())

def convert_ucsc(ifile):
	outfile = re.sub(".bdg", "_ucsc.bdg", ifile)
	output = open(outfile, "w")
	with open(ifile) as f:
		for line in  f:
			line = line.rstrip()
			word = line.split("\t")
			if "GL" in word[0] or "MG" in word[0] or "JH" in word[0]:
				pass
			elif word[0] == "MT":
				chromo = "chrM"
				output.write("{}\t{}\t{}\t{}\n".format(chromo, word[1], word[2], word[3])),
			else:
				chromo = "chr" + word[0]
				output.write("{}\t{}\t{}\t{}\n".format(chromo, word[1], word[2], word[3])),
	output.close()

def macs2_normalised_track(sample, control, chrlen):
	name = re.sub("_sort", "", sample)
	name = re.sub(".bam", "", name)
	command = "macs2 callpeak -t {} -c {}  -B --nomodel --SPMR -g mm -n {}".format(sample, control, name)
	#subprocess.call(command.split())
	convert_ucsc("{0}_treat_pileup.bdg".format(name))
	convert_ucsc("{0}_control_lambda.bdg".format(name))
	command2 = "macs2 bdgcmp -t {0}_treat_pileup_ucsc.bdg -c {0}_control_lambda_ucsc.bdg -o {0}_FE.bdg -m FE".format(name)
	command3 = "macs2 bdgcmp -t {0}_treat_pileup_ucsc.bdg -c {0}_control_lambda_ucsc.bdg -o {0}_logLR.bdg -m logLR -p 0.00001".format(name)
	command4 = "bdg2bw {0}_FE.bdg {1}".format(name, chrlen)
	command5 = "bdg2bw {0}_logLR.bdg {1}".format(name, chrlen)
	print "{}\n{}\n{}\n{}\n{}\n".format(command, command2, command3, command4, command5)
	subprocess.call(command2.split())
	subprocess.call(command3.split())
	subprocess.call(command4.split())
	subprocess.call(command5.split())

def deeptools_normalised_bw(sample, control, chrlen, name):
	command = "bamCompare -b1 {} -b2 {} -o {}.bw".format(sample, control, name)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Random Programs, first is ngs.plot.R\n')
	subparsers = parser.add_subparsers(help='Programs included',dest="subparser_name")

	ngsplot_parser = subparsers.add_parser('ngsplot', help="Runs ngs.plot.R")
	ngsplot_parser.add_argument('-input', help='Input file, either sam or bam', required=False)
	ngsplot_parser.add_argument('-format', help='Options are sam/bam', required=False)
	ngsplot_parser.add_argument('-s', action='store_true', help='If selected, bam files will be sorted and indexed', required=False)
	ngsplot_parser.add_argument('-region', help='Region options are: tss, tes, genebody, exon, cgi, enhancer, dhs', default="tss", required=False)
	ngsplot_parser.add_argument('-bed', help='Custom region in bed format', required=False)
	ngsplot_parser.add_argument('-size', help='Size around region of interest, default=2000', default=2000, required=False)

	macs2_parser = subparsers.add_parser('macs2_track', help="Generates input normalised bedgraph")
	macs2_parser.add_argument('-s', help='Sample bam file', required=True)
	macs2_parser.add_argument('-c', help='Control bam file', required=True)
	macs2_parser.add_argument('-g', help='Genome, either mm10/hg19', required=True)
	
	args = vars(parser.parse_args())
	if args["subparser_name"] == "ngsplot":
		name = args["input"]
		if args["format"] == "sam":
			convert_sam_bam(args["input"])
			name = re.sub(".sam", "_sort.bam", args["input"])

		if args["bed"]:
			ngs_plot(name,args["size"], bed=args["bed"])
		else:
			ngs_plot(name, args["size"], region=args["region"])
	elif args["subparser_name"] == "macs2_track":
		path = pkg_resources.resource_filename('pychiptools', 'data/')
		chrom = path + "/{}.chrom.sizes".format(args["g"])
		macs2_normalised_track(args["s"], args["c"], chrom)