#!/usr/bin/python

########################################################################
# 10 July 2014
# Patrick Lombard, Centre for Stem Stem Research
# Core Bioinformatics Group
# University of Cambridge
# All right reserved.
########################################################################

import sys, re, os
import csv
from collections import defaultdict
import subprocess
from operator import itemgetter
import HTSeq
import argparse

def chippeakanno(peaks, output_prefix):
	rscript = "library('ChIPpeakAnno')\n"
	rscript += "library(biomaRt)\n"
	rscript += "ensembl = useMart('ensembl')\n"
	rscript += "ensembl = useDataset('mmusculus_gene_ensembl', mart=ensembl)\n"
	rscript += "TSS.mouse.NCBIM38 = getAnnotation(mart=ensembl, featureType='TSS')\n"
	rscript += "peaks <- BED2RangedData('{}')\n".format(peaks)
	rscript += "annotatedPeaktss = annotatePeakInBatch(peaks, AnnotationData=TSS.mouse.NCBIM38, PeakLocForDistance='middle')\n"
	rscript += "annotatedPeaktrans = annotatePeakInBatch(peaks, AnnotationData=TSS.mouse.NCBIM38, FeatureLocForDistance='middle', PeakLocForDistance='middle', output='shortestDistance')\n"
	rscript += "annotatedPeaktss <- as.data.frame(annotatedPeaktss)\n"
	rscript += "annotatedPeaktrans <- as.data.frame(annotatedPeaktrans)\n"
	rscript += "C1BM <- getBM(c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'),filters = 'ensembl_gene_id', values = annotatedPeaktrans$feature, mart = ensembl)\n"
	rscript += "C2BM <- getBM(c('ensembl_gene_id', 'external_gene_name','description', 'gene_biotype'),filters = 'ensembl_gene_id', values = annotatedPeaktss$feature, mart = ensembl)\n"
	rscript += "annotrans <- cbind(annotatedPeaktrans, C1BM[match(annotatedPeaktrans$feature, C1BM[,1]), 2:4])\n"
	rscript += "annotss <- cbind(annotatedPeaktss, C2BM[match(annotatedPeaktss$feature, C2BM[,1]), 2:4])\n"
	rscript += "write.table(annotss, file='{}_nearest_tss.tsv', sep='\\t', quote=F, row.names=F)\n".format(output_prefix)
	rscript += "write.table(annotrans, file='{}_nearest_gene.tsv', sep='\\t', quote=F, row.names=F)\n".format(output_prefix)
	return rscript

def run_rcode(rscript, name):
	rcode = open(name, "w")
	rcode.write(rscript)
	rcode.close()
	try:
		subprocess.call(['Rscript', name])
	except:
		error("Error in running {}\n".format(name))
		error("Error: %s\n" % str(sys.exc_info()[1]))
		error( "[Exception type: %s, raised in %s:%d]\n" % ( sys.exc_info()[1].__class__.__name__, 
		os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
		traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
		sys.exit(1)
		
def main():
	parser = argparse.ArgumentParser(description='Annotation of Peaks.\n')
	parser.add_argument('-p', '--peak', help='Peak file', required=True)
	parser.add_argument('-o', '--output', help='Output prefix', required=True)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())
	rscript = chippeakanno(args["peak"], args["output"])
	run_rcode(rscript, "chipseqanno_rcode.R")
		