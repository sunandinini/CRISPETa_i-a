#!/usr/bin/env python3
import re, sys, os
from Bio.Seq import Seq
import subprocess

# Doench Score Rule_set2.0
script_dir = sys.argv[3]

RS2_path = script_dir + '/../Rule_Set_2_scoring_v1/analysis'

sys.path.insert(0,RS2_path)
import pandas as pd
import csv, argparse, sys
import pickle
import model_comparison
import rs2_score_calculator

global model

# read input and output files from command line
input_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[2], 'w')

# Functions
#def get_genome_coordinates(header):
#	h = header.split(":")[-1]
#	chrom = header.split(":")[2]
#	strand = h.split("(")[-1][0]
#	start = h.split("(")[0].split("-")[0]
#	end = h.split("(")[0].split("-")[1]
#
#	return chrom, start, end, strand

def load_model(path=script_dir):
    model_file = path + '/../Rule_Set_2_scoring_v1/saved_models/V3_model_nopos.pickle'
    try:
        with open(model_file, 'rb') as f:
            model= pickle.load(f)
            return model  
    except:
        raise Exception("could not find model stored to file %s" % model_file)

#Load model for scoring sgRNAs
model = load_model()

def calc_score_rule_set_2(seq, strand, model=model):
	if strand == "-":
		seq = str(Seq(seq).reverse_complement())

	if 'N' in seq:
		score = 0
	else:
		score = model_comparison.predict(seq, -1, -1, model)
	return(score)

# Print header in output
output_file.write("\t".join(["target_id", "chr", "sgRNA_start", "sgRNA_end", "sgRNA_context", "sgRNA_strand", "off0", "off1", "off2", "off3", "off4", "score_rule_set2", "\n"]))

# Scoring sgRNAs
for line in input_file:

	if line.startswith("target_id"): continue

	line = line.strip().split("\t")
	sgRNA = line[4]
	strand = line[5]

	score = calc_score_rule_set_2(sgRNA, strand)
	output_file.write("\t".join(line+[str(score)]+["\n"]))

input_file.close()
output_file.close()