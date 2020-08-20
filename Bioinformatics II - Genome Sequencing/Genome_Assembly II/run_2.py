# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 21:56:08 2020

@author: femiogundare
"""


import argparse
from codes_2 import *

OUTPUT_PATH = 'C:\\Users\\Dell\\Desktop\\Bioinformatics\\Bioinformatics II - Genome Sequencing\\Genome_Assembly II\\output'
DATASET_PATH = 'C:\\Users\\Dell\\Desktop\\Bioinformatics\\Bioinformatics II - Genome Sequencing\\Genome_Assembly II\\datasets'



""""
data = DATASET_PATH + '/dataset_6206_4 (1).txt'
with open(data, 'r') as file:
    lines = file.readlines()
    parameters, gapped_patterns = lines[0].strip().split(), lines[1:]
    k, d = parameters[0], parameters[1]
    
    gapped_patterns = [gapped_pattern.strip() for gapped_pattern in gapped_patterns]
    gapped_patterns = [gapped_pattern.split('|') for gapped_pattern in gapped_patterns]
    
    string_spelled = StringSpelledByGappedPatterns(gapped_patterns=gapped_patterns, k=int(k), d=int(d))
    
    output_path = OUTPUT_PATH + '\string_spelled_by_gapped_patterns.txt'
    f = open(output_path, 'w')
    f.write(string_spelled)
    f.close()
"""

data = DATASET_PATH + '/dataset_204_16.txt'
with open(data, 'r') as file:
    lines = file.readlines()
    parameters, read_pairs = lines[0].strip().split(), lines[1:]
    k, d = parameters[0], parameters[1]
    
    read_pairs = [read_pair.strip() for read_pair in read_pairs]
    read_pairs = [read_pair.split('|') for read_pair in read_pairs]
    
    stringReconstruction = StringReconstructionFromReadPairs(k=int(k), d=int(d), read_pairs=read_pairs)
    
    output_path = OUTPUT_PATH + '\string_reconstruction_from_read_pairs.txt'
    f = open(output_path, 'w')
    f.write(stringReconstruction)
    f.close()