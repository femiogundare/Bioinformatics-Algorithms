# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 13:20:56 2020

@author: femiogundare
"""

import argparse
from codes import *

OUTPUT_PATH = 'C:\\Users\\Dell\\Desktop\\Bioinformatics\\Bioinformatics II - Genome Sequencing\\Genome_Assembly I\\output'
DATASET_PATH = 'C:\\Users\\Dell\\Desktop\\Bioinformatics\\Bioinformatics II - Genome Sequencing\\Genome_Assembly I\\datasets'


"""
data = DATASET_PATH + '/dataset_197_3.txt'
with open(data, 'r') as file:
    lines = file.readlines()
    k, Text = lines[0].strip(), lines[1].strip()
    string_composition = StringComposition(Text, int(k))
    
    string_composition = '\n'.join(string_composition)
    output_path = OUTPUT_PATH + '\string_composition.txt'
    f = open(output_path, 'w')
    f.write(string_composition)
    f.close()


data = DATASET_PATH + '/dataset_198_3.txt'
with open(data, 'r') as file:
    paths = file.readlines()
    paths = [path.strip() for path in paths]
    genome = PathToGenome(paths)
    
    output_path = OUTPUT_PATH + '\path_to_genome.txt'
    f = open(output_path, 'w')
    f.write(genome)
    f.close()


data = DATASET_PATH + '/dataset_198_10.txt'
with open(data, 'r') as file:
    patterns = file.readlines()
    patterns = [pattern.strip() for pattern in patterns]
    overlap_graph = OverlapGraph(Patterns=patterns)
    
    output_path = OUTPUT_PATH + '\overlap_graph.txt'
    f = open(output_path, 'w')
    
    for key, val in overlap_graph.items():
        if val != []:
            values = ','.join(val)
            f.writelines(key + ' -> ' + values)
            f.writelines('\n')
    f.close()
    

data = DATASET_PATH + '/dataset_199_6.txt'
with open(data, 'r') as file:
    lines = file.readlines()
    k, Text = lines[0].strip(), lines[1].strip()
    debruijnGraph = DeBruijnGraph(Text=Text, k=int(k))
    
    output_path = OUTPUT_PATH + '\DeBruijnGraph_graph.txt'
    f = open(output_path, 'w')
    
    for key, val in debruijnGraph.items():
        if val != []:
            values = ','.join(val)
            f.writelines(key + ' -> ' + values)
            f.writelines('\n')
    f.close()
"""

data = DATASET_PATH + '/dataset_200_8.txt'
with open(data, 'r') as file:
    kmers = file.readlines()
    kmers = [kmer.strip() for kmer in kmers]
    debruijnGraph = DeBruijnGraphFromKmers(kmers=kmers)
    
    output_path = OUTPUT_PATH + '\DeBruijnGraph_fromKmers.txt'
    f = open(output_path, 'w')
    
    for key, val in debruijnGraph.items():
        if val != []:
            values = ','.join(val)
            f.writelines(key + ' -> ' + values)
            f.writelines('\n')
    f.close()