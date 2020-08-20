# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 23:11:52 2020

@author: Dell
"""

import argparse
from codes import *

OUTPUT_PATH = 'C:\\Users\\Dell\\Desktop\\Bioinformatics\\Bioinformatics II - Genome Sequencing\\Genome_Assembly II\\output'
DATASET_PATH = 'C:\\Users\\Dell\\Desktop\\Bioinformatics\\Bioinformatics II - Genome Sequencing\\Genome_Assembly II\\datasets'


"""
data = DATASET_PATH + '/dataset_203_2.txt' 
with open(data, 'r') as file:
    graph = dict((line.strip().split(' -> ') for line in file))
    for key in graph:
        graph[key] = graph[key].split(',')
        
    path_list = EulerianCircuit(graph)
    
    path_display = "->".join(str(elements) for elements in path_list)
    output_path = OUTPUT_PATH + '\Eulerian_Cycle.txt'
    f = open(output_path , 'w')
    f.write(path_display)
    f.close()


data = DATASET_PATH + '/dataset_203_6.txt' 
with open(data, 'r') as file:
    graph = dict((line.strip().split(' -> ') for line in file))
    for key in graph:
        graph[key] = graph[key].split(',')
        
    path_list = EulerianPath(directed_graph=graph)
    
    path_display = "->".join(str(elements) for elements in path_list)
    output_path = OUTPUT_PATH + '\Eulerian_Path.txt'
    f = open(output_path , 'w')
    f.write(path_display)
    f.close()


data = DATASET_PATH + '/dataset_203_7.txt'
with open(data, 'r') as file:
    lines = file.readlines()
    k, Patterns = lines[0].strip(), lines[1:]
    Patterns = [pattern.strip() for pattern in Patterns]
    
    string_reconstruction = StringReconstruction(patterns=Patterns)
    
    string_reconstruction = ''.join(string_reconstruction)
    output_path = OUTPUT_PATH + '\string_reconstruction.txt'
    f = open(output_path, 'w')
    f.write(string_reconstruction)
    f.close()


data = DATASET_PATH + '/dataset_203_11.txt'
with open(data, 'r') as file:
    lines = file.readlines()
    k = lines[0].strip()
    print(k)
    
    k_universal_circular_string = k_UniversalCircularString(k=int(k))
    
    #k_universal_circular_string = ''.join(k_universal_circular_string)
    output_path = OUTPUT_PATH + '\k-Universal_Circular_String.txt'
    f = open(output_path, 'w')
    f.write(k_universal_circular_string)
    f.close()

    
"""