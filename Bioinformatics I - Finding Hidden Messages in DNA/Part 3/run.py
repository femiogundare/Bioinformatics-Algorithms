# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 23:40:34 2020

@author: femiogundare
"""

from task_3 import *

"""
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_156_8.txt', 'r') as file:
    lines = file.readlines()
    parameters, dna_strings = lines[0].strip().split(), lines[1:]
    dna_strings = [dna.strip() for dna in dna_strings]
    k, d = parameters[0], parameters[1]
    print(MotifEnumeration(Dna_list=dna_strings, k=int(k), d=int(d)))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_156_8 (1).txt', 'r') as file:
    lines = file.readlines()
    parameters, dna_strings = lines[0].strip().split(), lines[1:]
    dna_strings = [dna.strip() for dna in dna_strings]
    k, d = parameters[0], parameters[1]
    print(MotifEnumeration(Dna_list=dna_strings, k=int(k), d=int(d)))
    
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_158_9.txt', 'r') as file:
    lines = file.readlines()
    k, text = lines[0].strip(), lines[1:]
    text = [t.strip() for t in text]
    print(MedianString(Dna_list=text, k=int(k)))
    
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_5164_1.txt', 'r') as file:
    lines = file.readlines()
    pattern, text = lines[0].strip(), lines[1]
    text = [t.strip() for t in text.split()]
    print(DistanceBetweenPatternAndStrings(Pattern=pattern, Dna_list=text))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_158_9.txt', 'r') as file:
    lines = file.readlines()
    k, text = lines[0].strip(), lines[1:]
    text = [t.strip() for t in text]
    print(EfficientMedianString(Dna_list=text, k=int(k)))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_159_3.txt', 'r') as file:
    lines = file.readlines()
    text, k, matrix = lines[0].strip(), lines[1].strip(), lines[2:]
    print(ProfileMostProbableKmer(Text=text, k=int(k), profile_matrix=matrix))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_159_5.txt', 'r') as file:
    lines = file.readlines()
    parameters, dna_strings = lines[0].strip().split(), lines[1:]
    dna_strings = [dna.strip() for dna in dna_strings]
    k, t = parameters[0], parameters[1]
    greedy_motifs = GreedyMotifSearch(k=int(k), t=int(t), Dna=dna_strings)
    with open('C:\\Users\\Dell\\Desktop\\sample.txt', 'w') as myfile:
        myfile.write('\n'.join(greedy_motifs))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_160_9.txt', 'r') as file:
    lines = file.readlines()
    parameters, dna_strings = lines[0].strip().split(), lines[1:]
    dna_strings = [dna.strip() for dna in dna_strings]
    k, t = parameters[0], parameters[1]
    greedy_motifs_with_pseudocount = GreedyMotifSearchWithPseudoCount(k=int(k), t=int(t), Dna=dna_strings)
    with open('C:\\Users\\Dell\\Desktop\\sample.txt', 'w') as myfile:
        myfile.write('\n'.join(greedy_motifs_with_pseudocount))
"""
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 3\datasets\dataset_161_5.txt', 'r') as file:
    lines = file.readlines()
    parameters, dna_strings = lines[0].strip().split(), lines[1:]
    dna_strings = [dna.strip() for dna in dna_strings]
    k, t = parameters[0], parameters[1]
    last_motifs = RandomizedMotifSearch(Dna_list=dna_strings, k=int(k), t=int(t))
    i = 0
    while i<1000:
        best_motifs = RandomizedMotifSearch(Dna_list=dna_strings, k=int(k), t=int(t))
    
        if Score(best_motifs) < Score(last_motifs):
            last_motifs = best_motifs
        
        i +=1
    with open('C:\\Users\\Dell\\Desktop\\sample.txt', 'w') as myfile:
        myfile.write('\n'.join(last_motifs))