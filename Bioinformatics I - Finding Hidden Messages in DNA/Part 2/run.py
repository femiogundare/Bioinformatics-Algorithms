# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 21:03:30 2020

@author: femiogundare
"""

from task_2 import *
"""
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 2\datasets\dataset_7_6.txt', 'r') as file:
    text = file.read().strip()
    print(MinimumSkew(Genome=text))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 2\datasets\dataset_9_3.txt', 'r') as file:
    pattern_1, pattern_2 = file.readlines()
    print(HammingDistance(p=pattern_1, q=pattern_2))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 2\datasets\dataset_9_4.txt', 'r') as file:
    lines = file.readlines()
    pattern, text, d = lines[0].strip(), lines[1].strip(), int(lines[2].strip())
    print(ApproximatePatternMatching(Text=text, Pattern=pattern, d=d))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 2\datasets\dataset_9_6.txt', 'r') as file:
    lines = file.readlines()
    pattern, text, d = lines[0].strip(), lines[1].strip(), int(lines[2].strip())
    print(ApproximatePatternCount(Text=text, Pattern=pattern, d=d))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 2\datasets\dataset_3014_4.txt', 'r') as file:
    text, d = file.readlines()
    print(GenerateMismatchedKmers(kmer=text, d=int(d)))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 2\datasets\dataset_9_7.txt', 'r') as file:
    lines = file.readlines()
    text, parameters =  lines[0].strip(), lines[1].strip().split()
    k, d = int(parameters[0]), int(parameters[1])
    print(FrequentWordsWithMismatches(Text=text, k=k, d=d))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 2\datasets\dataset_9_8.txt', 'r') as file:
    lines = file.readlines()
    text, parameters =  lines[0].strip(), lines[1].strip().split()
    k, d = int(parameters[0]), int(parameters[1])
    print(FrequentWordsWithMismatchesAndReverseComplement(Text=text, k=k, d=d))
"""

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 2\datasets\dataset_3014_4.txt', 'r') as file:
    text, d = file.readlines()
    print(GenerateMismatchedKmers(kmer=text, d=int(d)))