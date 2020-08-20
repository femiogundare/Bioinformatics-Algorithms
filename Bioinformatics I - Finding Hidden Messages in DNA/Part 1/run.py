# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 23:33:17 2020

@author: femiogundare
"""


import argparse
from task_1 import *
"""
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\dataset_2_7.txt', 'r') as file:
    lines = file.readlines()
    text, pattern = lines[0].strip(), lines[1].strip()
    PatternCount(text=text, pattern=pattern)

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\dataset_2_10.txt', 'r') as file:
    lines = file.readlines()
    text, k = lines[0].strip(), int(lines[1].strip())
    print(FrequentWords(text, k)(text=text, k=k))
    
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\dataset_3_2.txt', 'r') as file:
    text = file.read().strip()
    print(ReverseComplement(text))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\dataset_3_5.txt', 'r') as file:
    lines = file.readlines()
    pattern, text = lines[0].strip(), lines[1].strip()
    print(PatternMatch(text=text, pattern=pattern))
    
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\Vibrio_cholerae.txt', 'r') as file:
    text = file.read().strip()
    print(PatternMatch(text=text, pattern='CTTGATCAT'))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\dataset_2994_5.txt', 'r') as file:
    lines = file.readlines()
    text, k = lines[0].strip(), int(lines[1].strip())
    print(ComputingFrequencies(Text=text, k=k))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\dataset_3010_2.txt', 'r') as file:
    pattern = file.read().strip()
    print(PatternToNumber(Pattern=pattern))
    
with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\dataset_3010_5.txt', 'r') as file:
    lines = file.readlines()
    index, k = int(lines[0].strip()), int(lines[1].strip())
    print(NumberToPattern(index, k))

with open('C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics I - Finding Hidden Messages in DNA\Part 1\datasets\E_coli.txt', 'r') as file:
    text = file.read().strip()
    #print(ClumpFinding(Genome=text, k=k, t=t, L=L))
    print(BetterClumpFinding(Genome=text, k=9, t=3, L=500))
"""
print(PatternCount('ACTGTACGATGATGTGTGTCAAAG', 'TGT'))
print(PatternMatch('GACGATATACGACGATA', 'ATA'))
print(FasterFrequentWords('CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA', 3))