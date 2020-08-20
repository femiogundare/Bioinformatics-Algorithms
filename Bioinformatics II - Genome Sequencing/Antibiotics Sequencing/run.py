# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 13:11:02 2020

@author: femiogundare
"""

import argparse
from codes import *

OUTPUT_PATH = 'C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics II - Genome Sequencing\\Antibiotics Sequencing\\output'
DATASET_PATH = 'C:\\Users\\Dell\\Desktop\Bioinformatics\\Bioinformatics II - Genome Sequencing\\Antibiotics Sequencing\\datasets'


"""
data = DATASET_PATH + '/dataset_96_4.txt' 
with open(data, 'r') as file:
    
    Rna_string = file.read().strip('\n')
    
    protein_translation = ProteinTranslationProblem(Rna_string=Rna_string)
    
    output_path = OUTPUT_PATH + '\Protein_translation_problem.txt'
    f = open(output_path , 'w')
    f.write(protein_translation)
    f.close()


data = DATASET_PATH + '/dataset_96_7.txt' 
with open(data, 'r') as file:
    lines = file.readlines()
    text, peptide = lines[0].strip(), lines[1].strip()
    peptide_encoding = PeptideEncodingProblem(Text=text, Peptide=peptide)
    peptide_encoding = '\n'.join(peptide_encoding)
    
    output_path = OUTPUT_PATH + '\Peptide_encoding_problem.txt'
    f = open(output_path , 'w')
    f.write(peptide_encoding)
    f.close()


data = DATASET_PATH + '/dataset_98_4.txt' 
with open(data, 'r') as file:
    peptide = file.read().strip('\n')
    theoretical_spectrum = TheoreticalSpectrumProblem(Peptide=peptide)
    theoretical_spectrum = ' '.join([str(x) for x in theoretical_spectrum])
    
    output_path = OUTPUT_PATH + '/theoretical_spectrum_problem.txt'
    f = open(output_path , 'w')
    f.write(theoretical_spectrum)
    f.close()


data = DATASET_PATH + '/dataset_100_6.txt' 
with open(data, 'r') as file:
    spectrum = file.read().strip('\n')
    spectrum = [int(mass) for mass in spectrum.split(' ')]
    
    cyclopeptide_sequencing = CyclopeptideSequencing(Spectrum=spectrum)
    cyclopeptide_sequencing = ' '.join([str(x) for x in cyclopeptide_sequencing])
    
    output_path = OUTPUT_PATH + '/cyclopeptide_sequencing.txt'
    f = open(output_path , 'w')
    f.write(cyclopeptide_sequencing)
    f.close()


data = DATASET_PATH + '/dataset_100_6 (1).txt' 
with open(data, 'r') as file:
    spectrum = file.read().strip('\n')
    spectrum = [int(mass) for mass in spectrum.split(' ')]
    
    cyclopeptide_sequencing = CyclopeptideSequencing(Spectrum=spectrum)
    cyclopeptide_sequencing = ' '.join([str(x) for x in cyclopeptide_sequencing])
    
    output_path = OUTPUT_PATH + '/cyclopeptide_sequencing_1.txt'
    f = open(output_path , 'w')
    f.write(cyclopeptide_sequencing)
    f.close()


data = DATASET_PATH + '/dataset_102_3.txt' 
with open(data, 'r') as file:
    lines = file.readlines()
    peptide, spectrum = lines[0].strip(), lines[1]
    spectrum = spectrum.strip('\n')
    spectrum = [int(mass) for mass in spectrum.split(' ')]
    
    cyclopeptide_scoring = CyclopeptideScoring(Peptide=peptide, Spectrum=spectrum)
    
    output_path = OUTPUT_PATH + '/cyclopeptide_scoring.txt'
    f = open(output_path , 'w')
    f.write(str(cyclopeptide_scoring))
    f.close()


data = DATASET_PATH + '/dataset_4913_1.txt' 
with open(data, 'r') as file:
    lines = file.readlines()
    peptide, spectrum = lines[0].strip(), lines[1]
    spectrum = spectrum.strip('\n')
    spectrum = [int(mass) for mass in spectrum.split(' ')]
    
    linear_peptide_scoring = LinearPeptideScoring(Peptide=peptide, Spectrum=spectrum)
    
    output_path = OUTPUT_PATH + '/linear_peptide_scoring.txt'
    f = open(output_path , 'w')
    f.write(str(linear_peptide_scoring))
    f.close()


data = DATASET_PATH + '/dataset_4913_3.txt' 
with open(data, 'r') as file:
    lines = file.readlines()
    Leaderboard, Spectrum, N = lines[0].strip(), lines[1:-1], lines[-1].strip()
    Leaderboard = [peptide for peptide in Leaderboard.split(' ')]
    
    Spectrum = Spectrum[0].strip('\n')
    Spectrum = [int(spectrum) for spectrum in Spectrum.split(' ')]
    
    trim = Trim(Leaderboard, Spectrum, int(N))
    trim = ' '.join([str(x) for x in trim])
    
    output_path = OUTPUT_PATH + '/trim.txt'
    f = open(output_path , 'w')
    f.write(str(trim))
    f.close()


data = DATASET_PATH + '/dataset_102_8.txt' 
with open(data, 'r') as file:
    lines = file.readlines()
    N, Spectrum = lines[0].strip(), lines[1].strip('\n')
    
    Spectrum = [int(spectrum) for spectrum in Spectrum.split(' ')]
    
    leaderboard_cyclopeptide_sequencing = LeaderboardCyclopeptideSequencing(Spectrum, int(N))
    
    output_path = OUTPUT_PATH + '/leaderboard_cyclopeptide_sequencing.txt'
    f = open(output_path , 'w')
    f.write(str(leaderboard_cyclopeptide_sequencing))
    f.close()


#LEADERBOARD CYCLOPEPTIDE SEQUENCING FOR TYROCIDINE B1 25% MISSING/FALSE
data = DATASET_PATH + '/Tyrocidine_B1_Spectrum_25.txt' 
with open(data, 'r') as file:
    Spectrum = file.readlines()[0].strip()
    print(Spectrum)
    N = 1000
    
    Spectrum = [int(spectrum) for spectrum in Spectrum.split(' ')]
    
    leaderboard_cyclopeptide_sequencing = LeaderboardCyclopeptideSequencing(Spectrum, int(N))
    
    output_path = OUTPUT_PATH + '/leaderboard_cyclopeptide_sequencing_tyrocidine_b1_25_percent.txt'
    f = open(output_path , 'w')
    f.write(str(leaderboard_cyclopeptide_sequencing))
    f.close()


data = DATASET_PATH + '/dataset_104_4.txt' 
with open(data, 'r') as file:
    Spectrum = file.read().strip('\n') 
    
    Spectrum = [int(spectrum) for spectrum in Spectrum.split(' ')]
    
    spectral_convolution = SpectralConvolution(Spectrum)
    spectral_convolution = ' '.join([str(x) for x in spectral_convolution])
    
    output_path = OUTPUT_PATH + '/spectral_convolution.txt'
    f = open(output_path , 'w')
    f.write(str(spectral_convolution))
    f.close()
"""

data = DATASET_PATH + '/dataset_104_7.txt' 
with open(data, 'r') as file:
    lines = file.readlines()
    M, N, Spectrum = lines[0].strip(), lines[1].strip(), lines[2:]
    
    Spectrum = Spectrum[0].strip('\n')
    Spectrum = sorted([int(spectrum) for spectrum in Spectrum.split(' ')])
    
    convolutional_cyclopeptide = ConvolutionalCyclopeptideSequencing(int(M), int(N), Spectrum)
    #convolutional_cyclopeptide = ' '.join([str(x) for x in convolutional_cyclopeptide])
    
    output_path = OUTPUT_PATH + '/convolutional_cyclopeptide_sequencing.txt'
    f = open(output_path , 'w')
    f.write(str(convolutional_cyclopeptide))
    f.close()