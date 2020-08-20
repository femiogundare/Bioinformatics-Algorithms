# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 08:47:15 2020

@author: femiogundare
"""

import math
import random
import numpy as np


def NumberToPattern(index, k):
    """
    Sample Input: 45; 4
    Sample Output: AGTC
    """
    number_to_symbol = {0:'A', 1:'C', 2:'G', 3:'T'}
    
    if k==1:
        return number_to_symbol[index]
    
    prefixIndex = index//4
    remainder = index%4
    symbol = number_to_symbol[remainder]
    PrefixPattern = NumberToPattern(prefixIndex, k-1)
    return ''.join([PrefixPattern, symbol])


def HammingDistance(p, q):
    """
    Computes the hamming distance between two k-mers p and q
    
    Sample Input: GGGCCGTTGGT; GGACCGTTGAC
    Sample Output: 3
    """
    assert len(p)==len(q)
    k = len(p)
    num_of_mismatches = 0
    
    for i in range(k):
        if p[i] != q[i]:
            num_of_mismatches +=1
        else:
            continue
    hamming_distance = num_of_mismatches
    
    return hamming_distance


def GenerateMismatchedKmers(kmer, d):
    mismatched_list = []
    nucleotides = ['A', 'C', 'G', 'T']
    if d == 0:
        return [kmer]
    if len(kmer) == 1:
        return ["A", "C", "G", "T"]
    
    # generate neighbor mismatches
    for neighbor in GenerateMismatchedKmers(kmer[1:], d):
        if HammingDistance(kmer[1:], neighbor) < d:
            for nucleotide in nucleotides:
                mismatched_list.append(nucleotide+neighbor)
        else:
            mismatched_list.append(kmer[0]+neighbor)
    
    return mismatched_list


def MotifEnumeration(Dna_list, k, d):
    """
    Parameters: 
        Dna_list is a list of Dna strings
        k and d are the number of nucleotides and number of mismatches respectively
    
    Sample Input:
        3 1
        ATTTGGC; TGCCTTA; CGGTATC; GAAAATT
        
    Sample Output: ATA ATT GTT TTT
    
    The function outputs the motifs  present in the DNA strings found in the Dna list
    *** A k-mer is a (k, d)-motif if it appears in every Dna string in the Dna list with at most d mismatches
    
    The easiest way to solve this problem is:
        1. Get all the kmers in the first DNA string that is in the DNA list
        2. For each of the kmers generated in (1) above, generate its possible neighbors with at most d mismatches;
           then store the all the neighbors of all the k-mers in ONE list
        3. Repeat (1) & (2) for the second, third, fourth, ..., and last DNA strings
        4. Then append each of the list of neighbors from k-mers of the first, second, third, ..., and last DNA strings
           to a list such that the list is like:
               all_neighbors = [
                   [all_neighbors_of_first_dna_string],
                   [all_neighbors_of_second_dna_string],
                   [all_neighbors_of_third_dna_string],
                   [all_neighbors_of_fourth_dna_string],
                                  .
                                  .
                                  .
                   [all_neighbors_of_last_dna_string]
                   ]
        5. Select neighbors that are found in each of the lists and call them 'common neighbors' — these neighbors are
           the motifs BECAUSE they satisfy the conditions of appearing in every DNA string with at most d mismatches!
    """
    
    #list-in-list, where each of the lists in this list contains the all the kmers in a DNA string
    kmers_in_each_dna_string = []     
    #list-in-list, where each of the lists in this list contains the all the possible neighbors of each kmer in a DNA string
    neighbors_in_each_dna_string = []
    
    for dna_string in Dna_list:
        kmers_in_this_dna_string = []
        neighbors_in_this_dna_string = []
        
        for i in range(len(dna_string)-k+1):
            kmer = dna_string[i : i+k]
            kmers_in_this_dna_string.append(kmer)   #step 1
            
            possible_neighbors_of_kmer = GenerateMismatchedKmers(kmer=kmer, d=d)
            for each_neighbor in possible_neighbors_of_kmer:
                neighbors_in_this_dna_string.append(each_neighbor)   #step 2
            
        kmers_in_each_dna_string.append(kmers_in_this_dna_string)
        neighbors_in_each_dna_string.append(neighbors_in_this_dna_string)   #step 4
        
    all_neighbors = neighbors_in_each_dna_string
    common_neighbors = set(all_neighbors[0]).intersection(*all_neighbors)   #step 5
    common_neighbors = ' '.join(common_neighbors)
    
    return common_neighbors



def MedianString(Dna_list, k):
    """
    Input: A collection of strings Dna and an integer k.
    Output: A k-mer Pattern minimizing d(Pattern, Dna) among all k-mers Pattern in all strings in Dna list
    
    Sample Input: 3; AAATTGACGCAT, GACGACCACGTT, CGTCAGCGCCTG, GCTGAGCACCGG, AGTTCGGGACAG
    Sample Output: GAC
    
    The goal of the Median string problem is that in a Dna list with several dna strings, we should find a PATTERN
    whose hamming distance between any k-mer in the strings in the Dna list is minimum.
    The best way to solve this is:
        1. Obtain all the k-mers in the first Dna string and append them to a list
        2. For each k-mer in the first Dna string, compute its hamming distance with other k-mers in the first string
           and select the k-mer whose hamming distance is minimum with the other k-mer
        3. Repeat (1) and (2) for other Dna strings, but using the k-mers of the previous strings and the current 
           Dna string to find the hamming distance with the k-mers of the current string
    """
    
    distance = float('inf')
    patterns = []
    
    for dna_string in Dna_list:
        for i in range(len(dna_string)-k+1):
            pattern = dna_string[i : i+k]
            if pattern not in patterns:
                patterns.append(pattern)
        
        for pattern in patterns:
            for j in range(len(dna_string)-k+1):
                hamming_distance = HammingDistance(pattern, dna_string[j : j+k])
                if distance > hamming_distance:
                    distance = hamming_distance
                    Median = pattern
                
    return Median



def DistanceBetweenPatternAndStrings(Pattern, Dna_list):
    """
    Input: A string Pattern followed by a collection of strings Dna.
    Output: d(Pattern, Dna).
    
    Sample Input: AAA; TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT
    Sample Output: 5
    
    This function calculates the sum of hamming distance between Pattern and the k-mer with the minimum hamming distance
    in each of the strings
    """
    k = len(Pattern)
    distance = 0
    
    for dna_string in Dna_list:
        hamming_distance = float('inf')
        for i in range(len(dna_string)-k+1):
            pattern_prime = dna_string[i : i+k]
            hammingDistance = HammingDistance(Pattern, pattern_prime) 
            if hamming_distance > hammingDistance:
                hamming_distance = hammingDistance
        
        distance = distance + hamming_distance
    
    return distance



def EfficientMedianString(Dna_list, k):
    distance = float('inf')
    for i in range(4**k):
        Pattern = NumberToPattern(index=i, k=k)
        
        if distance > DistanceBetweenPatternAndStrings(Pattern=Pattern, Dna_list=Dna_list):
            distance = DistanceBetweenPatternAndStrings(Pattern=Pattern, Dna_list=Dna_list)
            Median = Pattern
    
    return Median



def ProfileMostProbableKmer(Text, k, profile_matrix):
    """
    Input: A string Text, an integer k, and a 4 × k matrix Profile.
    Output: A Profile-most probable k-mer in Text.
    
    Given a Dna string, k and a profile matrix, find the Profile-most probable k-mer i.e the k-mer with the 
    highest probability in Text
    """
    
    kmer_and_probability = {}
    #profile_matrix = {key: list(map(float, profile_matrix[num].split(' '))) for (num,key) in enumerate('ACGT')}
    
    for i in range(len(Text)-k+1):
        kmer = Text[i : i+k]
        #A_row, C_row, G_row, T_row = profile_matrix['A'], profile_matrix['C'], profile_matrix['G'], profile_matrix['T']
        A_row, C_row, G_row, T_row = profile_matrix[0], profile_matrix[1], profile_matrix[2], profile_matrix[3]
        probability_of_kmer = []
        
        for j in range(len(kmer)):
            if kmer[j] == 'A':
                probability_of_kmer.append(A_row[j])
            elif kmer[j] == 'C':
                probability_of_kmer.append(C_row[j])
            elif kmer[j] == 'G':
                probability_of_kmer.append(G_row[j])
            elif kmer[j] == 'T':
                probability_of_kmer.append(T_row[j])
        probability_of_kmer = np.prod(np.array(probability_of_kmer))
        kmer_and_probability[kmer] = probability_of_kmer
        
    maxProbability = max([prob for prob in kmer_and_probability.values()])
    profile_most_probable_kmer = [kmer for kmer, prob in kmer_and_probability.items() if prob==maxProbability]
    
    return ''.join(profile_most_probable_kmer[0])



def Profile(passed_motifs, k):
    "Creates a profile matrix for each passed profile matrix."
    # create empty matrix of floats
    matrix = []
    for i in range(4):
        # make k number of [0.0] entries in matrix for each spot in kmer
        matrix.append([0.0] * k)
    #print('K: ',k)
    #print(matrix)    
    # for each position in kmer, count bases
    number_motifs = len(passed_motifs)
    #print('number of motifs: ', number_motifs)
    for i in range(k):
        motif_count = {"A" : 0, "C" : 0, "G" : 0, "T" : 0}
        for motif in passed_motifs:
            #print('motif length: ', len(motif))
            motif_count["A"] += motif[i].count("A")
            motif_count["C"] += motif[i].count("C")
            motif_count["G"] += motif[i].count("G")
            motif_count["T"] += motif[i].count("T")

        # create matrix of profiles for each base
        matrix[0][i] = motif_count["A"] / number_motifs
        matrix[1][i] = motif_count["C"] / number_motifs
        matrix[2][i] = motif_count["G"] / number_motifs
        matrix[3][i] = motif_count["T"] / number_motifs
        
    return matrix

def profile(motifs):
    '''Returns the profile of the dna list motifs.'''
    columns = [''.join(seq) for seq in zip(*motifs)]
    return [[float(col.count(nuc)) / float(len(col)) for nuc in 'ACGT'] for col in columns]

def Score(motifs):
    """"
    Computes the score of the motifs by comparing the nucleotides of each of the motifs with the consensus string
    """
    k = len(motifs[0])
    #get the consensus string from the profile matrix
    consensus = []
    for i in range(k):
        freq = {"A" : 0, "C" : 0, "G" : 0, "T" : 0}
        for motif in motifs:
            freq["A"] += motif[i].count("A")
            freq["C"] += motif[i].count("C")
            freq["G"] += motif[i].count("G")
            freq["T"] += motif[i].count("T")
    
        # based on freq above, creat a consensus kmer to compare to passed motif
        max_freq = max(freq.values())
        for nt, count in freq.items():
            if count == max_freq:
                consensus.append(nt)
                break
                
    consensus = "".join(consensus)
    score_of_each_motif = []
    for motif in motifs:
        score_of_motif = HammingDistance(consensus, motif)
        score_of_each_motif.append(score_of_motif)
        
    score = sum(score_of_each_motif)
    
    return score



def GreedyMotifSearch(k, t, Dna):
    # create list of best motifs from out of the first DNA string given
    best_motifs = []
    for seq in Dna:
        first = seq[0:k]
        best_motifs.append(first)
    # best_motifs = [seq[:k] for seq in dna]
    
    # iterate over kmers in first Dna string, create a motif list for each kmer
    first_seq = Dna[0]
    for start in range(len(first_seq) - k + 1):
        kmer = first_seq[start : start + k]
        # start motif list based on this kmer from first_seq
        motif = [kmer]
        
        # iterate over subsequent Dna strings, make  profile from them based on first_seq kmer
        for i in range(1, t):
            matrix = Profile(motif, k)
            most_probable = ProfileMostProbableKmer(Dna[i], k, matrix)
            motif.append(most_probable)
            
        # score motif, replace if best
        if Score(motif) < Score(best_motifs):
            best_motifs = motif
            
    return best_motifs



def PseudocountProfile(passed_motifs, k):
    "Creates a profile matrix for each passed profile matrix."
    # create empty matrix of floats
    matrix = []
    for i in range(4):
        # make k number of [0.0] entries in matrix for each spot in kmer
        matrix.append([0.0] * k)
        
    #print('K: ',k)
    # for each position in kmer, count bases
    total_counts = len(passed_motifs) + 4
    for i in range(k):
        motif_count = {"A" : 1, "C" : 1, "G" : 1, "T" : 1}
        for motif in passed_motifs:
            #print('motif length: ', len(motif))
            motif_count["A"] += motif[i].count("A")
            motif_count["C"] += motif[i].count("C")
            motif_count["G"] += motif[i].count("G")
            motif_count["T"] += motif[i].count("T")

        # create matrix of profiles for each base
        matrix[0][i] = motif_count["A"] / total_counts
        matrix[1][i] = motif_count["C"] / total_counts
        matrix[2][i] = motif_count["G"] / total_counts
        matrix[3][i] = motif_count["T"] / total_counts
        
    return matrix


def GreedyMotifSearchWithPseudoCount(k, t, Dna):
    # create list of best motifs from out of the first DNA string given
    best_motifs = []
    for seq in Dna:
        first = seq[0:k]
        best_motifs.append(first)
    # best_motifs = [seq[:k] for seq in dna]
    
    # iterate over kmers in first Dna string, create a motif list for each kmer
    first_seq = Dna[0]
    for start in range(len(first_seq) - k + 1):
        kmer = first_seq[start : start + k]
        # start motif list based on this kmer from first_seq
        motif = [kmer]
        
        # iterate over subsequent Dna strings, make  profile from them based on first_seq kmer
        for i in range(1, t):
            matrix = PseudocountProfile(motif, k)
            most_probable = ProfileMostProbableKmer(Dna[i], k, matrix)
            motif.append(most_probable)
            
        # score motif, replace if best
        if Score(motif) < Score(best_motifs):
            best_motifs = motif
            
    return best_motifs



def RandomizedMotifSearch(Dna_list, k, t):
    #get a random kmer of length k from each string in Dna
    best_motifs = []
    
    for dna_string in Dna_list:
        overlap = len(dna_string)-k+1
        i = random.randint(0, overlap-1)
        kmer = dna_string[i : i+k]
        best_motifs.append(kmer)
        
    while True:
        profile_matrix = PseudocountProfile(passed_motifs=best_motifs, k=k)
        motifs = []
        
        for dna_string in Dna_list:
            most_probable_kmer = ProfileMostProbableKmer(dna_string, k, profile_matrix)
            motifs.append(most_probable_kmer)
        #print(motifs)  
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs[:]
        else:
            return best_motifs



def GibbsSampler(Dna_list, k, t, N):
    #get a random kmer of length k from each string in Dna
    best_motifs = []
    
    for dna_string in Dna_list:
        overlap = len(dna_string)-k+1
        i = random.randint(0, overlap-1)
        kmer = dna_string[i : i+k]
        best_motifs.append(kmer)
        
    motifs = best_motifs 
    
    for j in range(0, N):
        #select a random Dna string
        i = random.choice(range(0, t))
        
        #make a profile matrix for all strings t except i
        motifs_subset = motifs[:]
        #remove the position i
        motifs_subset.pop(i)
        # created a profile matrix based on subset
        matrix = PseudocountProfile(motifs_subset, k)
        # profile randomly generated a motif for i
        most_probable = ProfileMostProbableKmer(Dna_list[i], k, matrix)
        # put new kmer in place on original list
        motifs[i] = most_probable 
        
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
            
    return best_motifs


dna = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']



def RandomizedMotifSearch_k(Dna_list, k, t):
    #get a random kmer of length k from each string in Dna
    best_motifs = ['CCA',

'CCT',

'CTT',

'TTG']
    
        
    while True:
        profile_matrix = PseudocountProfile(passed_motifs=best_motifs, k=k)
        motifs = []
        
        for dna_string in Dna_list:
            most_probable_kmer = ProfileMostProbableKmer(dna_string, k, profile_matrix)
            motifs.append(most_probable_kmer)
        #print(motifs)  
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs[:]
        else:
            return best_motifs
        
dna_strings = ['AAGCCAAA',

'AATCCTGG',

'GCTACTTG',

'ATGTTTTG']        
last_motifs = RandomizedMotifSearch_k(Dna_list=dna_strings, k=3, t=4)
i = 0
while i<1:
    best_motifs = RandomizedMotifSearch_k(Dna_list=dna_strings, k=3, t=4)
    
    if Score(best_motifs) < Score(last_motifs):
        last_motifs = best_motifs
        
    i +=1
print(last_motifs)