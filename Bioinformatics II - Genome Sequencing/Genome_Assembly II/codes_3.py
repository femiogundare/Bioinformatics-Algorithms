# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:26:47 2020

@author: femiogundare
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 22:43:25 2020

@author: femiogundare
"""


from codes import PathToGenome, EulerianPath, StringSpelledByGappedPatterns



def DeBruijnGraphFromReadPairs(ReadPairs):
    """
    Construct the De Bruijn Graph of a set of read-pairs
    
    Input: A collection of (k, d)-mers
    Output: The adjacency list of the de Bruijn graph DeBruijn(k_d_mers).
    
    Sample Input:
        ACGT|CGTT
        TCGT|GTGA
        
    Sample Output:
        ACG|CGT -> CGT|GTT
        TCG|GTG -> CGT|TGA
    """
    
    prefix_of_read_1, prefix_of_read_2 = [], []
    suffix_of_read_1, suffix_of_read_2 = [], []
    
    for read_pair in ReadPairs:
        read_1, read_2 = read_pair[0], read_pair[1]
        
        prefix_of_read_1.append(read_1[:-1])
        prefix_of_read_2.append(read_2[:-1])
        
        suffix_of_read_1.append(read_1[1:])
        suffix_of_read_2.append(read_2[1:])
        
    prefix = [i + '|' + j for i, j in zip(prefix_of_read_1, prefix_of_read_2)]
    suffix = [i + '|' + j for i, j in zip(suffix_of_read_1, suffix_of_read_2)]
    
    adjacent_list = {x:[y] for x, y in zip(prefix, suffix)}
    
    return adjacent_list    

read_pairs = [['GAGA', 'TTGA'], ['TCGT', 'GATG'], ['CGTG', 'ATGT'], ['TGGT', 'TGAG'], ['GTGA', 'TGTT'], 
              ['GTGG', 'GTGA'], ['TGAG', 'GTTG'], ['GGTC', 'GAGA'], ['GTCG', 'AGAT']]





def StringSpelledByGappedPatterns(gapped_patterns, k, d):
    """
    Implement StringSpelledByGappedPatterns.
    
    Input: Integers k and d followed by a sequence of (k, d)-mers (a1|b1), … , (an|bn) such that Suffix(ai|bi) = Prefix(ai+1|bi+1) for 1 ≤ i ≤ n-1.
    Output: A string Text of length k + d + k + n - 1 such that the i-th (k, d)-mer in Text is equal to (ai|bi)  for 1 ≤ i ≤ n (if such a string
     exists).
    
    Sample Input:
        4 2
        GACC|GCGC
        ACCG|CGCC
        CCGA|GCCG
        CGAG|CCGG
        GAGC|CGGA
    
    Sample Output:
        GACCGAGCGCCGGA
    """
    
    first_patterns = []
    second_patterns = []
    
    for gapped_pattern in gapped_patterns:
        first_patterns.append(gapped_pattern[0])
        second_patterns.append(gapped_pattern[1])
        
    PrefixString = PathToGenome(paths=first_patterns)
    SuffixString = PathToGenome(paths=second_patterns)
    
    length_of_path = k+d+1
    
    for i in range(length_of_path, len(PrefixString)):
        if PrefixString[i] != SuffixString[i-k-d]:
            return 'there is no string spelled by the gapped patterns'
    
    return PrefixString + SuffixString[-(k+d):]




def StringReconstructionFromReadPairs(k, d, read_pairs):
    """
    Solve the String Reconstruction from Read-Pairs Problem.
    
    Input: Integers k and d followed by a collection of paired k-mers PairedReads.
    Output: A string Text with (k, d)-mer composition equal to PairedReads.
     
    Sample Input:
        4 2
        GAGA|TTGA
        TCGT|GATG
        CGTG|ATGT
        TGGT|TGAG
        GTGA|TGTT
        GTGG|GTGA
        TGAG|GTTG
        GGTC|GAGA
        GTCG|AGAT
         
    Sample Output:
        GTGGTCGTGAGATGTTGA
    """
    
    deBruijn = DeBruijnGraphFromReadPairs(ReadPairs=read_pairs)
    eulerian_path = EulerianPath(directed_graph=deBruijn)
    
    gapped_patterns = [gapped_pattern.split('|') for gapped_pattern in eulerian_path]
    
    string_spelled = StringSpelledByGappedPatterns(gapped_patterns=gapped_patterns, k=k, d=d)
    
    return string_spelled



t = ['ACC|ATA', 'ACT|ATT', 'ATA|TGA', 'ATT|TGA', 'CAC|GAT', 'CCG|TAC', 'CGA|ACT', 'CTG|AGC', 'CTG|TTC',
     'GAA|CTT', 'GAT|CTG', 'GAT|CTG', 'TAC|GAT', 'TCT|AAG', 'TGA|GCT', 'TGA|TCT', 'TTC|GAA'
     
     ]





print(StringReconstructionFromReadPairs(k=3, d=1, read_pairs=t))