# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 13:06:40 2020

@author: femiogundare
"""


def StringComposition(Text, k):
    """
    Generate the k-mer composition of a string
    
    Input: A string Text and an integer k.
    Output: COMPOSITIONk(Text), where the k-mers are arranged in lexicographic order.
    
    Sample Input: 5; CAATCCAAC
    Sample Output: CAATC; AATCC; ATCCA; TCCAA; CCAAC
    """
    
    kmers = []
    
    for i in range(len(Text)-k+1):
        kmers.append(Text[i : i+k])
        
    string_composition = sorted(kmers)
    
    return string_composition


def PathToGenome(paths):
    """
    Output the Genome given its paths
    """
    first_path = paths[0]
    other_paths = paths[1:]
    
    last_nucleotide_of_other_paths = [path[-1] for path in other_paths]
    
    genome = first_path + ''.join(last_nucleotide_of_other_paths)
    
    return genome


def OverlapGraph(Patterns):
    """
    Construct the overlap graph of a collection of k-mers
    
    Input: A collection Patterns of k-mers.
    Output: The overlap graph Overlap(Patterns)
    
    Sample Input: 
        ATGCG
        GCATG
        CATGC
        AGGCA
        GGCAT
        GGCAC
    Sample Output:
        CATGC -> ATGCG
        GCATG -> CATGC
        GGCAT -> GCATG
        AGGCA -> GGCAC,GGCAT
    
    To solve this problem, form a node for each k-mer in Patterns and connect k-mers Pattern and Pattern' by a 
    directed edge if Suffix(Pattern) is equal to Prefix(Pattern'). Th resulting graph is the OVERLAP GRAPH
    """
    
    adjacent_list = {}
    
    for node in Patterns:
        adjacent_list[node] = []
        
    for node in adjacent_list.keys():
        for kmer in Patterns:
            if node[1:]==kmer[:-1]:
                adjacent_list[node].append(kmer)
                
    return adjacent_list


def DeBruijnGraph(Text, k):
    """
    Construct the De Bruijn Graph of a string Text
    
    Input: An integer k and a string Text
    Output: DeBruijnk(Text)
    
    Sample Input:
        4
        AAGATTCTCTAAGA
        
    Sample Output:
        AAG -> AGA,AGA
        AGA -> GAT
        ATT -> TTC
        CTA -> TAA
        CTC -> TCT
        GAT -> ATT
        TAA -> AAG
        TCT -> CTA,CTC
        TTC -> TCT
    """
    
    adjacent_list = {}
    edges = []
    
    for i in range(len(Text)-k+1):
        edge = Text[i : i+k]
        edges.append(edge)
    
    edges = sorted(edges)
    nodes = [pattern[:-1] for pattern in edges]
    
    for node in nodes:
        adjacent_list[node] = []
        
    for node in adjacent_list.keys():
        for edge in edges:
            if node==edge[:-1]:
                adjacent_list[node].append(edge[1:])
    
    return adjacent_list


def DeBruijnGraphFromKmers(kmers):
    """
    Construct the De Bruijn Graph of a set of k-mers
    
    Input: A collection of k-mers Patterns.
     Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).
    
    Sample Input:
        GAGG
        CAGG
        GGGG
        GGGA
        CAGG
        AGGG
        GGAG
        
    Sample Output:
        AGG -> GGG
        CAG -> AGG,AGG
        GAG -> AGG
        GGA -> GAG
        GGG -> GGA,GGG
    """
    
    adjacent_list = {}
    
    edges = sorted(kmers)
    nodes = [pattern[:-1] for pattern in edges]
    
    for node in nodes:
        adjacent_list[node] = []
        
    for node in adjacent_list.keys():
        for edge in edges:
            if node==edge[:-1]:
                adjacent_list[node].append(edge[1:])
    
    return adjacent_list