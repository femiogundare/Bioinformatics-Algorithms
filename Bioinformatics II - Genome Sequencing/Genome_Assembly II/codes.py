# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 23:11:45 2020

@author: femiogundare
"""

import random
import itertools
from random import shuffle

"""
https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
"""



def PathToGenome(paths):
    """
    Output the Genome given its paths
    """
    first_path = paths[0]
    other_paths = paths[1:]
    
    last_nucleotide_of_other_paths = [path[-1] for path in other_paths]
    
    genome = first_path + ''.join(last_nucleotide_of_other_paths)
    
    return genome


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





def EulerianCircuit(directed_graph):
    """
    Returns the Eulerian Circuit in a given directed graph.
    
    Input: The adjacency list of an Eulerian directed graph.
    Output: An Eulerian cycle in this graph.
    
    Sample Input:
        0 -> 3
        1 -> 0
        2 -> 1,6
        3 -> 2
        4 -> 2
        5 -> 4
        6 -> 5,8
        7 -> 9
        8 -> 7
        9 -> 6
    
    Sample Output:
        6->8->7->9->6->5->4->2->1->0->3->2->6
    """
    
    all_start_nodes = [k for k, v in directed_graph.items()]
    all_end_nodes = [c for k, v in directed_graph.items() for c in v]
    
    monitor = 0
    for node in list(set(all_start_nodes + all_end_nodes)):
        count_in = all_end_nodes.count(node);
        try:
            count_out = len(directed_graph[node]);
        except KeyError:
            count_out = 0
            directed_graph[node] = []
        if count_in == count_out:    #check if the in-degree of node is equal to the out-degree of node
            monitor +=1
        else:
            print('Graph lacks a Eulerian Circuit!')
            break
    
    #ensures that all nodes have same in-degree as out-degree
    assert monitor == len(list(set(all_start_nodes + all_end_nodes)))   
    
    
    stack = []   #to keep the nodes
    circuit = []
    
    nodes = list(directed_graph.keys())
    random_node = random.choice(nodes)   #we can start from any node, so I select a node randomly
    stack = [random_node]    #add the randomly selected node to the stack
    
    while stack:
        current_node = stack[-1]    #pick the last vertex and set it as the current one
        
        #if the current vertex has outgoing edges(i.e neighbors) remaining
        if directed_graph[current_node]:    
            #find and remove the next vertex that is adjacent to the current node
            next_node = directed_graph[current_node].pop()
            #add the new node to the stack
            stack.append(next_node)
        
        #if the current node has no outgoing edges(i.e neighbors) remaining
        else:
            #add it to the circuit
            circuit.append(stack.pop())   
    
    return circuit[::-1]





def EulerianPath(directed_graph):
    """
    Returns the Eulerian Path in a given directed graph.
    
    Input: The adjacency list of a directed graph that has an Eulerian path.
    Output: An Eulerian path in this graph.
    
    Sample Input:
        0 -> 2
        1 -> 3
        2 -> 1
        3 -> 0,4
        6 -> 3,7
        7 -> 8
        8 -> 9
        9 -> 6
    
    Sample Output:
        6->7->8->9->6->3->0->2->1->3->4
    """
    
    all_start_nodes = [k for k, v in directed_graph.items()]
    all_end_nodes = [c for k, v in directed_graph.items() for c in v]
    
    start_node, end_node = None, None
    
    for node in list(set(all_start_nodes + all_end_nodes)):
        count_in = all_end_nodes.count(node);
        try:
            count_out = len(directed_graph[node]);
        except KeyError:
            count_out = 0
            directed_graph[node] = []
        if count_in > count_out:
            end_node = node
            print('End node: ', end_node)
            #triggered = True
        elif count_in < count_out:
            start_node = node
            print('Start node: ', start_node)
            #trigger = True
    #if trigger == False:
        #print('This Graph is perfectly balanced!')
        
        
    stack = []   #to keep the nodes
    path = []
    
    nodes = list(directed_graph.keys())
    random_node = random.choice(nodes)   #we can start from any node, so I select a node randomly
    stack = [start_node]    #add the randomly selected node to the stack
    
    while stack:
        current_node = stack[-1]    #pick the last vertex and set it as the current one
        
        #if the current vertex has outgoing edges(i.e neighbors) remaining
        if directed_graph[current_node]:    
            #find and remove the next vertex that is adjacent to the current node
            next_node = directed_graph[current_node].pop()
            #add the new node to the stack
            stack.append(next_node)
        
        #if the current node has no outgoing edges(i.e neighbors) remaining
        else:
            #add it to the circuit
            path.append(stack.pop())   
    
    return path[::-1]





def StringReconstruction(patterns):
    """
    Solve the String Reconstruction Probelem
    
    Input: An integer k followed by a list of k-mers Patterns.
    Output: A string Text with k-mer composition equal to Patterns. (If multiple answers exist, you may return any one.)
    
    Sample Input:
        4
        CTTA
        ACCA
        TACC
        GGCT
        GCTT
        TTAC
        
    Sample Output:
        GGCTTACCA
    """
    
    deBruijn = DeBruijnGraphFromKmers(kmers=patterns)    #constructs the De Bruijn graph of the set of kmers
    path = EulerianPath(deBruijn)                        #outputs the Eulerian Path in the De Bruijn graph
    Text = PathToGenome(paths=path)                      #converts the Eulerian Path to the Genome string
    
    return Text

#t = ['AAAT', 'AATG', 'ACCC', 'ACGC', 'ATAC', 'ATCA', 'ATGC', 'CAAA', 'CACC', 'CATA', 'CATC', 'CCAG', 'CCCA',
# 'CGCT', 'CTCA', 'GCAT', 'GCTC', 'TACG', 'TCAC', 'TCAT', 'TGCA'
# ]
#print(StringReconstruction(patterns=t))


def k_UniversalCircularString(k):
    """
    Find a k-universal circular string.
    
    Input: An integer k.
    Output: A k-universal circular string.
    """
    
    data = itertools.product('01', repeat=k)
    l_data = [''.join(value) for value in list(data)]
    
    shuffle(l_data)
    
    deBruijn = DeBruijnGraphFromKmers(kmers=l_data)
    circuit = EulerianCircuit(deBruijn)[:-(k-1)]
    Text = PathToGenome(paths=circuit)
    
    return Text




def StringSpelledByGappedPatterns(k, d, gapped_patterns):
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
    
    read_1 = ''
    read_2 = ''
    
    #Add all first letters from the Left and Right sides
    #Do not include last line from both Left and Right sides
    for k_d_mer in gapped_patterns[:-1]:
        read_1 += k_d_mer[0][0]
        read_2 += k_d_mer[1][0]
    
    #add only last kmers (whole kmer) from both left and right sides
    read_1 += gapped_patterns[-1][0]
    read_2 += gapped_patterns[-1][1]
    
    #store the length of the string to be outputed (i.e string to be expected by combining the reads) ---- (k+d+k+n-1)
    str_length = k + d + k + len(gapped_patterns) -1
    
    #store the length of each read ===> Read_1 = Read_2
    read_length = len(read_1)
    
    #store the difference between the length of the read and the string to be expected by combining the reads
    match_length = str_length - read_length
    
    #search for overlap between the two strings
    if read_1[match_length:] == read_2[:-match_length]:
        return read_1 + read_2[-match_length:]
    
    return None