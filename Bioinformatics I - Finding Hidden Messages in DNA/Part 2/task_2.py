# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 20:47:20 2020

@author: femiogundare
"""

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


def ReverseComplement(text):
    mapper = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ''.join(mapper.get(ch, ch) for ch in text)
    reverse_complement = complement[::-1]
    
    return reverse_complement

def PatternToNumber(Pattern):
    """
    Sample Input: AGT
    Sample Output: 11
    """
    if len(Pattern)==0:
        return 0
    
    symbol_to_number = {'A':0, 'C':1, 'G':2, 'T':3}
    last_symbol = Pattern[-1]
    prefix = Pattern[:-1]
    return 4*PatternToNumber(prefix) + symbol_to_number[last_symbol]


#*******************************************************************************************************

def MinimumSkew(Genome):
    """
    Finds a position in a genome where the skew diagram attains a minimum
    Input: A DNA string Genome
    Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|)
    
    Sample Input: TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT
    Sample Output: 11 24
    """
    position_of_minimum_skew = []
    skew_of_each_nucleotide = {}
    
    skew_of_each_nucleotide[0] = 0
    mapper = {'A':0, 'C':-1, 'G':1, 'T':0}
    
    for i in range(len(Genome)):
        if Genome[i]=='C':
            skew_of_each_nucleotide[i+1] = skew_of_each_nucleotide[i] + mapper['C']
        elif Genome[i]=='G':
            skew_of_each_nucleotide[i+1] = skew_of_each_nucleotide[i] + mapper['G']
        else:
            skew_of_each_nucleotide[i+1] = skew_of_each_nucleotide[i]
    
    skew_values = skew_of_each_nucleotide.values()
    minimum_skew = min(skew_values)
    position_of_minimum_skew = [str(p) for p, v in skew_of_each_nucleotide.items() if v==minimum_skew]
    
    return ' '.join(position_of_minimum_skew)


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


def ApproximatePatternMatching(Text, Pattern, d):
    """
    Find all approximate occurrences of a pattern in a string
    
    Sample Input:
        ATTCTGGA
        CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
        3
    Sample Output: 6 7 26 27
    """
    positions = []
    
    for i in range(len(Text)-len(Pattern)+1):
        pattern = Text[i : i+len(Pattern)]
        
        if HammingDistance(p=pattern, q=Pattern) <= d:
            positions.append(i)
            
    positions = ' '.join(str(x) for x in positions)
    
    return positions


"""
The difference between PatternMatch(text, pattern) and ApproximatePatternMatching(Text, Pattern, d) is 
quite simple.

* PatternMatch takes in a DNA string and a Pattern and looks for the positions where that Pattern can be
  found in the DNA string
* ApproximatePatternMatching takes in a DNA string, a Pattern and a mismatch value d. It also looks for
  the positions where the given Pattern can be found in the DNA string. However, it takes cognizance of
  MISMATCH...
  Let's say we our Pattern is ATGCCGTA and on searching the DNA string, we find a pattern ATGCCGAA, our
  ApproximatePatternMatching will include the position of pattern. It all depends on the value of d
  d=1 means 1 mismatch i.e we are to allow our algorithm to store positions where the input Pattern corres
      ponds with pattern with only 1 mismatch/error
  and so on
"""


def ApproximatePatternCount(Text, Pattern, d):
    """
    Computes the number of times a PARTICULAR k-mer appears in a DNA string, taking mismatch into consideration
    
    Sample Input: GAGG; TTTAGAGCCTTCAGAGG; 2
    Sample Output: 4
    """
    count = 0
    
    for i in range(len(Text)-len(Pattern)+1):
        pattern = Text[i : i+len(Pattern)]
        
        if HammingDistance(p=pattern, q=Pattern) <= d:
            count += 1
        
    return count


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



def ComputingFrequenciesWithMismatches(Text, k, d):
    """
    Computes the number of times every mismatch k-mer of each k-mer in a DNA string appears in that string
    """
    FrequencyArray = {}
    
    for i in range(4**k):   
        FrequencyArray[i] = 0
    for j in range(len(Text)-k+1):
        Pattern = Text[j : j+k]
        mismatch_patterns = GenerateMismatchedKmers(kmer=Pattern, d=d)
        
        for mismatch_pattern in mismatch_patterns:
            j = PatternToNumber(mismatch_pattern)
            FrequencyArray[j] = FrequencyArray[j] + 1
            
    FrequencyArray = [str(k) for k in FrequencyArray.values()]
    FrequencyArray = ' '.join(FrequencyArray)
    
    return FrequencyArray


"""
ComputingFrequencies(Text, k) is very simlar to ComputingFrequenciesWithMismatches(Text, k, d) except that:
    *the former computes the frequency of each k-mer in a DNA string;
        WHILE
    *the latter computes the frequency of the mismatch of each k-mer in a DNA string
"""



def FrequentWordsWithMismatches(Text, k, d):
    """
    Gets the mismatch k-mer(s) that occur the most  in a DNA string
    """
    FrequentPatterns = []
    #computes the number of times every k-mer in a DNA string, Text appears in Text
    FrequencyArray = ComputingFrequenciesWithMismatches(Text, k, d)  
    #puts the integers in the string into a list
    FrequencyArray = [int(x) for x in FrequencyArray.split()]   #puts the integers in the string into a list
    maxCount = max(FrequencyArray)
    
    for i in range(4**k):
        if FrequencyArray[i]==maxCount:   
            pattern = NumberToPattern(index=i, k=k)   
            FrequentPatterns.append(pattern)     
    FrequentPatterns = list(set(FrequentPatterns))  
    
    return ' '.join(FrequentPatterns)

"""
FasterFrequentWords(Text, k) is very simlar to FrequentWordsWithMismatches(Text, k, d) except that:
    *the former uses ComputingFrequencies(Text, k) while the latter uses 
     ComputingFrequenciesWithMismatches(Text, k, d)
"""



def ComputingFrequenciesWithMismatchesAndReverseComplement(Text, k, d):
    """
    Computes the number of times every mismatch k-mer and reverse complement of each k-mer 
    in a DNA string appears in that string
    """
    FrequencyArray = {}
    
    for i in range(4**k):   
        FrequencyArray[i] = 0
    for j in range(len(Text)-k+1):
        Pattern = Text[j : j+k]
        ReversedPattern = ReverseComplement(Pattern)
        mismatch_patterns = GenerateMismatchedKmers(kmer=Pattern, d=d)
        mismatch_reversed_patterns = GenerateMismatchedKmers(kmer=ReversedPattern, d=d)
        
        mismatch_patterns.extend(mismatch_reversed_patterns)
        
        for mismatch_pattern in mismatch_patterns:
            j = PatternToNumber(mismatch_pattern)
            FrequencyArray[j] = FrequencyArray[j] + 1
            
    FrequencyArray = [str(k) for k in FrequencyArray.values()]
    FrequencyArray = ' '.join(FrequencyArray)
    
    return FrequencyArray



def FrequentWordsWithMismatchesAndReverseComplement(Text, k, d):
   
    FrequentPatterns = []
    #computes the number of times every k-mer in a DNA string, Text appears in Text
    FrequencyArray = ComputingFrequenciesWithMismatchesAndReverseComplement(Text, k, d)  
    #puts the integers in the string into a list
    FrequencyArray = [int(x) for x in FrequencyArray.split()]   #puts the integers in the string into a list
    maxCount = max(FrequencyArray)
    
    for i in range(4**k):
        if FrequencyArray[i]==maxCount:   
            pattern = NumberToPattern(index=i, k=k)   
            FrequentPatterns.append(pattern)     
    FrequentPatterns = list(set(FrequentPatterns))  
    
    return ' '.join(FrequentPatterns)
