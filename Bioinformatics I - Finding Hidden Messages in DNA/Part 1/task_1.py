# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 23:28:56 2020

@author: femiogundare
"""


def PatternCount(text, pattern):
    """
    Computes the number of times a PARTICULAR k-mer appears in a DNA string
    Sample Input: GCGCG
    Sample Output: GCG  
    """
    count = 0
    
    for i in range(len(text)-len(pattern)+1):
        if (text[i : i+len(pattern)] == pattern):
            count +=1
        
    return count


def FrequentWords(text, k):
    """
    Gets the k-mer(s) that occur the most (i.e k-mers with the highest frequency) in a DNA string
    """
    frequent_patterns = []
    counts = []
    
    for i in range(len(text)-k+1):
        pattern = text[i : i+k]
        counts.append(PatternCount(text, pattern))  #computes the number of times the k-mer (pattern) appears in the string
    
    maxCount = max(counts)
    
    for i in range(len(text)-k+1):
        if counts[i] == maxCount:   #checks if count of the k-mer corresponds to the highest count in counts
            pattern = text[i : i+k]
            frequent_patterns.append(pattern)   #if yes, it appends that k-mer to the frequent_patterns list
    
    #selcts the k-mers of unqiue values cos repitition is bound to occur in the frequent_patterns list
    frequent_patterns = list(set(frequent_patterns))   
    
    return frequent_patterns


def ReverseComplement(text):
    mapper = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ''.join(mapper.get(ch, ch) for ch in text)
    reverse_complement = complement[::-1]
    
    return reverse_complement


def PatternMatch(text, pattern):
    """
    Find all occurrences of a pattern in a string
    """
    positions = []
    
    for i in range(len(text)-len(pattern)+1):
        if text[i : i+len(pattern)] == pattern:
            positions.append(i)  
    positions = ' '.join(str(x) for x in positions)
    
    return positions


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
        


def ComputingFrequencies(Text, k):
    """"
    Computes the number of times every k-mer in a DNA string appears in that string
    Sample Input: ACGCGGCTCTGAAA; 2
    Sample Output: 2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0
    """
    FrequencyArray = {}
    
    for i in range(4**k):   
        FrequencyArray[i] = 0
    for j in range(len(Text)-k+1):
        Pattern = Text[j : j+k]
        j = PatternToNumber(Pattern)
        FrequencyArray[j] = FrequencyArray[j] + 1
    FrequencyArray = [str(k) for k in FrequencyArray.values()]
    FrequencyArray = ' '.join(FrequencyArray)
    
    return FrequencyArray



def FasterFrequentWords(Text, k):
    """
    Gets the k-mer(s) that occur the most (i.e k-mers with the highest frequency) in a DNA string
    """
    FrequentPatterns = []
    #computes the number of times every k-mer in a DNA string, Text appears in Text
    FrequencyArray = ComputingFrequencies(Text, k)  
    #puts the integers in the string into a list
    FrequencyArray = [int(x) for x in FrequencyArray.split()]   #puts the integers in the string into a list
    maxCount = max(FrequencyArray)
    
    for i in range(4**k):
        if FrequencyArray[i]==maxCount:   #checks if count of the k-mer at position i corresponds to the highest count in FrequencyArray
            pattern = NumberToPattern(index=i, k=k)   #get the k-mer that corresponds to that particular index
            FrequentPatterns.append(pattern)     #adds the k-mer to FrequentPatterns
    FrequentPatterns = list(set(FrequentPatterns))  #selects only the unique values in FrequentPatterns
    
    return FrequentPatterns


"""
Differences between FrequentWords and FasterFrequentWords are:
    1. FasterFrequentWords utilizes the ComputingFrequencies(Text, k) fucntion which calculates the number of times
       every k-mer in a DNA string appears in that string, and then gets the maximum frequency outputed by
       ComputingFrequencies. It then checks which k-mer(s) corresponds to that maximum frequency.
                                           WHILE
        FrequentWords calls the function PatternCount(Text, k) |Text|-k+1 times, and has computational power of:
            |Text|-k+1 * (|Text|-k+1)*k
    2. FrequentWords is way more computationally expensive than FasterFrequentWords
"""


def ClumpFinding(Genome, k, t, L):
    """
    Gets the k-mers that form a (L, t)-clumps in the Genome i.e
    it gets the k-mers that appear at least t times in any window L of a long Genome string
    since these k-mers are appear in a window L of a Genome string at least t times, we say they form a CLUMP
    """
    genome_length = len(Genome)
    FrequentPatterns = []
    CLUMP = {}
    
    for i in range(4**k):    #taking into account that there are 4^k possible k-mers
        CLUMP[i] = 0         #setting the number of clumps formed by each of the 4^k possible k-mers to 0
        
    for i in range(genome_length-L+1):
        Text = Genome[i : i+L]    #obtain a string of window L of the total DNA string
        #computes the number of times every k-mer in Text appears in Text
        frequency_of_every_kmer = ComputingFrequencies(Text=Text, k=k)  
        frequency_of_every_kmer = [int(x) for x in frequency_of_every_kmer.split()]
        
        for index in range(4**k):
            if frequency_of_every_kmer[index]>=t:   #for every k-mer whose occurrence in any window L of the total DNA string is greater than t,
                #set the value of clump to 1 i.e the k-mer at theparticular index forms a clump in the window L of the total DNA string
                #k-mers that do not form clumps will still have their clumps tied to 0
                #so it's like if clump is formed, set k-mer position in clump to 1, else 0
                CLUMP[index] = 1                    
                
    for i in range(4**k):   #looping through the range of the possible k-mers (i.e 4^k)
        if CLUMP[i]==1:     #if the k-mer actually formed a clump in any window L of the DNA string,
            pattern = NumberToPattern(i, k)    #get the pattern/k-mer itself by converting its index
            #add the k-mer to the list of frequent k-mers since it satisfies the condition that 
            #its frequency in each of the window L must be greater than t
            FrequentPatterns.append(pattern)   
    FrequentPatterns = list(set(FrequentPatterns))  
    
    return ' '.join(FrequentPatterns)


def BetterClumpFinding(Genome, k, t, L):
    genome_length = len(Genome)
    FrequentPatterns = []
    CLUMP = {}
    
    for i in range(4**k):    #taking into account that there are 4^k possible k-mers
        CLUMP[i] = 0  
    
    Text = Genome[:L]
    frequency_of_every_kmer = ComputingFrequencies(Text=Text, k=k)
    frequency_of_every_kmer = [int(x) for x in frequency_of_every_kmer.split()]
    
    for i in range(4**k):
        if frequency_of_every_kmer[i]>=t:
            CLUMP[i] = 1
            
    for i in range(1, genome_length-L+1):
        FirstPattern = Genome[i-1 : i-1+k]
        index = PatternToNumber(FirstPattern)
        frequency_of_every_kmer[index] = frequency_of_every_kmer[index] - 1
        LastPattern = Genome[i+L-k : i+L-k+k]
        index = PatternToNumber(LastPattern)
        frequency_of_every_kmer[index] = frequency_of_every_kmer[index] + 1
        
        if frequency_of_every_kmer[index]>=t:
                CLUMP[index] = 1
                
    for i in range(4**k):
        if CLUMP[i]==1:
            pattern = NumberToPattern(i, k)
            FrequentPatterns.append(pattern)
    FrequentPatterns = list(set(FrequentPatterns))  
    
    return len(FrequentPatterns)
        