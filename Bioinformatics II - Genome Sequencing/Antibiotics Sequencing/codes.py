# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 13:10:47 2020

@author: femiogundare
"""

import copy
from collections import Counter


GeneticCodeArray = {
                        "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
                        "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
                        "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
                        "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
                        "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
                        "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                        "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
                        "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                        "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
                        "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                        "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                        "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                        "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
                        "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                        "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                        "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"
                        }


def ReverseComplement(text):
    mapper = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ''.join(mapper.get(ch, ch) for ch in text)
    reverse_complement = complement[::-1]
    
    return reverse_complement



def ProteinTranslationProblem(Rna_string):
    """
    Input: An RNA string Pattern and the array GeneticCode.
    Output: The translation of Pattern into an amino acid string Peptide.
    
    Sample Input:
        AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
    Sample Output:
        MAMAPRTEINSTRING
    """
    
    codon_length = 3
    amino_acid_of_each_codon = []
    
    for i in range(0, len(Rna_string), codon_length):
        codon = Rna_string[i : i+codon_length]
        
        if GeneticCodeArray[codon] != 'STOP':
            amino_acid_of_each_codon.append(GeneticCodeArray[codon])
        else:
            continue
        
    return ''.join(amino_acid_of_each_codon)



def PeptideEncodingProblem(Text, Peptide):
    """
    Find substrings of a genome encoding a given amino acid sequence.
    
    Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
    Output: All substrings of Text encoding Peptide (if any such substrings exist).
    
    Sample Input:
        ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
        MA
    Sample Output:
        ATGGCC
        GGCCAT
        ATGGCC
    """
    
    substring_length = 3*len(Peptide)   
    all_possible_substrings = []
    
    substrings_encoding_the_given_peptide_Peptide = []
    
    for i in range(len(Text)-substring_length+1):
        substring = Text[i : i+substring_length]
        all_possible_substrings.append(substring)
        
    for substring in all_possible_substrings:
        complement = substring
        reverse_complement = ReverseComplement(text=substring)
        
        rna_of_complement = complement.replace('T', 'U')
        rna_of_reverse_complement = reverse_complement.replace('T', 'U')
        
        peptide_of_rna_of_complement = ProteinTranslationProblem(Rna_string=rna_of_complement)
        peptide_of_rna_of_reverse_complement = ProteinTranslationProblem(Rna_string=rna_of_reverse_complement)
        
        if (peptide_of_rna_of_complement==Peptide) or (peptide_of_rna_of_reverse_complement==Peptide):
            substrings_encoding_the_given_peptide_Peptide.append(substring)
        else:
            continue
        
    return substrings_encoding_the_given_peptide_Peptide


#print(PeptideEncodingProblem(Text='ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', Peptide='MA'))



Alphabet = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
#Alphabet = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

AminoAcidMass = {
    'G': 57,
    'A': 71,
    'S': 87,
    'P': 97,
    'V': 99,
    'T': 101,
    'C': 103,
    'I': 113,
    'L': 113,
    'N': 114,
    'D': 115,
    'K': 128,
    'Q': 128,
    'E': 129,
    'M': 131,
    'H': 137,
    'F': 147,
    'R': 156,
    'Y': 163,
    'W': 186
}





def Subpeptides(peptide):
    l = len(peptide)
    ls = []
    looped = peptide + peptide
    for start in range(0, l):
        for length in range(1, l):
            ls.append((looped[start:start + length]))
    ls.append(peptide)
    return ls



def TheoreticalSpectrumProblem(Peptide):
    """
    Generate the theoretical spectrum of a cyclic peptide.
    
    Input: An amino acid string Peptide.
    Output: Cyclospectrum(Peptide).
    
    Sample Input:
        LEQN
    Sample Output:
        0 113 114 128 129 227 242 242 257 355 356 370 371 484
    """
    
    theoretical_spectrum = [0]
    
    subpeptides = Subpeptides(Peptide)
    
    for subpeptide in subpeptides:
        subpeptide_theoretical_spectrum = 0
        subpeptide = list(subpeptide)
        
        for s in subpeptide:
            subpeptide_theoretical_spectrum += AminoAcidMass[s]
        
        theoretical_spectrum.append(subpeptide_theoretical_spectrum)
        
    return sorted(theoretical_spectrum)



def LinearSpectrum(Peptide):
    """
    Implement Linear Spectrum.
    
    Input: An amino acid string Peptide.
    Output: The linear spectrum of Peptide.
    """
    
    PrefixMass = {}
    PrefixMass[0] = 0
    
    for i in range(len(Peptide)):
        for alphabet in Alphabet:
            if alphabet==Peptide[i]:
                PrefixMass[i+1] = PrefixMass[i] + AminoAcidMass[alphabet]
                
    linearSpectrum = [0]
                
    for i in range(len(PrefixMass)):
        for j in range(i+1, len(PrefixMass)):
            linearSpectrum.append(PrefixMass[j]-PrefixMass[i])
            
    return sorted(linearSpectrum)



def CyclicSpectrum(Peptide):
    """
    Implement Cyclic Spectrum.
    
    Input: An amino acid string Peptide.
    Output: The cyclic spectrum of Peptide.
    """
    
    PrefixMass = {}
    PrefixMass[0] = 0
    
    for i in range(len(Peptide)):
        for alphabet in Alphabet:
            if alphabet==Peptide[i]:
                PrefixMass[i+1] = PrefixMass[i] + AminoAcidMass[alphabet]
                
    peptideMass = PrefixMass[len(Peptide)]
    cyclicSpectrum = [0]
    
    for i in range(len(PrefixMass)):
        for j in range(i+1, len(PrefixMass)):
            cyclicSpectrum.append(PrefixMass[j]-PrefixMass[i])
            
            if i>0 and j<len(Peptide):
                cyclicSpectrum.append(peptideMass-(PrefixMass[j]-PrefixMass[i]))
                
    return sorted(cyclicSpectrum)



def expand(a, b):
    output = []
    for i in a:
        for j in b:
            output.append(i+j)
    return output
            
def counterSubset(list1, list2):
    c1, c2 = Counter(list1), Counter(list2)
    for k, n in c1.items():
        if n > c2[k]:
            return False
    return True
    
    
    
    
alphabets_to_remove = ['I', 'K']
Alphabet_modified = [alphabet for alphabet in Alphabet if alphabet not in alphabets_to_remove]

AminoAcidMassCopy = AminoAcidMass.copy()
del AminoAcidMassCopy['I'], AminoAcidMassCopy['K']


    
def CyclopeptideSequencing(Spectrum):
    """
    Implement Cyclopeptide Sequencing.
    Note: Multiple solutions may exist. You may return any one.
    
    Sample Input:
        0 113 128 186 241 299 314 427
    Sample Output:
        186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186
    """

    CandidatePeptides = [k for k,v in AminoAcidMassCopy.items() if v in Spectrum]
    print(len(CandidatePeptides))
    #print(CandidatePeptides)
    FinalPeptides = []
    
    while CandidatePeptides:
        CandidatePeptides = expand(CandidatePeptides, Alphabet_modified)
        #print(CandidatePeptides)
        CandidatePeptides_copy = copy.deepcopy(CandidatePeptides)
        #print(CandidatePeptides_copy)
        
        for Peptide in CandidatePeptides_copy:
            individual_amino_acids = list(Peptide)
            mass_of_peptide = 0
            for amino_acid in individual_amino_acids:
                mass_of_peptide += AminoAcidMassCopy[amino_acid]
            
            if mass_of_peptide == Spectrum[-1]:
                if (counterSubset(CyclicSpectrum(Peptide), Spectrum)==True) and Peptide not in FinalPeptides:
                    FinalPeptides.append(Peptide)
                CandidatePeptides.remove(Peptide)
                
            elif (counterSubset(LinearSpectrum(Peptide), Spectrum)==False):
                CandidatePeptides.remove(Peptide)
    
    output = []
    for peptide in FinalPeptides:
        individual_amino_acids = list(peptide)
        mass_of_individual_amino_acids = [AminoAcidMassCopy[amino_acid] for amino_acid in individual_amino_acids]
        joined_masses = '-'.join([str(m) for m in mass_of_individual_amino_acids])
        output.append(joined_masses)
        
                
    return output


#print(CyclopeptideSequencing([0, 113, 128, 186, 241, 299, 314, 427]))



def CyclopeptideScoring(Peptide, Spectrum):
    """
    Solve the Cyclopeptide Scoring Problem.
    
    Input: An amino acid string Peptide and a collection of integers Spectrum.
    Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
    
    Sample Input:
        NQEL
        0 99 113 114 128 227 257 299 355 356 370 371 484
    Sample Output:
        11
    """
    
    cyclicSpectrum = CyclicSpectrum(Peptide)
    score = 0
    
    for mass in Spectrum:
        if mass in cyclicSpectrum:
            score +=1
            cyclicSpectrum.remove(mass)
        
    return score



def LinearPeptideScoring(Peptide, Spectrum):
    """
    Solve the Linear Peptide Scoring Problem.
    
    Input: An amino acid string Peptide and a collection of integers Spectrum.
    Output: The score of Linear Peptide against Spectrum, Score(Peptide, Spectrum).
    
    Sample Input:
        NQEL
        0 99 113 114 128 227 257 299 355 356 370 371 484
    Sample Output:
        8
    """
    
    linearSpectrum = LinearSpectrum(Peptide)
    score = 0
    
    for mass in Spectrum:
        if mass in linearSpectrum:
            score +=1
            linearSpectrum.remove(mass)
        
    return score



def Trim(Leaderboard, Spectrum, N):
    """
    Implement Trim.
    
    Input: A collection of peptides Leaderboard, a collection of integers Spectrum, and an integer N.
    Output: The N highest-scoring linear peptides on Leaderboard with respect to Spectrum.
    
    Sample Input:
        LAST ALST TLLT TQAS
        0 71 87 101 113 158 184 188 259 271 372
        2
    Sample Output:
        LAST ALST
    """
    
    LinearScores = []
    for j in range(len(Leaderboard)):
        Peptide = Leaderboard[j]
        LinearScores.append(LinearPeptideScoring(Peptide=Peptide, Spectrum=Spectrum))
    
    Leaderboard = [x for _,x in sorted(zip(LinearScores, Leaderboard), reverse=True)]
    LinearScores = sorted(LinearScores, reverse=True)
    #print(Leaderboard)
    #print(LinearScores)
    
    
    for j in range(N, len(Leaderboard)):
        if LinearScores[j] < LinearScores[N-1]:
            del Leaderboard[j:]
            return Leaderboard
    
    return Leaderboard




def LeaderboardCyclopeptideSequencing(Spectrum, N):
    """
    Implement LeaderboardCyclopeptideSequencing.
    
    Input: An integer N and a collection of integers Spectrum.
    Output: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).
    
    Sample Input:
        10
        0 71 113 129 147 200 218 260 313 331 347 389 460
        
    Sample Output:
        113-147-71-129
    """
    
    Leaderboard = [k for k,v in AminoAcidMassCopy.items() if v in Spectrum]
    LeaderPeptide = ''
    
    while Leaderboard:
        Leaderboard = expand(Leaderboard, Alphabet_modified)
        Leaderboard_copy = copy.deepcopy(Leaderboard)
        
        for Peptide in Leaderboard_copy:
            individual_amino_acids = list(Peptide)
            mass_of_peptide = 0
            for amino_acid in individual_amino_acids:
                mass_of_peptide += AminoAcidMassCopy[amino_acid]
            
            if mass_of_peptide == Spectrum[-1]:
                if CyclopeptideScoring(Peptide, Spectrum) > CyclopeptideScoring(LeaderPeptide, Spectrum):
                    LeaderPeptide = Peptide
                
            elif mass_of_peptide > Spectrum[-1]:
                Leaderboard.remove(Peptide)
                
        Leaderboard = Trim(Leaderboard, Spectrum, N)
        
    #LeaderPeptide for sample dataset is LFAE
    LeaderPeptide = [AminoAcidMassCopy[amino_acid] for amino_acid in LeaderPeptide]
    LeaderPeptide = '-'.join([str(x) for x in LeaderPeptide])
    
    return LeaderPeptide


#print(LeaderboardCyclopeptideSequencing(
#    Spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460], 
#    N=10))


#xxx = '0 97 99 114 128 147 147 163 186 227 241 242 244 260 261 262 283 291 333 340 357 385 389 390 390 405 430 430 447 485 487 503 504 518 543 544 552 575 577 584 632 650 651 671 672 690 691 738 745 747 770 778 779 804 818 819 820 835 837 875 892 917 932 932 933 934 965 982 989 1030 1039 1060 1061 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1225 1322'
#Spectrum = [int(spectrum) for spectrum in xxx.split(' ')]
#print(LeaderboardCyclopeptideSequencing(Spectrum, 1000))

#print(CyclopeptideScoring(Peptide='VYQNFWPFLK', Spectrum=Spectrum))



def SpectralConvolution(Spectrum):
    """
    Compute the convolution of a spectrum.
    
    Input: A collection of integers Spectrum in increasing order..
    Output: The list of elements in the convolution of Spectrum. If an element has multiplicity k, it should appear exactly k times; you may return the elements in any order.
    
    Sample Input:
        0 137 186 323
    Sample Output:
        137 137 186 186 323 49
    """
    
    Spectrum = sorted(Spectrum)
    elements = []
    
    for i in range(len(Spectrum)):
        for j in range(len(Spectrum)):
            elements.append(Spectrum[i]-Spectrum[j])
            
    elements = [x for x in elements if x>0]
    
    return elements


#print(SpectralConvolution(Spectrum=[0, 137, 186, 323]))
    


def ConvolutionalCyclopeptideSequencing(M, N, Spectrum):
    """
    Implement ConvolutionCyclopeptideSequencing.
    
    Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
    Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements (and ties) of the convolution of Spectrum that fall between 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).

    Sample Input:
        20
        60
        57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493 
    Sample Output:
        99-71-137-57-72-57
    """
    
    spectral_convolution = SpectralConvolution(Spectrum)
    spectral_convolution = [k for k in spectral_convolution if k>=57 and k<200]
    
    List = Counter(spectral_convolution).most_common(M)
    frequency = List[-1][1]

    List1 = dict(Counter(spectral_convolution))
    ExtendedAlphabetMass = [k for k, v in List1.items() if v >= frequency]
    
    ExtendedAlphabet = []
    for i in range(0,len(ExtendedAlphabetMass)):
        ExtendedAlphabet.append(chr(ExtendedAlphabetMass[i]))
        
    global AminoAcid
    global AminoAcidMass
    Alphabet_modified = ExtendedAlphabet
    AminoAcidMassCopy = ExtendedAlphabetMass
    
    leader_peptide = LeaderboardCyclopeptideSequencing(Spectrum, N)
    
    return leader_peptide


#print(ConvolutionalCyclopeptideSequencing(
#    M=20, N=60, 
#    Spectrum=[57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]
#    ))


#print(Counter(SpectralConvolution(Spectrum=[0, 57, 118, 179, 236, 240, 301])))
print(CyclopeptideScoring(Peptide=' MAMA', Spectrum=[0, 71, 178, 202, 202, 202, 333, 333, 333, 404, 507, 507]))
#print(LinearPeptideScoring(Peptide='PEEP', Spectrum=[0, 97, 97, 129, 194, 196, 226, 226, 244, 258, 323, 323, 452]))