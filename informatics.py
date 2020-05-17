import argparse
import re
import math
from collections import Counter
from itertools import zip_longest, permutations

mass_table = {
    'A':71.03711,
    'C':103.00919,
    'D':115.02694,
    'E':129.04259,
    'F':147.06841,
    'G':57.02146,
    'H':137.05891,
    'I':113.08406,
    'K':128.09496,
    'L':113.08406,
    'M':131.04049,
    'N':114.04293,
    'P':97.05276,
    'Q':128.05858,
    'R':156.10111,
    'S':87.03203,
    'T':101.04768,
    'V':99.06841,
    'W':186.07931,
    'Y':163.06333,
}

def count_nucleotides(dna_string):
    c = Counter(dna_string)
    return c.get('A'), c.get('C'), c.get('G'), c.get('T')

def transcribe_dna(dna_string):
    return dataset.replace('T', 'U')

def reverse_complement(dna_string):
    result = ""
    complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    dataset = list(dataset)
    dataset.reverse()
    for char in dataset:
        result += complements[char]
    return result

def rabbit_pairs(time, litter):
    if time == 0:
        return 0
    elif time == 1:
        return 1
    else:
        return rabbit_pairs(time-1, litter) + litter*rabbit_pairs(time-2, litter)

def gc_content(dna_string):
    length = len(dna_string)
    c = Counter(dna_string)
    return ((c.get('G') + c.get('C')) / length) * 100

def FASTA_parser(dataset):
    dict_ = {}

    with open(dataset) as file_:
        data = file_.read().split()
    for line in data:
        if line[0] == ">":
            id = line[1:]
            dict_[line[1:]] = ""
        else:
            dict_[id] += line
    return dict_

def hamming_distance(dna_string1, dna_string2):
    hamming_distance = 0
    for n1, n2 in zip(dna_string1, dna_string2):
        if n1 != n2:
            hamming_distance += 1
    return hamming_distance

def chunks(list_, chunk_size):
    for i in range(0, len(list_), chunk_size):
        yield list_[i:i+chunk_size]

def rna_translation(rna_string):
    translated_rna = ""
    codons = [n1 + n2 + n3 for n1 in "UCAG" for n2 in "UCAG" for n3 in "UCAG"]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))

    for amino_acid in chunks(rna_string, 3):
        if codon_table[amino_acid] == "*":
            break
        translated_rna += codon_table[amino_acid]
    return translated_rna

def find_motif(dna_string, sub_string):
    lookahead = '(?='+ sub_string + ')'
    motifs = [str(m.start() + 1) for m in re.finditer(lookahead, dna_string)]
    return " ".join(motifs)

def protein_mass(protein_string):
    protein_mass = 0
    for amino_acid in protein_string:
        protein_mass += mass_table[amino_acid]
    return protein_mass

def consensus_string(dataset):
    data_ = FASTA_parser(dataset)
    matrix = [array for array in data_.values()]
    length = len(matrix[0])
    A = [0 for _ in range(length)]
    C = [0 for _ in range(length)]
    G = [0 for _ in range(length)]
    T = [0 for _ in range(length)]
    for i in range(length):
        for array in matrix:
            if array[i] == 'A':
                A[i] += 1
            if array[i] == 'C':
                C[i] += 1
            if array[i] == 'G':
                G[i] += 1
            if array[i] == 'T':
                T[i] += 1

    nac = ['A', 'C', 'G', 'T']
    values = []
    consensus_string = ""
    for i in range(length):
        for array in [A, C, G, T]:
            values.append(array[i])
        dict_ = dict(zip(nac, values))
        consensus_string += max(dict_, key=dict_.get)
        values = []
    return consensus_string

def gene_order(n):
    print(int(math.factorial(len(range(n))) / math.factorial(len(range(n)) - n)))
    for i in permutations(range(1, n+1), n):
        print(" ".join(map(str, i)))

