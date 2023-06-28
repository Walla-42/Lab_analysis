# This program is to showcase the programs that I have developed over the past week. This is a complex program that
# demonstrates strings, functions, and loops


# imports for program
from collections import Counter
import sys


# dictionary for functions
def count_bases_RNA():
    seq = RNA
    return Counter(seq)


def count_bases_DNA():
    nuc = DNA
    return Counter(nuc)


def reverse_sequencing(rna_sequence):
    reverse_sequence = rna_sequence.replace('A', 't').replace('U', 'a').replace(
        'G', 'c').replace('C', 'g').upper()
    return reverse_sequence


def reverse_complement_DNA(dna_sequence):
    reverse_complement = dna_sequence.replace('A', 't').replace(
        'T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
    return reverse_complement


def Transcription(dna_sequence):
    mRNA = dna_sequence.replace('A', 'u').replace('T', 'a').replace('C', 'g').replace(
        'G', 'c').upper()
    return mRNA


def translate(seq):
    table = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
             "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
             "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
             "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
             "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
             "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
             "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
             "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
             "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
             "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
             "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
             "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
             "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
             "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
             "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
             "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


# start of Program
while True:
    Sequence_type = input('Is your sequence DNA, or RNA? ')
    if Sequence_type == "DNA":
        break
    elif Sequence_type == "RNA":
        break
    else:
        print("Please enter valid sequence type.")

while True:
    if Sequence_type == "DNA":
        DNA = input("What is your DNA Sequence? ")
        Work_type = input('Would you like a base count, transcription, or the reverse compliment strand? ')
        if (Work_type == "base count" or Work_type == "Base Count"
                or Work_type == "base Count" or Work_type == "Base count"):
            print("This is the count for each base in your sequence:\n " + str(Counter(DNA)))
            break
        elif Work_type == "transcription" or Work_type == "Transcription":
            print("RNA Sequence:\n" + Transcription(DNA))
            break
        elif (Work_type == "Reverse Compliment Strand" or Work_type == "reverse compliment strand"
              or Work_type == "Reverse compliment Strand" or Work_type == "Reverse Compliment strand"
              or Work_type == "reverse Compliment Strand"):
            print("Reverse compliment DNA: \n" + reverse_complement_DNA(DNA))
            break
        else:
            print('Please enter a valid response.')
    else:
        break
while True:
    if Sequence_type == "RNA":
        RNA = input('What is your RNA Sequence? ')
        Work_type = input('Would you like a base count, translation, or the DNA Sequence? ')
        if (Work_type == 'Base Count' or Work_type == 'Base count' or Work_type == 'base count'
                or Work_type == 'base Count'):
            print("This is the count for each base in your sequence:\n " + str(count_bases_RNA()))
            break
        elif Work_type == 'Translation' or Work_type == 'translation':
            print("Translated RNA:\n " + translate(RNA))
            break
        elif Work_type == 'DNA Sequence' or Work_type == 'dna Sequence' or Work_type == 'DNA sequence':
            print("DNA Sequence:\n " + reverse_sequencing(RNA))
            break
        else:
            print("Please enter a valid response.")
    else:
        sys.exit()
