from collections import Counter
import sys

# /home/walla42/python/Bioinformatics/DNA_pipeline/DNA-pipeline/selected_chromosome.fna
# Dictionary for functions
def frequency_map(text, k):
    """Defines frequency of repeating sequences of length 'k' in sequence 'text'.
    Can also be used on normal text, not just sequence data. """
    freq = {}
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i:i + k]
        freq[pattern] = 0
        for j in range(n - k + 1):
            if text[j:j + k] == pattern:
                freq[pattern] += 1
    return freq

def frequent_words(text, k):
    words = []
    freq = frequency_map(text, k)
    max_freq = max(freq.values())
    for key, value in freq.items():
        if value == max_freq:
            words.append(key)
    return words

def count_bases(seq):
    """Counts the number of each base in a sequence for RNA or DNA."""
    return Counter(seq)

def reverse_complement_dna(dna_sequence):
    """Gives the reverse complement sequence for DNA sequences."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna_sequence))

def transcription(dna_sequence):
    """Transcribes DNA to RNA."""
    return dna_sequence.replace('A', 'U').replace('T', 'A').replace('C', 'G').replace('G', 'C')

def translate(seq):
    """Translates RNA to Protein. Does not account for start and stop codons or different sequence frames."""
    table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
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
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table.get(codon, "")
    return protein

def open_and_parse_fasta(filepath):
    """Opens and parses a .fna or .fasta file, returning a dictionary with labels as keys and sequences as values."""
    fasta_dict = {}
    with open(filepath, "r") as fasta_file:
        label = ""
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                label = line[1:]
                fasta_dict[label] = ""
            else:
                fasta_dict[label] += line
    return fasta_dict

def process_file():
    """Handles user input for file directory and processes the file accordingly."""
    while True:
        directory = input('What is your file directory: ').strip()

        if directory.endswith((".fna", ".FASTA", ".seq")):
            return open_and_parse_fasta(directory)

        elif directory.endswith(".txt"):
            print('\033[1;31m' + "File type not supported..." + '\033[0m')

        else:
            print('\033[1;31m' + "File type not supported, please use another file.\n Supported file types include: .fna, .FASTA, .seq" + '\033[0m')

def main():
    """Main function to handle user interaction and sequence analysis."""
    sequence_data = process_file()

    print('\033[1;33m' + "Sequence Received..." + '\033[0m')

    for label in sequence_data:
        print("This is your reference sequence name: \n" + '\033[1;32m' + label + '\033[0m')
        header = label
        break

    sequence = sequence_data[header]

    while True:
        sequence_type = input('Is your sequence DNA or RNA? ').upper()
        if sequence_type in ["DNA", "RNA"]:
            break
        else:
            print('\033[1;31m' + "Please enter a valid response..." + '\033[0m')

    while True:
        if sequence_type == "DNA":
            work_type = input('Would you like a base count, kmer analysis, transcription, reverse complement strand, or would you like to try a new sequence? Else type END. \n').lower()

            if work_type == "base count":
                print("This is the count for each base in your sequence:\n " + '\033[1;32m' + str(count_bases(sequence)) + '\033[0m')

            elif work_type == "kmer analysis":
                k = int(input('How many bases should be in your pattern? '))
                print('\033[1;33m' + "One moment....processing request" + '\033[0m')
                print("These are the most frequent repeats: " + '\033[1;32m' + str(frequent_words(sequence, k)) + '\033[0m')

            elif work_type == "transcription":
                rna = transcription(sequence)
                print('\033[1;33m' + "Transcribed RNA sequence:\n" + '\033[0m' + rna)
                save_sequence(header, rna, "RNA Sequence")

            elif work_type == "reverse complement strand":
                rev_complement = reverse_complement_dna(sequence)
                print('\033[1;33m' + "Reverse complement sequence:\n" + '\033[0m' + rev_complement)
                save_sequence(header, rev_complement, "Reverse Complement Strand")

            elif work_type in ["end", "new sequence"]:
                break

            else:
                print('\033[1;31m' + "Please enter a valid response..." + '\033[0m')

        elif sequence_type == "RNA":
            work_type = input('Would you like a base count, translation, or the DNA Sequence? ').lower()

            if work_type == 'base count':
                print("This is the count for each base in your sequence:\n " + count_bases(sequence))

            elif work_type == 'translation':
                print("Translated RNA:\n " + translate(sequence))

            elif work_type == 'dna sequence':
                rev_dna = reverse_complement_dna(sequence)
                print('\033[1;33m' + "DNA sequence from RNA:\n" + '\033[0m' + rev_dna)
                save_sequence(header, rev_dna, "DNA Sequence")

            elif work_type in ['stop', 'end']:
                break

            else:
                print('\033[1;31m' + "Please enter a valid response..." + '\033[0m')

def save_sequence(header, sequence, description):
    """Prompts user to save the sequence to a file."""
    while True:
        answer = input("Would you like to save this sequence? ").lower()
        if answer == "yes":
            file_name = input("What should I name the file? ") + ".fna"
            file_path = input("Where would you like to save the file? ")
            with open(f"{file_path}/{file_name}", 'a') as file:
                file.write(f">{header} {description}\n{sequence}\n")
                print('\033[1;43m' + "Sequence has been saved in the specified directory." + '\033[0m')
            break
        elif answer == "no":
            print(sequence)
            print('\033[1;43m' + "Sequence not saved..." + '\033[0m')
            break
        else:
            print('\033[1;31m' + "Please enter a valid response..." + '\033[0m')

if __name__ == "__main__":
    main()
