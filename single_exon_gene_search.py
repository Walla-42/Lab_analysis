from Bio.SeqUtils import seq3
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
import re

def find_orfs(dna_seq):
    """Function to find all open reading frames in a DNA sequence"""
    standard_table = CodonTable.unambiguous_dna_by_id[1]
    start_codons = ["ATG", "TTG", "GTG"]  # ATG is the main start codon but TTG, GTG are possible start codons in virus, parasites and bacteria
    stop_codons = standard_table.stop_codons
    orfs = []
    seq_len = len(dna_seq)
    
    for frame in range(3):
        i = frame
        while i < seq_len - 2:
            codon = dna_seq[i:i+3]
            if codon in start_codons:
                for j in range(i, seq_len - 2, 3):
                    stop_codon = dna_seq[j:j+3]
                    if stop_codon in stop_codons:
                        orf = dna_seq[i:j+3]
                        orfs.append(orf)
                        break
            i += 3
    return orfs


def translate_dna(dna_seq):
    """Function to translate DNA sequence to protein sequence. Input DNA sequence you would like to translate. Output is a string of amino acids. """
    standard_table = CodonTable.unambiguous_dna_by_id[1]
    return dna_seq.translate(table=standard_table, to_stop=True)



def find_all_positions(string, substring):
    """ Function to find all positions of a substring in a string"""
    start = 0
    while start < len(string):
        start = string.find(substring, start)
        if start == -1:
            break
        yield start
        start += 1


def search_gene(dna_seq, protein_seq):
    """Function for searching genomes for specific genes via their proteins. Only works on genes with no exons/introns or genes with only one exon. """
    protein_regex = re.compile(protein_seq.replace('N', '.'))
    orfs = find_orfs(dna_seq)
    match_found = False  

    for orf in orfs:
        translated_seq = translate_dna(Seq(orf))
        possible_protein_seqs = [str(translated_seq)]

        for protein_seq in possible_protein_seqs:
            if protein_regex.search(protein_seq):
                print(f"Found gene at position: {dna_seq.find(orf)}")
                match_found = True  

    if not match_found:
        print("No matching gene sequences found.")


def select_chromosome(fasta_file):    
    from Bio import SeqIO

    """User guided sequence that parses a .fna or .fasta file according to how the user would like it analyzed. This program was tested on the 
    genome of the plasmodium falciparum to select specific chromosomes for analysis and export them to a seperate .fna file. """

    while True:
        # fasta_file = input('Enter file path: ')
        
        if fasta_file[-4:] == '.fna' or fasta_file[-6:] == '.fasta':
            print('\033[1;33m' + 'Opening and parsing file...' + '\033[0m')
            break

        else:
            print('\033[1;31m' + "file type not supported. Please select a .fna or .fasta file to read." + '\033[0m')

    for record in SeqIO.parse(fasta_file, 'fasta'):
        print(f"Description: {record.description}")

    while True:
        select = input('Which chromosome would you like to analyze? ')
        file_name = 'selected_chromosome.fna'
        match_found = False

        for record in SeqIO.parse(fasta_file, 'fasta'):
            if select in record.description:
                SeqIO.write(record, file_name, 'fasta')

                match_found =True

                confirm = '\033[1;33m' + f"The sequence above has been selected and saved as {file_name}" + '\033[0m'
                print(f"ID: {record.id}")
                print(f"Name: {record.name}")
                print(f"Description: {record.description}")
                print(f"Length: {len(record.seq)}")
                break

        if match_found:
            break

        else:
            print('\033[1;31m' + "No matching records found. Please select a valid record." + '\033[0m')
    return confirm


print(select_chromosome('plasmodiumfalciparum_genome.fna'))

Target=  "MASAPGLAFANITLMLDLPQLPAIFFVNVRNNFKIFMNEIKQKTVEGEDIFYPHNRINLQNKQINKMGRTRKYSNNKEWIFGNPF"
# This target protein is a protein sequence in the above genome coding for an ATP synthase-associated protein. It has 1 exon and will work with this code, returning 'Found gene at position: 302825
# The above gene is located on chromosome 3. Just copy and paste the NC_000521.4 Plasmodium falciparum 3D7 genome assembly, chromosome: 3 when prompted. 

fasta_file = 'selected_chromosome.fna'


for record in SeqIO.parse(fasta_file, 'fasta'):
    print('\033[1;33m' + f"Record being analyzed: {record.description}" + '\033[0m')
    search_gene(record.seq, Target)