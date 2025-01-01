from Bio import SeqIO

"""User guided sequence that parses a .fna or .fasta file according to how the user would like it analyzed. This program was tested on the 
genome of the plasmodium falciparum to select specific chromosomes for analysis and export them to a seperate .fna file. This code is part 
of a much larger program to search for specific sequences in a genome which I am currently working on. """

while True:
    fasta_file = input('Enter file path: ')
   
    if fasta_file[-4:] == ('.fna' or '.fasta'):
        print('\033[1;33m' + 'Opening and parsing file...' + '\033[0m')
        break

    else:
        print('\033[1;31m' + "file type not supported. Please select a .fna or .fasta file to read." + '\033[0m')

for record in SeqIO.parse(fasta_file, 'fasta'):
    print(f"Description: {record.description}")

while True:
    select = input('Which chromosome would you like to analyze? ')
    file_name = input('What would you like to name your sequence? ') + '.fna'
    match_found = False

    for record in SeqIO.parse(fasta_file, 'fasta'):
        if select in record.description:
            SeqIO.write(record, file_name, 'fasta')

            match_found =True

            print('\033[1;33m' + f"The sequence below has been selected and saved as {file_name}" + '\033[0m')
            print(f"ID: {record.id}")
            print(f"Name: {record.name}")
            print(f"Description: {record.description}")
            print(f"Length: {len(record.seq)}")
            break

    if match_found:
        break

    else:
        print('\033[1;31m' + "No matching records found. Please select a valid record." + '\033[0m')
