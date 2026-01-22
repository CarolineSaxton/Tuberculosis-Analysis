from Bio import Entrez, SeqIO
Entrez.email = "caroline.saxtonrowe@gmail.com" 

# Fetch the DNA sequence from the NCBI database
accession_number = "NC_000962.3" 
with Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text") as handle:
    record = SeqIO.read(handle, "fasta")
    dna_sequence = record.seq

print()
print(f"Gene Name: Mycobacterium tuberculosis H37Rv strain")
print(f"RefSeq Accession Number: {record.id}")
print(f"Total Nucleotides Length: {len(record.seq)}")
print(f"First 50 Nucleotides: {dna_sequence[:50]}")
print()
# Import a Counter to find and graph total number of adenine, thymine, guanine and cytosine nucleotides
from collections import Counter
nucleotide_counts = Counter(dna_sequence)
print(f"Total Count of Nucleotides: {nucleotide_counts}")

import matplotlib.pyplot as plt  # Import matplotlib for plotting
total_nucleotides = sum(nucleotide_counts.values())
nucleotide_percentages = {key: (value / total_nucleotides) * 100 for key, value in nucleotide_counts.items()}

plt.bar(
    nucleotide_percentages.keys(),  # x-axis: Nucleotide types (A, T, C, G)
    nucleotide_percentages.values(),  # y-axis: Percentage of each nucleotide
    color=['yellow', 'red', 'blue', 'green']  # Colors for each nucleotide bar
)
plt.title("Nucleotide Frequency")  # Provides context for the visualization
plt.xlabel("Nucleotide")  # Indicates the nucleotides being counted (A, T, C, G)
plt.ylabel("Frequency")  # Shows the frequency of each nucleotide
plt.show()  # Renders the bar chart

# Calculate GC Content
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

gc_content = calculate_gc_content(dna_sequence)

print(f"GC Content: {gc_content:.2f}%")  # Output example: "GC Content: 46.67%"
print()
# Codon to Amino Acid Conversions
amino_acid_codes = {
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",  # Alanine
    "TGT": "Cys", "TGC": "Cys",                          # Cysteine
    "GAT": "Asp", "GAC": "Asp",                          # Aspartic acid
    "GAA": "Gln", "GAG": "Gln",                          # Glutamic acid
    "TTT": "Phe", "TTC": "Phe",                          # Phenylalanine
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",  # Glycine
    "CAT": "His", "CAC": "His",                          # Histidine
    "ATA": "Ile", "ATT": "Ile", "ATC": "Ile",              # Isoleucine
    "AAA": "Lys", "AAG": "Lys",                          # Lysine
    "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",  # Leucine
    "ATG": "Met",                                      # Methionine (START codon)
    "AAT": "Asn", "AAC": "Asn",                          # Asparagine
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",  # Proline
    "CAA": "Glu", "CAG": "Glu",                          # Glutamine
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",  # Arginine
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",  # Serine
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",  # Threonine
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",  # Valine
    "TGG": "Trp",                                      # Tryptophan
    "TAT": "Tyr", "TAC": "Tyr",                          # Tyrosine
    "TAA": "STOP", "TAG": "STOP", "TGA": "STOP"               # STOP codons
}
# Find nucleotide and its coordinate, as well as  the cooresponding codon and amino acid
complete_genome_analysis = input("Would you like perform complete genome analysis? (yes/no): ")
if complete_genome_analysis == "yes":
    choice1 = input("Would you like to use a nucleotide coordinate (input:'n') or codon number (input:'c')?  ")
    if choice1 == "n":
        snp_position = int(input("Nucleotide Coordinate: ")) - 1
        print(f"Nucleotide: {dna_sequence[snp_position]}")
        codon_size = 3
        codon_sequence = [dna_sequence[i:i + codon_size] for i in range(0, len(dna_sequence), codon_size)]
        codon_number = (snp_position // 3) + 1
        print(f"Codon Number: {codon_number}")
        codon = codon_sequence[codon_number - 1]
        print(f"Codon: {codon}")
    if choice1 == "c":
        codon_number = int(input("Codon Number: "))
        codon_size = 3
        codon_sequence = [dna_sequence[i:i + codon_size] for i in range(0, len(dna_sequence), codon_size)]
        codon = codon_sequence[codon_number - 1]
        print(f"Codon: {codon}")
    amino_acid = amino_acid_codes.get(codon)
    print(f"Amino Acid: {amino_acid}")

print("Mutation Simulation-")
mutation = input("Choose a New Nucleotide Replacement(uppercase): ")
mutated_sequence = dna_sequence[:snp_position] + mutation + dna_sequence[snp_position+1:]
codon_size = 3
mcodon_sequence = [mutated_sequence[i:i + codon_size] for i in range(0, len(mutated_sequence), codon_size)]
mcodon_number = (snp_position // 3) + 1
print(f"Mutated Codon Number: {mcodon_number}")
mcodon = mcodon_sequence[mcodon_number - 1]
print(f"Mutated Codon: {mcodon}")
mamino_acid = amino_acid_codes.get(mcodon)
print(f"Mutated Amino Acid: {mamino_acid}")

print()

# rpoB gene analysis for MDR
rpoB_analysis = input("Would you like to perform rpoB gene analysis? (yes/no): ")
if rpoB_analysis == "yes":
    rpoB_gene = dna_sequence[759806:763325]
    rpoB_codon_sequence = [rpoB_gene[i:i + codon_size] for i in range(0, len(rpoB_gene), codon_size)]
    rpoB_codon_number = int(input("rpoB Codon Number: ")) -1
    rpoB_codon = rpoB_codon_sequence[rpoB_codon_number]
    print(rpoB_codon)
    amino_acid = amino_acid_codes.get(rpoB_codon)
    print(amino_acid)
    if 426 < rpoB_codon_number < 452:
        print("The entered codon cooresponds to an amino acid that is within the Rifampicin Resistance Determining Region (RRDR) of the rpoB gene. Mutations in this region may cause Rifampicin anibiotic-resistance strains in the bacterium.")
