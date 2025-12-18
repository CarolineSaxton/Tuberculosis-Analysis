'''from Bio import SeqIO
file_path = "/Users/csaxtonrowe1329/Downloads/OCA2_datasets/ncbi_dataset/data/gene.fna"
file_format = "fasta"
for seq_record in SeqIO.parse(file_path, file_format):
        # Access the sequence and other information
        sequence_id = seq_record.id
        dna_sequence = seq_record.seq  # This is a Bio.Seq.Seq object
        description = seq_record.description

        print(f"ID: {sequence_id}")
        print(f"Sequence: {dna_sequence}")
        print(f"Description: {description}")
        print("-" * 20)
        print(f"Sequence Length: {len(dna_sequence)}")  # Display the total length of the sequence
        print(dna_sequence[:100])  # Display the first 100 bases of the sequence
from collections import Counter  # Import Counter for easy counting of elements in a sequence

# Use Counter to count each nucleotide in the DNA sequence
nucleotide_counts = Counter(dna_sequence)

# Print the counts of each nucleotide
print(nucleotide_counts)

# Function to calculate the GC content of a DNA sequence
def calculate_gc_content(sequence):
    # Count the number of 'G' and 'C' bases in the sequence
    gc_count = sequence.count('G') + sequence.count('C')
    
    # Calculate GC content as a percentage of the total sequence length
    return (gc_count / len(sequence)) * 100

# Call the function and calculate GC content of OCA2 genome
gc_content = calculate_gc_content(dna_sequence)

# Print the GC content formatted to two decimal places
print(f"GC Content: {gc_content:.2f}%")  # Output example: "GC Content: 46.67%"'''

sequence = "AAACCCGGGTTT"
group_size = 3
codon_sequence = [sequence[i:i + group_size] for i in range(0, len(sequence), group_size)]
inputted = int(input("Codon Number: ")) -1
codon = codon_sequence[inputted]
print(codon)

dna_codons = {
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
amino_acid = dna_codons.get(codon)
print(amino_acid)