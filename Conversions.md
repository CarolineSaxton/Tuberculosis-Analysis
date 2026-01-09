# Conversions to get amino acid/codon
Uses a nucleotide coordinate or codon number to determine the 3 bases in the codon and its cooresponding amino acid.

1. Codon to Amino Acid Conversion Chart
```
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
```

2. Nucleotide Coordinate or Codon Number entered, code converts sequence into groups of 3 (codons), then into amino acids
```
print("Performing Broad Genomic Analysis...")
choice1 = input("Would you like to use a nucleotide coordinate (input:'n') or codon number (input:'c')?  ")
if choice1 == "n":
    coordinate = int(input("Nucleotide Coordinate: ")) - 1
    print(f"Nucleotide: {dna_sequence[coordinate]}")
    codon_size = 3
    codon_sequence = [dna_sequence[i:i + codon_size] for i in range(0, len(dna_sequence), codon_size)]
    codon_number = (coordinate // 3) + 1
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
```
