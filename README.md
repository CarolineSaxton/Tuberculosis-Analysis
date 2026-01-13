# Tuberculosis Gene Analysis
My project was to create a code to perform basic microbiology gene analysis, using python. The analysis involves amino acid conversions, nucleotide sequence strength determintation, and nucleotide mutation simulation. The knowldge that this information provides about the gene is important for scientists to be able to access easily for research on mutations and drug resistance in bacteria. It is currently tailored for analysis of the gene _M. Tuberculosis_ H37Rv, a common lab isolate of Mycobacterium Tuberculosis, the bacteria that causes tuberculosis infection. All referencing of the gene and its genome come from the NCBI database, accession number NC_000962.3. 

Tuberculosis is a serious infection that claims roughly a million lives globally each year, cause by the bacteria _Mycobacterium Tuberculosis_. While most cases can currently be treated with antibiotics, more and more antibiotic-resistant strains are appearing. This is a result of random mutations in the bacteria's DNA. Therefore, it is extremely impotant that as research on resistance continues to expand, basic genetic information is quickly acccessible and easy to understand, which can be done with code like this.


## Part 1: Counting and Graphing Total Nucleotide Counts
Knowing the total counts of each nucleotide and percent of each out of the total of any gene or sequence is important to understanding its strength. This is all determined by a gene's "GC Content", meaning its percent of Guanine or Cytosine nucleotides.  Guanine and Cytosine form triple hydrogen bonds between each other, stronger than those formed by the two other nucloetides (Adenine/Thymine). Therefore, the higher the GC content of a gene, the stronger its DNA is considered because it contains a greater number of (stronger) triple bonds.

1. Imports a sequence counter, adding up the totals of each A (Adenine), T (Thymine), C (Cytosine), and G (Guanine) nucleotide.
```from collections import Counter
nucleotide_counts = Counter(dna_sequence)
print(f"Total Count of Nucleotides: {nucleotide_counts}")
```

2. Create a bar graph to visualize the data from sequence counter for users.
```import matplotlib.pyplot as plt  # Import matplotlib for plotting
total_nucleotides = sum(nucleotide_counts.values())
nucleotide_percentages = {key: (value / total_nucleotides) * 100 for key, value in nucleotide_counts.items()}
plt.bar(  
    nucleotide_percentages.keys(),  # x-axis: Nucleotide types (A, T, C, G)
    nucleotide_percentages.values(),  # y-axis: Percentage of each nucleotide
    color=['yellow', 'red', 'blue', 'green']  # Colors for each nucleotide bar
)  
plt.title("Nucleotide Frequency")
plt.xlabel("Nucleotide")  # Indicates the nucleotides being counted (A, T, C, G)
plt.ylabel("Frequency")
plt.show()
```

Code output will look something like this:

<img width="630" height="474" alt="image" src="https://github.com/user-attachments/assets/1d1e974b-8358-4012-a37b-87ddb876d389" />


3. Determine the Guanine/Cytosine content
```
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100
gc_content = calculate_gc_content(dna_sequence)
print(f"GC Content: {gc_content:.2f}%")  # Output example: "GC Content: 46.67%"
```


## Part 2: Conversions to get amino acid/codon
Next, the code asks for either a nucleotide coordinate or codon number, used to determine the 3 bases in the codon and its cooresponding amino acid.
Allowing the user to enter either a nucleotide or codon allows for more seamless use, gien that both are very common in research and both are not always explicity given. Being able to easily determine the 3 bases in a codon or its amino acid is imporant to determine a mutation's affect on a sequence, simulated in the next part.

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


## Part 3: Simulating the affect of a mutation on the codon/amino acid
1. 

## Part 4: _rpoB_ gene analysis for Rifampicin resistance

1. Sets location of the rpoB gene in terms of the nucelotide range within the _Mycobacterium Tuberculosis_ H37Rv complete genome
```
rpoB_analysis = input("Would you like to perform 'rpoB gene' specific analysis? (yes/no): ")
if rpoB_analysis == "yes":
    rpoB_gene = dna_sequence[759806:763325]
```

2. Splits the rpoB DNA sequence into codons and identifies a specific codon within the rpoB sequence based on the entered number
```
    rpoB_codon_sequence = [rpoB_gene[i:i + codon_size] for i in range(0, len(rpoB_gene), codon_size)]
    rpoB_codon_number = int(input("rpoB Codon Number: ")) -1
    rpoB_codon = rpoB_codon_sequence[rpoB_codon_number]
    print(rpoB_codon)
```

3. Converts the identified codon into an amino acid and determines if it falls within the RRDR of the rpoB gene
```
    amino_acid = amino_acid_codes.get(rpoB_codon)
    print(amino_acid)
    if 426 < rpoB_codon_number < 452:
        print("The entered codon cooresponds to an amino acid that is within the Rifampicin Resistance Determining Region (RRDR) of the rpoB gene.")
```
