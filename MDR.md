# _rpoB_ gene analysis for Rifampicin resistance

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
