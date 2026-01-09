# rpoB gene analysis and general gene function analysis

1. Rifampicin Resistance Gene
```
rpoB_analysis = input("Would you like to perform 'rpoB gene' specific analysis? (yes/no): ")
if rpoB_analysis == "yes":
    rpoB_gene = dna_sequence[759806:763325]
    rpoB_codon_sequence = [rpoB_gene[i:i + codon_size] for i in range(0, len(rpoB_gene), codon_size)]
    rpoB_codon_number = int(input("rpoB Codon Number: ")) -1
    rpoB_codon = rpoB_codon_sequence[rpoB_codon_number]
    print(rpoB_codon)
    amino_acid = amino_acid_codes.get(rpoB_codon)
    print(amino_acid)
    if 426 < rpoB_codon_number < 452:
        print("part of rrdr segment")
```
