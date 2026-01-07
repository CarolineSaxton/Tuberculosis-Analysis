# Counting and Graphing Total Nucleotide Counts
Knowing the total counts of each nucleotide and percent of each out of the total of any gene or sequence is important to its analysis. 

`from collections import Counter
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
plt.title("Nucleotide Frequency")
plt.xlabel("Nucleotide")  # Indicates the nucleotides being counted (A, T, C, G)
plt.ylabel("Frequency")  # Shows the frequency of each nucleotide
plt.show()  # Renders the bar chart

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

gc_content = calculate_gc_content(dna_sequence)

print(f"GC Content: {gc_content:.2f}%")  # Output example: "GC Content: 46.67%"

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
        print("part of rrdr segment")`
