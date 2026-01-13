# Read Me for Tuberculosis Gene Analysis
This code can perform basic microbiology gene analysis involving amino acid conversions and nucleotide sequence strength determintation. It is currently tailored for analysis of the gene _M. Tuberculosis_ H37Rv, a common lab isolate of Mycobacterium Tuberculosis, the bacteria that causes tuberculosis infection. All refrencing of the gene and its genome come from the NCBI database, accession number NC_000962.3. 







t




## Part 1: Counting and Graphing Total Nucleotide Counts
Knowing the total counts of each nucleotide and percent of each out of the total of any gene or sequence is important to its analysis. (Add more later)

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

