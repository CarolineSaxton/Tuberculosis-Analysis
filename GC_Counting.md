# Steps to Counting and Graphing Total Nucleotide Counts
__Knowing the total counts of each nucleotide and percent of each out of the total of any gene or sequence is important to its analysis. (Add more later)__

1. Imports a sequence counter, adding up the totals of each A (Adenine), T (Thymine), C (Cytosine), and G (Guanine) nucleotide.
> `from collections import Counter`  
`nucleotide_counts = Counter(dna_sequence)`  
`print(f"Total Count of Nucleotides: {nucleotide_counts}")`

2. Create a bar graph to visualize the data from sequence counter for users.
>`import matplotlib.pyplot as plt  # Import matplotlib for plotting`  
`total_nucleotides = sum(nucleotide_counts.values())`  
`nucleotide_percentages = {key: (value / total_nucleotides) * 100 for key, value in nucleotide_counts.items()}`  
`plt.bar(`  
    `nucleotide_percentages.keys(),  # x-axis: Nucleotide types (A, T, C, G)`  
    `nucleotide_percentages.values(),  # y-axis: Percentage of each nucleotide`  
    `color=['yellow', 'red', 'blue', 'green']  # Colors for each nucleotide bar`  
`)`  
`plt.title("Nucleotide Frequency")`  
`plt.xlabel("Nucleotide")  # Indicates the nucleotides being counted (A, T, C, G)`  
`plt.ylabel("Frequency")`  
`plt.show()`
Code output will look something like this:

![Image of Graph](file:///Users/csaxtonrowe1329/Desktop/Screenshot%202026-01-07%20at%202.54.30%E2%80%AFPM.png)

3. Determine the Guanine/Cytosine content
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

gc_content = calculate_gc_content(dna_sequence)

print(f"GC Content: {gc_content:.2f}%")  # Output example: "GC Content: 46.67%"
