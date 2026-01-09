
from Bio import Entrez, SeqIO
Entrez.email = "add your email here" 

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
plt.title("Nucleotide Frequency")
plt.xlabel("Nucleotide")  # Indicates the nucleotides being counted (A, T, C, G)
plt.ylabel("Frequency")
plt.show()  # Renders the bar chart

# Calculate GC Content
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

gc_content = calculate_gc_content(dna_sequence)

print(f"GC Content: {gc_content:.2f}%")  # Output example: "GC Content: 46.67%"

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
        print("part of rrdr segment")

if gene_functions_analysis == "yes":
    import itertools
    virulence = {"ranges": itertools.chain(range(28362,29205), range(71587,71826), range(71821,72222), range(152324, 154129), range(154232,155599), range(161771,162673), range(196861,197658)), "Name": "virulence"}
    inform_path = {"ranges": itertools.chain(range(1,1524), range(2025,3260), range(3280,4437), range(5240,7267), range(7302,9818), range(12468,13016), range(43562,46471), range(58192,58482), range(58586,59080), range(59122,59374), range(59409,59867), range(60396,63020), range(79486,80193), range(278585,279529), range(499713,500366), range(661295,663373), range(734254,734970), range(734254,734970), range(735517,736224), range(748276,748812), range(748849,749241), range(759807,763325), range(763370,767320), range(769792,770550), range(781560,782404), range(782485,784590), range(784821,786011), range(800487,800792), range(800809,802435), range(802528,803370), range(803411,805106), range(805110,806166), range(811373,811741), range(811742,812622), range(812627,812812), range(812976,813374), range(813398,813937), range(813940,814308), range(814328,814990), range(814993,815630), range(827543,828076), range(828140,828892), range(1046136,1048415), range(1053765,1054241), range(1058260,1060575), range(1094886,1095059), range(1111612,1112223), range(1129152,1130105), range(1138967,1142671), range(1286595,1287020), range(1294168,1296054), range(1332092,1332964), range(1353522,1354136), range(1364413,1365186), range(1365344,1365808), range(1399970,1401661), range(1407339,1408238), range(1446379,1448031), range(1453204,1455012), range(1455163,1455405)), "Name": "information pathways"}
    cell_wall = {"ranges": itertools.chain(range(9914,10828), range(14089,14877),  range(53663,55699), range(55696,57378), range(68620,71559), range(80624,81673), range(81676,82668), range(98480,99250)), "Name": "cell wall"}
    rna = {"ranges": itertools.chain(range(10887,10960), range(11112,11184), range(25644,25726), range(80240,80440)), "Names": "rna"}
    inserts_phages = {"ranges": itertools.chain(range(33582,33794), range(400192,401703), range(472781,474106), range(475816,476184), range(606551,608062), range(701406,702014), range(702016,702759), range(832534,832848), range(889072,889398), range(889347,890333), range(890388,891482), range(921575,921865), range(947312,947644), range(1027104,1027685), range(1027685,1029337), range(1169423,1170670), range(1176928,1177242), range(1177239,1177373), range(1277893,1278300), range(1278269,1278820), range(1779314,1779724), range(1779930,1780241), range(1780199,1780699), range(1926202,1927137), range(1996152,1997413), range(2195989,2197353), range(2260665,2261688), range(2343027,2343916), range(2343994,2344224), range(2358389,2360041), range(2365465,2366726), range(2550065,2551326), range(2583435,2583779), range(2635628,2636889), range(2970551,2971549), range(2971659,2972027), range(2972160,2973421), range(2983071,2983874), range(3116818,3118227), range(3288464,32905606), range(3313283,3313672), range(3481451,3482698), range(3551281,3552542), range(3552764,3554025), range(3710433,3711694), range(3711749,3713461), range(3753765,3754256), range(3800092,3801463), range(3883525,3884917), range(3890830,3892091), range(4075752,4076099), range(4076484,4077730), range(4198874,4199089), range(4252993,4254327), range(4318775,4319266)), "Names": "inserts and phages"}
    met_resp = {"ranges": itertools.chain(range(14914,15612), range(34295,36610), range(65552,66694), range(66923,98362), range(75301,76212), range(83996,85168), range(89924,90403), range(90400,92322), range(92328,93278), range(93289,93951), range(93951,95417), range(95414,96892), range(97758,98351), range(99684,100451)), "Name": "met_resp"}
    regulatory = {"ranges": itertools.chain(range(27595,28365), range(86528,87133), range(89575,89919)), "Names": "regulatory"}
    conserved = {"ranges": itertools.chain(range(4434,4997), range(29245,29607), range(29722,31068), range(31189,31506), range(31514,31819), range(32057,33154), range(33224,33553), range(36867,37262), range(41304,41912), range(52831,53244), range(57410,57973), range(59896,60417), range(63200,63892), range(63909,64967), range(82478,83983), range(88204,89025), range(89022,89480)), "Names": "conserved"}
    lipids = {"ranges": itertools.chain(range(36607,36870), range(37259,38947), range(96927,97601), range(107600,108151), range(108156,109778), range(110001,117539), range(144049,145626), range(171215,172168), range(172211,173143), range(194993,196657), range(256064,257677), range(264067,265476), range(265507,266295), range(276058,277764), range(292171,293493), range(324567,326249), range(340998,341906), range(483977,485734), range(485731,489939), range(558895,559755), range(559888,560748), range(648536,649672), range(746363,747037), range(771484,773112), range(773123,774061), range(774783,775574), range(921970,922875), range(948559,949395), range(951632,952711), range(955077,956288), range(956293,958455), range(970505,972457), range(997782,999299), range(1008207,1008938), range(1131625,1133259), range(1180684,1182315), range(1221959,1222786), range(1264314,1264556), range(1264606,1264947), range(1270062,1271144), range(1313725,1315191), range(1315234,1319982), range(1320035,1321453), range(1335794,1337215), range(1349332,1351125), range(1349332,1351125), range(1485862,1487031), range(1508968,1509288), range(1509281,1510846), range(1510846,1512006), range(1517491,1518234), range(1599658,1601037), range(1659763,1660620), range(1673440,1674183), range(1674202,1675011), range(1682157,1686257), range(1712302,1714053)), "Names": "lipids"}

    print("This nucloetide is part of a gene that may affect the H37Rv mycobacterium's:")
    if snp_position in inform_path["ranges"]:
        print("Information Pathways")
    if snp_position in cell_wall["ranges"]:
        print("Cell Wall and Cell Processes")
    if snp_position in virulence["ranges"]:
        print("Virulence, Detoxification, and Adaption Genes")
    if snp_position in rna["ranges"]:
        print("Stable RNA's")
    if snp_position in inserts_phages["ranges"]:
        print("Insertion Sequences and Phages")
    if snp_position in met_resp["ranges"]:
        print("Intermediary Metabolism and Respiration")
    if snp_position in regulatory["ranges"]:
        print("Regulatory Proteins")
    if snp_position in conserved["ranges"]:
        print("Conserved Hypotheticals (unknown functions)")
    if snp_position in lipids["ranges"]:
        print("Lipid Metabolism")
