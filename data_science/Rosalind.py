seq_rev = "ATCTAATCTCTCCCCAGGAGAGAATGTCGGCATCGGGGGAGGACCGGTCTAACTGCGAGATCTACACGTGCAATTTTGAGAGGTAGTCAGCACCTCGAGAATAAACGGACAACAGTGTCATTGGGCATTTGCAGACGAAACCTGGTGGGTTTAGCTAACAACATATTTAATATGAGTGTATGGCCAATAAGGAATTGTGGCGGTCGACTAGGCATGTTGCTTTCCGCATTGTCGTTAAGTATTTCTGATGCAAAGGGAAGGGTCTAAGGGTCTTTAGGAGACAACCGGTTCGCTATAGCAAGTGCGCTCAGCGCTGATCAGACAGCAGCATGATGGATATTCGGCCGGCTTCTGATAGCGGGAGCGTCTTTGCAACCACCCTCTCAAGTTCTCTTCCGCCGTCGTTCTCGTGCTACACCCAGAAGGTGGACGGCAAACCTCCCTTTTCACACTACTTGCGAGACCGTTGAAAAACGCTAGAAATCATCCCTATGGACCTTACGCTAATGTGGGTTCCTGAGCCCAAGAGAAACCCTTCGTGACCTAAATACCGCAGGAGCTACACTAACACATCAAATACGTCTCGGCCGTGAATCTAAATTGAATATTGTCGGAATCGTGTCTGATTACCTATTGGATGCGACCTCTAGTAGTCCGCCTGCCCTGCTTTGACCTCGGTCGCTTGCACCCCGTTGGATTGCTCAGCTGGTCCGTGCGCACTGGTGTCCAAAGAGCCAGATGTTAAGAAGGGCTGGAGTGCAATCACATCGGCATGAAAGTTTCTCGAATAGCGCCTCAGTTGCGGAAAGCGTCATATGACTGATGTACTTTCTGCTTTTCTACAGCAGCAATGATAGTCTTAGACGATGTACTCAGACGTATGAACGCACCA"

def reverse_complement(seq):
    nucleotides = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(nucleotides[base] for base in reversed(seq))

print(reverse_complement(seq_rev))

nucleotides = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
seq_rev_fin = ''.join(nucleotides[base] for base in reversed(seq_rev))
print("\nThe final sequence is:") 
print(seq_rev_fin)


DD = 1
Dd = 0.5
dd = 0

v_1=30
v_2=21
v_3=29
v_tot = v_1+v_2+v_3
value_1 = v_1/v_tot
value_2 = v_2/v_tot
value_3 = v_3/v_tot
P = value_1 * DD + value_2 * Dd + value_3 * dd
P

def dominant_phenotype_prob(k, m, n):
    total = k + m + n
    # total_pairs = total * (total - 1)
    prob = 0
    prob += (k / total) * ((k - 1) / (total - 1)) * 1        # AA x AA
    prob += (k / total) * (m / (total - 1)) * 1 * 2          # AA x Aa and Aa x AA
    prob += (k / total) * (n / (total - 1)) * 1 * 2          # AA x aa and aa x AA
    prob += (m / total) * ((m - 1) / (total - 1)) * 0.75     # Aa x Aa
    prob += (m / total) * (n / (total - 1)) * 0.5 * 2        # Aa x aa and aa x Aa
    return round(prob, 5)

print(dominant_phenotype_prob(18, 27, 16))  # Output: 0.78333

"https://chatgpt.com/c/68079a1b-d5a4-8010-94f4-2543e568d6ef"

val_1, val_2, val_3, val_4, val_5, val_6 = 16703, 19371, 16863, 19248, 16919, 16659
output = 2*(val_1 * 1 + val_2 * 1 + val_3 * 1 + val_4 * 0.75 + val_5 * 0.5 + val_6 * 0)
output

seq1= "GAACTTATCGTGTATGCGACCAACGTATTGACCGTACTGAAATTCTGAAACGGGACCGTGGTCCCCGTATTAGCTGGTACAGACCTTACCTCCAGCCTTTCCACCGATTGGTCCAGCGTTAGGCCACAGATGTTCCAGTCCTGTCGGCGAACCTCTCAACGAGAAATTAACTGAGGCCGCGAATCCTGTTACCTTACTCCAGAGACTCGTTTTGTGGCGGGGGGAGAAACATTCGCCTTTACGCGCAGACATACGATTCGAGGGTGGCTAGGTTCGGAGAAATGTAAAGTTTATTGATTCGATTGGGACAATCACAGATTAGACGCCATCGAGTCGATTATGGCGTTACTTCCATAGAAACCTCAGCACCCCCACCTGTGCTACCCGGCGGACCCTCTGGGAAAGCTTGTCCATTCCTCTGTTTTCCCTCCCCTACTCCTCGGATGTCCCAGCTAGAATGCGCGCCACCTGGCTAACCGCGATAGATGGGATGCTGCGCTCTGACCTCGACGCGGCTTATACACGTACAGGTGCAATGCCTTTACCGATCGTAGCGTTTCCGAATGACATAAATTACTTCTGAGGCCTAAATGTTTAGAGTTTAGACCCGCCTTTTACATTTTAATTCCAGGCCAAGACCTCCAACGCTTTCTTCCCTTATGTATGGGGCTCAGGGATGCAAGGGTACCACGACATACCCACCAGACAGTAGCGTGATCTAGGAGCAAATATGTCCAAGTGGCCAGTCGCTAACCGTGCCTAGCGCGCATAGTGCACCTTACGACACCGGACGATTCCATGGGGGACTACTGTTTCTCTGGATGTAACTTGGTATCGTGCGATACGCCAGTGGGATATCATACACCTGTTTTAGAGGGCCCGTTATCTCCTACAGACGATGTCCGGCATACGAGTACGAAAAATGACG"
seq2= "GGATTAAACGATATCGAGACAAGGCGCTTGACCGTTGTGCGAGGAGTATCCGGTGCGTTAGTCTACCAGCTTCCTTGGACAGTACCCGCCTCCGGCTGTCGTACCGACCGGTGTGCGTTTCGGAGCAATGATGATGATCCTAAGATGGGAACTCGTGAACAAGCTATTTACTGAGGCACCGATTAACGGTCCCAATATACGTTGAGGCGTTATGTGGTGGGATGATAAGCAGTCGGCATTCCGCACAGAAGTCACATTCAAGAATCTCAGGAGACGGTATAATGGAATCTATTCTGATTTACGTCAGCTTTTCGCAAACTAGTGGCGATCACCTAGAATATGACGTAAGTTATATCCAAGAGTGAACACGCTCACAGTTATGATCCCGGCCACAGGGAGGTTTTCCTTATTCGCTCCTCTCTTGTCCCAAGATAATTGTTCAGAGGCCGAAGGTAATGTATCCATCCATTGGCTCACTCCATCCCCTGAAATATTATGCGCCCTAGTGGAGGCCAGTGCTATGCTACTTGCACACAATACTTTCCTCATCCTGGACACTTTTAATGATATCATTTAGGTCTGTGGCCTATGCGCTTCGATATTTGAGCACCCAAGAACATATCATTTACACGTTACATACAGGGTTCGAATATCACCTCAGAATGCCGCCTGTGTTCTTCCATCATAACATAGGTTTGCCAACATCCAGAGGGGCGAACATTTGTTTCGTTTTTCATAGTGTTAAGTGGGGCACTTACCGCAGCGAGCGAACAGAACCCCACTGCACGTGACACCTCCTGGGGAGGTGACTTAATCCGGGGCAACCCCTGTAAATGCTGCGATCCCCAGGTTGGGTCCTCCACGCCCGTTCTCGAGCCATCGAGAACAGCTATGGAACATGAGATCATTACGAGAAAGACAAATTCAC"
mismatch = 0
for base1, base2 in zip(seq1, seq2):
    if base1 != base2:
        mismatch += 1

tot = sum(1 for base1, base2 in zip(seq1, seq2) if base1 != base2)

def hamming_distance(s, t):
    return sum(1 for a, b in zip(s, t) if a != b)

whole_string = "GATATATGCATATACTT"
whole_substring = "ATAT"

for i in whole_string:
    if  whole_string == i.index("ATAT"):
        print(f"Substring found at index: {i}")






whole_string = "GGATGGGGCCAAGACACTTCGGCTCCTTGAGACTTCGGAATAACTTCGGACTTCGGCCACTTCGGGGTGGACGACTTCGGGAACCACTTCGGACTTCGGACTTCGGACTTCGGGCGACTTCGGAACACTTCGGAGGAACTTCGGGACTTCGGAGACTTCGGTACTTCGGACTTCGGACTTCGGGAAGGACTTCGGTCACTTCGGATAACTTCGGGATGAACTTCGGTGAAGACTTCGGACTTCGGGTACCGTGACTTCGGCTTGTTTACGTACTTCGGAAACTTCGGTTCGCGGACTTCGGCCCCACTTCGGTTGCTCACCGGGACTTCGGGACTTCGGACTTCGGCTAGACTTCGGTCCCATATTGACTTCGGTGGGACTTCGGGTGCCTACTTCGGGGGAGATCGACTTCGGATAACTTCGGGTACTTCGGAACTTCGGATAACTTCGGACTTCGGCAACTTCGGATACATACTTCGGCCTAAACTTCGGGACTTCGGATGACTTCGGCGCGCATGACTTCGGGTGACTTCGGTACTTCGGATACTTCGGCCACTTCGGCCTCACTTCGGTACTTCGGATATTACTTCGGAACTTCGGCCACTTCGGTTGGTACTTCGGGACTTCGGAATACTTCGGCGTAACTTCGGGATACGCACTTCGGTAGACTTCGGACTTCGGACTTCGGCACTTCGGACTTCGGAACTTCGGGCACTTCGGCACATGATGATTACTTCGGGAGCACACTTCGGTAACTTCGGTACTTCGGTCACTTCGGAACTTCGGACTTCGGTAACTTCGGCGGAGACTTCGGGAACTTCGGGACTTCGGAACTTCGGGAATTAGACTTCGGAAAGGTACTTCGGACTTCGGACTTCGGAGTCCTCTGCGACTTCGGACTTCGGGTGTAAATACTTCGGACTTCGGAAAGTCTTACTTCGG"
whole_substring = "ACTTCGGAC"

for i in range(len(whole_string) - len(whole_substring) + 1):
    if whole_string[i:i+len(whole_substring)] == whole_substring:
        print(i + 1)


seq_revp = "ATATCCTAAGTAGACTCATTCTTCACACGGCGAAGACAGAGTGAATGAAACCATATTGCCAAAGTACTGATAGCTGCCAAAGAGTCTTAATGACAGCTGGTAATTCAGCGTCATCCTACACTGCTGCGTCGCTAGAAGAGGGTAGAAAGATCGTCATTGTGGATCTCCAATGGGCGTTTAGGTTTAGGTGTTCCTTGTAAGGGACGGTCATGATTCGCTTTGGGAAATCGGCGACTAAACAATGCTACAATTCGGAGCGCACACGCGCCATGTGTGTATACGGTAGTCCTGATCTTCGTCGAGGTTTCCTCAGAGATTTGTATCATTTCGTGCAAATGGCCGTCGCAAGTAACATAGTAAGTGGCGGGCTAAAACATTGATCCTGACATACGACTCTGCGGGTCCGCACCCTCTTGCACGATCGCTCCGGTCAAGCCCTGCAGGGCTTTCGGAGCCTCCGCGGGGTATTCCGGCTTATTGTTATTCAACTTCTCATCGGTATTTCCCTTCACACCAAAGTCAGCTGCAACTCTCAAAACAGGTCCATTCTTGGAGATTCGAGAACCCTTGTATGCTCTCAGACTCAAACACTCCGTGACTAATTAGTTACACGCGTCTACAACGCCAATAACGCTCCAATGGATTTTAACTTGTTAAGCTTATTAGGCGCAAGTAGTCTGCATCAGCGGGGATAGATGGGTCTATGACGGAGAAACGCCTCATTATATGACCGAAACCTCTGTCAACCAGGCGCCTATGACCCGGGCAAGATACATAGGTGTACTTGTGCTTATATCTAATGCGCGTCATCCGTTCTAAATCTTCAATGCTACGACTTATAGGCCCCTTACATACAGTATTCAGAATTGAGCGATGCGCG"
bases= {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

for i in range(len(seq_revp)):
    # We loop through every starting index i in the sequence.
    # We’re going to try extracting substrings starting from each of these positions.
    for l in range(4, 12 + 1):
    # For each starting index i, we try substrings of different lengths, from 4 to 12.
        if i+l <= len(seq_revp):
        # This makes sure the substring doesn’t go past the end of the sequence.
        # If i = 10 and l = 5, we want seq[10:15].
        # But if i + l > len(seq), it would go out of bounds — so we skip those.
            sub = seq_revp[i:i+l]
            # Here we extract the substring of length l starting at index i.
            # For example, if i = 3 and l = 6, we get seq[3:9].
            rev_comp = ''.join(bases[b] for b in reversed(sub))
            # This line computes the reverse complement of the substring:
            # reversed(sub) gives us the reverse of the substring.
            # bases[b] finds the complement of each base.
            # ''.join(...) puts it back together as a string.
            if sub == rev_comp:
                print(i+1, l)

rev_comp = ''.join(bases[base] for base in reversed(sub))

#################################################################

with open("your_fasta_file.txt") as f:
    lines = f.readlines()

sequence = ''
for line in lines:
    if not line.startswith('>'):
        sequence += line.strip()

# Now `sequence` contains the full DNA string (IDs excluded)
print(sequence)

codons = {
"UUU" : "F",
"CUU" : "L",
"AUU" : "I",
"GUU" : "V",
"UUC" : "F",
"CUC" : "L",
"AUC" : "I",
"GUC" : "V",
"UUA" : "L",
"CUA" : "L",
"AUA" : "I",     
"GUA" : "V",
"UUG" : "L",
"CUG" : "L",
"AUG" : "M",
"GUG" : "V",
"UCU" : "S",
"CCU" : "P",
"ACU" : "T",
"GCU" : "A",
"UCC" : "S",
"CCC" : "P",
"ACC" : "T",
"GCC" : "A",
"UCA" : "S",
"CCA" : "P",
"ACA" : "T",
"GCA" : "A",
"UCG" : "S",
"CCG" : "P",
"ACG" : "T",
"GCG" : "A",
"UAU" : "Y",
"CAU" : "H",
"AAU" : "N",
"GAU" : "D",
"UAC" : "Y",
"CAC" : "H",
"AAC" : "N",
"GAC" : "D",
"UAA" : "Stop",
"CAA" : "Q",
"AAA" : "K",
"GAA" : "E",
"UAG" : "Stop",
"CAG" : "Q",
"AAG" : "K",
"GAG" : "E",
"UGU" : "C",
"CGU" : "R",
"AGU" : "S",
"GGU" : "G",
"UGC" : "C",
"CGC" : "R",
"AGC" : "S",
"GGC" : "G",
"UGA" : "Stop",
"CGA" : "R",
"AGA" : "R",
"GGA" : "G",
"UGG" : "W",
"CGG" : "R",
"AGG" : "R",
"GGG" : "G"
}

with open("rosalind_prot.txt") as f:
    rna_to_prot = f.read().strip()
    
prot = ""
for i in range(0, len(rna_to_prot) -2, 3):
    # So to avoid going out of bounds, you stop at length - 2
    new_cod = rna_to_prot[i:i+3]
    amino = codons[new_cod]
    # new_cod = "AUG"
    # amino = codons[new_cod]   
    # codons["AUG"] returns "M"
    if amino == "Stop":
        break
    prot += amino

########################################################
lines = fasta.strip().split('\n')                      
dna = ''
for line in lines:
    if not line.startswith('>'):
        dna += line.strip()
########################################################

########################################################
fasta_dict = {}
current_id = ""

for line in fasta.strip().split('\n'):
    if line.startswith('>'):
        current_id = line[1:]  # remove '>'
        fasta_dict[current_id] = ""
    else:
        fasta_dict[current_id] += line.strip()
########################################################

with open("rosalind_splc1.txt") as f:
    content = f.read().strip().split(">")

# Remove any empty strings from splitting   
entries = [e for e in content if e]
dna_splicing = []
for entry in entries:
    lines = entry.splitlines()
    sequence = "".join(lines[1:])
    dna_splicing.append(sequence)

dna_seq = dna_splicing[0]
introns = dna_splicing[1:]

for intron in introns:
    dna_seq = dna_seq.replace(intron, "")

rna_seq = dna_seq.replace("T", "U")
    # So to avoid going out of bounds, you stop at length - 2
prot_splice=""
for i in range(0, len(rna_seq) - 2, 3):
    cod_splice = rna_seq[i:i+3]
    amino = codons[cod_splice]
    if amino == "Stop":
        break
    prot_splice += amino

prot_splice

################################################################################

from itertools import permutations

n = 5

perms = list(permutations(range(1, n+1)))

"""
Access elements by index

Use len() to count them

Reuse the data multiple times without re-generating it

If you don’t use list(), you'd only be able to loop through the permutations once.

Recap:
permutations(...) → lazy generator of all permutations

list(...) → makes it a real list you can store, count, loop multiple times
"""

# Print the total number of permutations
print(len(perms))

# Print each permutation
for p in perms:
    print(" ".join(map(str, p)))
    # This converts each integer in the tuple to a string, because join() only works with strings.
    # This line is printing one permutation (which is a tuple of integers) as a space-separated string. 
    # map() applies a function to every item in an iterable (like a list or tuple).

#################################################
""" map(str, iterable)
map() applies a function to every item in an iterable (like a list or tuple).

str is a function that turns things into strings.

So:


map(str, [1, 2, 3])
means:

“Take each number in the list and turn it into a string.”

Result:

['1', '2', '3']
✅ It's just like doing:

[str(1), str(2), str(3)]
"""
#######################################
""" ".join(list_of_strings)
" ".join(...) takes a list of strings and joins them into one single string, with spaces in between.

# Example:

" ".join(['1', '2', '3'])
Result:

'1 2 3'
"""

##################################
"""
p = (1, 2, 3)


map(str, p)      →  ['1', '2', '3']
" ".join(...)    →  '1 2 3'


print(" ".join(map(str, p)))   →  prints: 1 2 3
"""

#####################################################################
from itertools import permutations, product

n = 5

perms = list(permutations(range(1, n+1)))
# permutations(nums) gives all orderings of [1, 2] → [(1, 2), (2, 1)]
# Print the total number of permutations

# Print each permutation
signed_perms=[]
for p in perms:
    for sign in product([-1, 1], repeat=n):
    # product([-1, 1], repeat=n) gives all sign combos → [(-1, -1), (-1, 1), ...]
        s_perm = [a * s for a, s in zip(p, sign)]
    # for a, s in zip(...) This loops through each pair. At each loop: a is one number from permssign is the matching s from sign
    # zip(perms, sign) This pairs each element of perms with a corresponding element in sign.
    # This pairs elements from two lists together.
        signed_perms.append(s_perm)

print(len(signed_perms))
for p in signed_perms:
    print(" ".join(map(str, p)))

#####################################################################
n = int(open("./Data/rosalind_SIGN.txt", "r").read())

lst = [x for x in permutations([i for i in range(-n,n+1) if (i != 0)], n)]
# Give me all permutations of length n from the list of numbers from -n to n (excluding 0)
print(len(lst))
for i in lst:
    print(" ".join(map(str,i)))
#####################################################################

from itertools import permutations

list_input = [-2, -1, 1, 2]
n = 3

# permutations(list_input, 2)
for p in permutations(list_input, n):
    print(p)
