#!/usr/bin/python -w
# Usage: python mpileup2cons.py reference.fasta
# mpileup result is read from stdin

import csv
import sys

# reference.fasta
original_ref = sys.argv[1]

# consider allowing standard input or reading from a file, or at least require '-' to be passed as an argument to denote standard input when calling the script
with sys.stdin as f:
    ## Parse the output of mpileup and bcftools view -cg to obtain positions 
    ## where majority nucleotide != original reference.
    all_changes = []
    for line in csv.reader(f,delimiter='\t'):
        if line[0][0] == "#": # skip comments
            continue
        elif line[4] not in ["A","C","G","T"]: # skip if not recognised
            continue 
        else:
            # Determine frequency of variant
            info = line[7].split(";")
            # frequencies is [ref-forward , ref-reverse , alt-forward , alt-reverse]
            frequencies = map(int, filter(lambda x: 'DP4=' in x, info)[0][4:].split(','))
            ref_freq = frequencies[0] + frequencies[1] 
            var_freq = frequencies[2] + frequencies[3]
            
            if var_freq > ref_freq:
                change = [line[1],line[3],line[4]] # position,reference,variant
                all_changes.append(change)

# obtain the original reference sequence from the fasta file provided
with open(original_ref) as f:
    reference_data = f.read().rstrip().split('\n')
    assert reference_data[0].startswith('>'), "The reference file is not in fasta format (does not start with '>')"
    seq_name = reference_data[0]
    original_seq = list("".join(reference_data[1:])) # as a single string then as a list of single nucleotides
    
new_seq = original_seq[:]
print("Processing sequence {0}".format(seq_name))
for entry in all_changes:
    pos,old_nuc,new_nuc = entry
    assert original_seq[int(pos)-1] == old_nuc, "The positions of the required changes are not aligned with the original reference sequence."
    print("Changing position {0} from nucleotide {1} to {2}.").format(pos,old_nuc,new_nuc)
    new_seq[int(pos)-1] = new_nuc # replace nucleotides with changes identified from mpileup and bcftools

with open("new_reference.fasta","w") as f:
    f.write(seq_name + "\n" + "".join(new_seq) + "\n")
