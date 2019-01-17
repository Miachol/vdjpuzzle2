#!/usr/bin/python -w
# Usage: python extract_fasta.py <NameOfSequenceWith> <FastaOfSequences>

import sys

if len(sys.argv) != 3:
	print "Incorrect number of arguments."
	print "Usage: python extract_fasta.py <NameOfSequenceWith> <FastaOfSequences>"
	sys.exit

input_name = sys.argv[1]
input_file = sys.argv[2]
with open(input_file) as f:
	all_names = []
	all_sequences = []
	for line in f:
		if line.startswith(">"):
			all_names.append(line.strip())
			all_sequences.append("")
		else:
			all_sequences[-1] += line.strip()
	everything = dict(zip(all_names, all_sequences))
	
with open("reference.fasta", "w") as f:
	f.write(input_name + '\n' + everything[input_name])

