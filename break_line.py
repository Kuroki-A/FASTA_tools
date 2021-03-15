
import sys

fasta_file = sys.argv[1]

name = "" #fasta header
seq = "" #seqence

line_cnt = 0

#Reading fasta file
with open(fasta_file) as f:
    for line in f:
        line = line.rstrip()
        line_cnt += 1

        if line_cnt == 1:
            name = line
        else:
            seq += line

#output result
print(name)
print(seq)
