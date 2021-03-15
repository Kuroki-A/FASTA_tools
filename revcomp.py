
import sys

fasta_file = sys.argv[1]

name = ""
seq = ""

line_cnt = 0

with open(fasta_file) as f:
    for line in f:
        line = line.rstrip()
        line_cnt += 1

        if line_cnt == 1:
            name = line
        else:
            seq += line

#make revcomp seq
rev_comp_seq = seq.translate(
                 str.maketrans(
                    {'A':'T',
                     'T':'A',
                     'G':'C',
                     'C':'G',
                     'a':'t',
                     't':'a',
                     'g':'c',
                     'c':'g'
                     }
                 )
              )[::-1]

print(name)
print(rev_comp_seq)
