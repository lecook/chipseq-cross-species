#!usr/bin/env python3

import string
import random
import sys
import os
from Bio import SeqIO


file = sys.argv[1]

# Loop through all the files in the variable
with open(file, 'r') as fasta:

    # for each record (TWAR header) in the file parse it as a FASTA
    for record in SeqIO.parse(fasta, "fasta"):
        # set each fasta to a seq
        seq = record.seq
        print("Count of A: ", seq.count("A"))
        print("Count of C: ", seq.count("C"))
        print("Count of G: ", seq.count("G"))
        print("Count of T: ", seq.count("T"))
fasta.close()

#GCATGC
