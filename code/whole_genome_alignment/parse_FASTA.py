#!usr/bin/env python3

import string
import random
import sys
import os
from Bio import SeqIO

## usage python3 parse_FASTA.py [input.fasta] > [output.fasta]

file = sys.argv[1]

# Loop through all the files in the variable
with open(file, 'r') as all_TWARs:

    # for each record (TWAR header) in the file parse it as a FASTA
    for record in SeqIO.parse(all_TWARs, "fasta"):
        # set the record ID to a variable
        id = record.id
        # set the TWAR sequence to a variable
        seq = record.seq

        print(">" + str(id) + "\n" + str(seq) + "*")
all_TWARs.close()
