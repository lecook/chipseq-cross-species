#! /usr/bin/env python

## Author: Laura E Cook, University of Melbourne

import sys
import os
import subprocess

input1 = sys.argv[1]
output = sys.argv[2]
reads =  sys.argv[3]    

cmd1 = 'frac=$( samtools idxstats {} | cut -f3 | awk -v ct={} \'BEGIN {{total=0}} {{total += $1}} END {{print ct/total}}\' )'
cmd2 = cmd1.format(
    input1,
    reads
)

cmd3 = 'samtools view -bs $frac {} > {}'
cmd4 = cmd3.format(
    input1,
    output
)
os.system(cmd2)
os.system(cmd4)
