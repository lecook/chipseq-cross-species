#! /usr/bin/env python

## Author: Laura E Cook, University of Melbourne

import sys
import os
import subprocess

input1 = sys.argv[1]
output = sys.argv[2]
        
reads=10000000     

cmd1 = 'frac=$(samtools idxstats {} | cut -f3 | awk -v ct={} \'BEGIN {{total=0}} {{total += $1}} END {{print ct/total}}')'
cmd1 += 'samtools view -bs $frac {} > {}'
cmd2 = cmd1.format(
    input1,  # peak_pooled
    reads,  # peak1
    input1,
    output)

os.system(cmd2)
