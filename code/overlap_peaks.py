#! /usr/bin/env python

## Author: Laura E Cook, University of Melbourne
## Purpose:
#   1. Import replicate 1, replicate 2, and output file name from snakemake rule
#   2. Concatonate peak 1 adn peak 2 to get a pooled_peaks file
#   3. Intersect pooled peaks with peak 1
#   4. Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported.
#   5. Calculate
#   5.

import sys
import os
import subprocess

# replicate 1
peak1 = sys.argv[1]

# replicate 2
peak2 = sys.argv[2]

# pooled peaks
pooled = sys.argv[3]

output = sys.argv[4]

awk_param = '{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'
cut_param = '1-10'

cmd1 = 'intersectBed -wo '
cmd1 += '-a {} -b {} | '
cmd1 += 'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {}\' | '
cmd1 += 'cut -f {} | sort | uniq | '
cmd1 += 'intersectBed -wo '
cmd1 += '-a stdin -b {} | '
cmd1 += 'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {}\' | '
cmd1 += 'cut -f {} | sort | uniq > {}'
cmd2 = cmd1.format(
    pooled,  # peak_pooled
    peak1,  # peak1
    awk_param,
    cut_param,
    peak2,  # peak2
    awk_param,
    cut_param,
    output)

os.system(cmd2)
