##!/usr/bin/env bash

## This script loops through all scaffolds in the dunnart genome and
## generates a separate RepeatMasker command for each scaffold
## These commands are then used in a slurm array script to run jobs in parallel

TRA=($(for file in *.fa; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo 'RepeatMasker -q -xsmall smiCra1/'${tr}'.fa -default_search_engine hmmer -trf_prgm /home/lecook/.conda/envs/wga/bin/trf -hmmer_dir /home/lecook/.conda/envs/wga/bin/'


done
