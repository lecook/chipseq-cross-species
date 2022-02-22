#!/usr/bin/env bash

## Author: Laura E Cook, University of Melbourne
## Purpose: This script loops through all files in a directory
## These commands are then used in a slurm array script to run jobs in parallel

TRA=($(for file in *.maf; do echo $file |cut -d "." -f 1-2;done)) # change this to whatever your file prefix is

echo ${TRA[@]}

for tr in ${TRA[@]};

do

# generate lastz commands
echo 'lastz_32 /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.fa[multiple] /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1_RM/'${tr}.fa 'H=2000 K=2400 L=3000 Y=3400 --scores=/data/projects/punim0586/lecook/chipseq-pipeline/cross_species/bin/GenomeAlignmentTools/HoxD55.q --format=maf > /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/maf/'${tr}_mm10.smiCra1.maf

# generate convert maf to axt format commands
echo 'maf-convert psl '${tr}.maf' > '${tr}.psl

done