#!/usr/bin/env bash

## This script loops through all scaffolds in the dunnart genome and
## generates a separate command for each scaffold
## These commands are then used in a slurm array script to run jobs in parallel

TRA=($(for file in *.maf; do echo $file |cut -d "." -f 1-2;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

# generate lastz commands
#echo 'lastz_32 /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.fa[multiple] /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1_RM/'${tr}.fa 'H=2000 K=2400 L=3000 Y=3400 --scores=/data/projects/punim0586/lecook/chipseq-pipeline/cross_species/bin/GenomeAlignmentTools/HoxD55.q --format=maf > /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/maf/'${tr}_mm10.smiCra1.maf



# generate convert maf to axt format commands
echo 'maf-convert psl '${tr}.maf' > '${tr}.psl

## build co-linear alignment chains
#echo 'axtChain -linearGap=loose -scoreScheme=../../bin/GenomeAlignmentTools/HoxD55.q '${tr}'.axt /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.2bit /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.2bit '${tr}'.chain'

#echo 'chainSort '${tr}'.chain > '${tr}'.sorted.chain'

done
