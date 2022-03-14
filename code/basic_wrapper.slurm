#!/bin/bash

### Purpose: Basic SLURM wrapper script

### Author: Laura E Cook, University of Melbourne, 05/02/2020

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=10:00:00
#SBATCH --job-name=job
#SBATCH -A punim0586
#SBATCH -p mig
#SBATCH --error=slurmjob.stderr
#SBATCH --output=slurmjob.stdout
# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=lecook@student.unimelb.edu.au

conda create --name chip bowtie2 samtools picard deeptools multiqc phantompeakqualtools preseq macs2 bedtools fastqc pybedtools
