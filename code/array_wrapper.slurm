#!/bin/bash

### Author: Laura E Cook, University of Melbourne, 22/01/2020
### Last update: 07/02/2020
### Purpose: basic array wrapper script

# Array set up:
#SBATCH --array=1-10

# Partition for the job:

# The project ID which this job should run under:
#SBATCH -A punim0586

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --mem=50000

#SBATCH --partition mig

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# ends successfully
#SBATCH --mail-type=END

#SBATCH --mail-user=lecook@student.unimelb.edu.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=100:00:00

# The name of the job:
#SBATCH --job-name=array

# Output control:
#SBATCH --error="/data/projects/punim0586/lecook/chipseq-pipeline/cross_species/logs/mafToAxt_%a.stderr"
#SBATCH --output="/data/projects/punim0586/lecook/chipseq-pipeline/cross_species/logs/mafToAxt_%a.stdout"

# Load modules:

module load foss
module load lastz
#conda activate wga
module load last/1066
module load python/2.7.16
module load ucsc/21072020
PATH=$PATH:/data/projects/punim0586/lecook/chipseq-pipeline/cross_species/bin/GenomeAlignmentTools/src
PATH=$PATH:/data/projects/punim0586/lecook/chipseq-pipeline/cross_species/bin/GenomeAlignmentTools/kent/bin;export PATH

export KENTSRC_DIR=/data/projects/punim0586/lecook/chipseq-pipeline/cross_species/bin/GenomeAlignmentTools/kent/src/

# Run the simulations:
command=`head -n $SLURM_ARRAY_TASK_ID /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/scripts/jobList | tail -n 1`
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo $command
eval $command
