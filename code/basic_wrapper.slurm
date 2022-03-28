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

#conda activate chip
#multiBamSummary bins -p 32 -b ENCFF002HZV_E14.5_H3K4me3_q30.sorted.dedup.bam ENCFF045IPK_E10.5_H3K4me3_q30.sorted.dedup.bam ENCFF058AUT_E12.5_H3K4me3_q30.sorted.dedup.bam ENCFF117QRC_E13.5_H3K4me3_q30.sorted.dedup.bam ENCFF124TAB_E13.5_H3K4me3_q30.sorted.dedup.bam ENCFF124UYX_E10.5_H3K4me3_q30.sorted.dedup.bam ENCFF157KEH_E10.5_H3K4me3_q30.sorted.dedup.bam ENCFF182XFG_E15.5_H3K4me3_q30.sorted.dedup.bam ENCFF182ZPF_E12.5_H3K4me3_q30.sorted.dedup.bam ENCFF184CUE_E11.5_H3K4me3_q30.sorted.dedup.bam ENCFF203JQV_E12.5_H3K4me3_q30.sorted.dedup.bam ENCFF248PGK_E13.5_H3K4me3_q30.sorted.dedup.bam ENCFF258KCR_E15.5_H3K4me3_q30.sorted.dedup.bam ENCFF376FGM_E11.5_H3K4me3_q30.sorted.dedup.bam ENCFF401BKM_E15.5_H3K4me3_q30.sorted.dedup.bam ENCFF485UDC_E13.5_H3K4me3_q30.sorted.dedup.bam ENCFF665QBJ_E14.5_H3K4me3_q30.sorted.dedup.bam ENCFF717QDV_E11.5_H3K4me3_q30.sorted.dedup.bam ENCFF724DMU_E14.5_H3K4me3_q30.sorted.dedup.bam ENCFF727QTS_E15.5_H3K4me3_q30.sorted.dedup.bam ENCFF760QYZ_E11.5_H3K4me3_q30.sorted.dedup.bam ENCFF784ORI_E14.5_H3K4me3_q30.sorted.dedup.bam ENCFF825AVI_E10.5_H3K4me3_q30.sorted.dedup.bam ENCFF941QJZ_E12.5_H3K4me3_q30.sorted.dedup.bam --centerReads -out H3K4me3_multibamsum_mouse.npz --outRawCounts H3K4me3_multibamsum_mouse.txt

#multiBamSummary bins -p 32 -b ENCFF002HZV_E14.5_H3K27ac_q30.sorted.dedup.bam ENCFF011NFM_E12.5_H3K27ac_q30.sorted.dedup.bam ENCFF058AUT_E12.5_H3K27ac_q30.sorted.dedup.bam ENCFF117QRC_E13.5_H3K27ac_q30.sorted.dedup.bam ENCFF157KEH_E10.5_H3K27ac_q30.sorted.dedup.bam ENCFF182XFG_E15.5_H3K27ac_q30.sorted.dedup.bam ENCFF184CUE_E11.5_H3K27ac_q30.sorted.dedup.bam ENCFF194ORC_E13.5_H3K27ac_q30.sorted.dedup.bam ENCFF203JQV_E12.5_H3K27ac_q30.sorted.dedup.bam ENCFF213EBC_E10.5_H3K27ac_q30.sorted.dedup.bam ENCFF248PGK_E13.5_H3K27ac_q30.sorted.dedup.bam ENCFF290ZNF_E13.5_H3K27ac_q30.sorted.dedup.bam ENCFF327VAO_E14.5_H3K27ac_q30.sorted.dedup.bam ENCFF376FGM_E11.5_H3K27ac_q30.sorted.dedup.bam ENCFF394TZN_E12.5_H3K27ac_q30.sorted.dedup.bam ENCFF512SFE_E11.5_H3K27ac_q30.sorted.dedup.bam ENCFF515PKL_E11.5_H3K27ac_q30.sorted.dedup.bam ENCFF548BRR_E10.5_H3K27ac_q30.sorted.dedup.bam ENCFF584JFB_E15.5_H3K27ac_q30.sorted.dedup.bam ENCFF707WKL_E15.5_H3K27ac_q30.sorted.dedup.bam ENCFF727QTS_E15.5_H3K27ac_q30.sorted.dedup.bam ENCFF784ORI_E14.5_H3K27ac_q30.sorted.dedup.bam ENCFF825AVI_E10.5_H3K27ac_q30.sorted.dedup.bam ENCFF902HAR_E14.5_H3K27ac_q30.sorted.dedup.bam --centerReads -out H3K27ac_multibamsum_mouse.npz --outRawCounts H3K27ac_multibamsum_mouse.txt
#plotCorrelation --corData H3K4me3_multibamsum_mouse.npz --plotFile H3K4me3_pearsoncor_multibamsum_mouse.pdf --outFileCorMatrix H3K4me3_pearsoncor_multibamsum_matrix_mouse.txt --corMethod pearson --whatToPlot heatmap --skipZeros --plotNumbers --colorMap RdYlBu
#plotCorrelation --corData H3K27ac_multibamsum_mouse.npz --plotFile H3K27ac_pearsoncor_multibamsum_mouse.pdf --outFileCorMatrix H3K27ac_pearsoncor_multibamsum_matrix_mouse.txt --corMethod pearson --whatToPlot heatmap --skipZeros --plotNumbers --colorMap RdYlBu

plotFingerprint -p 32 -b ENCFF002HZV_E14.5_H3K27ac_q30.sorted.dedup.bam ENCFF011NFM_E12.5_H3K27ac_q30.sorted.dedup.bam ENCFF058AUT_E12.5_H3K27ac_q30.sorted.dedup.bam ENCFF117QRC_E13.5_H3K27ac_q30.sorted.dedup.bam ENCFF157KEH_E10.5_H3K27ac_q30.sorted.dedup.bam ENCFF182XFG_E15.5_H3K27ac_q30.sorted.dedup.bam ENCFF184CUE_E11.5_H3K27ac_q30.sorted.dedup.bam ENCFF194ORC_E13.5_H3K27ac_q30.sorted.dedup.bam ENCFF203JQV_E12.5_H3K27ac_q30.sorted.dedup.bam ENCFF213EBC_E10.5_H3K27ac_q30.sorted.dedup.bam ENCFF248PGK_E13.5_H3K27ac_q30.sorted.dedup.bam ENCFF290ZNF_E13.5_H3K27ac_q30.sorted.dedup.bam ENCFF327VAO_E14.5_H3K27ac_q30.sorted.dedup.bam ENCFF376FGM_E11.5_H3K27ac_q30.sorted.dedup.bam ENCFF394TZN_E12.5_H3K27ac_q30.sorted.dedup.bam ENCFF512SFE_E11.5_H3K27ac_q30.sorted.dedup.bam ENCFF515PKL_E11.5_H3K27ac_q30.sorted.dedup.bam ENCFF548BRR_E10.5_H3K27ac_q30.sorted.dedup.bam ENCFF584JFB_E15.5_H3K27ac_q30.sorted.dedup.bam ENCFF707WKL_E15.5_H3K27ac_q30.sorted.dedup.bam ENCFF727QTS_E15.5_H3K27ac_q30.sorted.dedup.bam ENCFF784ORI_E14.5_H3K27ac_q30.sorted.dedup.bam ENCFF825AVI_E10.5_H3K27ac_q30.sorted.dedup.bam ENCFF902HAR_E14.5_H3K27ac_q30.sorted.dedup.bam --plotFile H3K27ac_multiBAM_fingerprint_mouse.pdf --outQualityMetrics H3K27ac_multiBAM_fingerprint_metrics_mouse.txt --outRawCounts H3K27ac_multiBAM_fingerprint_rawcounts_mouse.txt --minMappingQuality 30 --skipZeros --centerReads

#plotFingerprint -p 32 -b ENCFF002HZV_E14.5_H3K4me3_q30.sorted.dedup.bam ENCFF045IPK_E10.5_H3K4me3_q30.sorted.dedup.bam ENCFF058AUT_E12.5_H3K4me3_q30.sorted.dedup.bam ENCFF117QRC_E13.5_H3K4me3_q30.sorted.dedup.bam ENCFF124TAB_E13.5_H3K4me3_q30.sorted.dedup.bam ENCFF124UYX_E10.5_H3K4me3_q30.sorted.dedup.bam ENCFF157KEH_E10.5_H3K4me3_q30.sorted.dedup.bam ENCFF182XFG_E15.5_H3K4me3_q30.sorted.dedup.bam ENCFF182ZPF_E12.5_H3K4me3_q30.sorted.dedup.bam ENCFF184CUE_E11.5_H3K4me3_q30.sorted.dedup.bam ENCFF203JQV_E12.5_H3K4me3_q30.sorted.dedup.bam ENCFF248PGK_E13.5_H3K4me3_q30.sorted.dedup.bam ENCFF258KCR_E15.5_H3K4me3_q30.sorted.dedup.bam ENCFF376FGM_E11.5_H3K4me3_q30.sorted.dedup.bam ENCFF401BKM_E15.5_H3K4me3_q30.sorted.dedup.bam ENCFF485UDC_E13.5_H3K4me3_q30.sorted.dedup.bam ENCFF665QBJ_E14.5_H3K4me3_q30.sorted.dedup.bam ENCFF717QDV_E11.5_H3K4me3_q30.sorted.dedup.bam ENCFF724DMU_E14.5_H3K4me3_q30.sorted.dedup.bam ENCFF727QTS_E15.5_H3K4me3_q30.sorted.dedup.bam ENCFF760QYZ_E11.5_H3K4me3_q30.sorted.dedup.bam ENCFF784ORI_E14.5_H3K4me3_q30.sorted.dedup.bam ENCFF825AVI_E10.5_H3K4me3_q30.sorted.dedup.bam ENCFF941QJZ_E12.5_H3K4me3_q30.sorted.dedup.bam --plotFile H3K4me3_multiBAM_fingerprint_mouse.pdf --outQualityMetrics H3K4me3_multiBAM_fingerprint_metrics_mouse.txt --outRawCounts H3K4me3_multiBAM_fingerprint_rawcounts_mouse.txt --minMappingQuality 30 --skipZeros --centerReads
