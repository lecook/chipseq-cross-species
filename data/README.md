# Data
Folder contains both raw reads for the dunnart ChIP-seq project and genomic data for mapping and performing whole genome alignments.

1. [Raw reads](#raw-reads)
2. [Genomic data](#genomic-data)

## Raw reads
Raw data files for dunnart craniofacial ChIP-seq experiment.

Experiment performed according to manufacturers instructions for the __MAGnify™ Chromatin Immunoprecipitation System (492024)__.

H3K4me3 antibody (ab8580) - https://www.abcam.com/histone-h3-tri-methyl-k4-antibody-chip-grade-ab8580.html

H3K27ac antibody (ab4729) - https://www.abcam.com/histone-h3-acetyl-k27-antibody-chip-grade-ab4729.html

DNA fragments isolated after pull down were sequenced by GENEWIZ (China). Paired-end sequencing at a depth of 40 million reads.

</br>

sample ID | replicate pool  | pull down        | concentration (ng/uL)| amount (ug) |
----------|-----------------|------------------|----------------------|-------------|
  A1      | A               | input control    | 0.39                 | 58.5        |
  A2      | A               | H3K4me3          | 0.258                | 38.7        |
  A3      | A               | H3K27ac          | 1.61                 | 241.5       |
  B1      | B               | input control    | 0.295                | 44.25       |
  B2      | B               | H3K4me3          | 0.242                | 36.3        |
  B3      | B               | H3K27ac          | 2.33                 | 349.5       |

</br>

Library construction and sequencing performed by GENEWIZ.
Raw data (Pass Filter Data) was processed by adapter trimming and low quality read removal using NGS quality control software Cutadapt (v1.9.1) to generate clean data for subsequent analysis. software:Cutadapt(version 1.9.1) The process includes the following steps:/
1. remove the adapter sequences;
2. remove the 5 'or 3' end bases of quality scores below 20;
3. remove the reads in which ‘N’ is above 10%;
4. remove reads that are less than 75 bp long after trimming.


The statistics of raw data is summarized in the table below.

Table. Raw data statistics

| Sample | length | Reads     | Bases       | Q20(%) | Q30(%) | GC(%) | N(ppm) |
|--------|--------|-----------|-------------|--------|--------|-------|--------|
| A-1    | 150.00 | 130587134 | 19588070100 | 97.95  | 95.35  | 36.18 | 1.79   |
| A-2    | 150.00 | 103324742 | 15498711300 | 96.83  | 93.28  | 51.24 | 1.65   |
| A-3    | 150.00 | 131071676 | 19660751400 | 97.04  | 93.71  | 50.57 | 1.67   |
| B-1    | 150.00 | 111574640 | 16736196000 | 97.27  | 93.12  | 36.35 | 0.94   |
| B-2    | 150.00 | 114146802 | 17122020300 | 95.59  | 89.93  | 52.60 | 0.90   |
| B-3    | 150.00 | 104714846 | 15707226900 | 95.89  | 90.50  | 50.02 | 0.88   |


The statistics of processed data is summarized in the table below.

Table. Filtered data statistics

| Sample | length | Reads     | Bases       | Q20(%) | Q30(%) | GC(%) | N(ppm) |
|--------|--------|-----------|-------------|--------|--------|-------|--------|
| A-1    | 147.44 | 130494346 | 19239781217 | 98.30  | 95.83  | 36.06 | 1.55   |
| A-2    | 147.09 | 103188262 | 15177774139 | 97.58  | 94.25  | 51.19 | 1.41   |
| A-3    | 147.00 | 130913300 | 19244653341 | 97.76  | 94.65  | 50.54 | 1.44   |
| B-1    | 146.92 | 111529316 | 16386102716 | 97.63  | 93.63  | 36.16 | 0.82   |
| B-2    | 146.50 | 114021348 | 16704280185 | 96.49  | 91.10  | 52.51 | 0.77   |
| B-3    | 146.92 | 104607414 | 15368815013 | 96.72  | 91.59  | 49.94 | 0.77   |


Column explain:\
(1) Sample: Sample name\
(2) length: Average length of the reads\
(3) Reads: Read count\
(4) Bases: Base count\
(5) Q20, Q30: The percentage of bases with quality scores (Qphred) higher than 20 or 30\
(6) GC%: The percentage of G+C in the reads\
(7) N(ppm): The number of base ‘N’ per million bases.\

Raw files not in github but can be found in NCBI after publication through GEO Series accession number GSE188990: 
- A-1_input_R1.fastq.gz
- A-1_input_R2.fastq.gz
- A-2_H3K4me3_R1.fastq.gz
- A-2_H3K4me3_R2.fastq.gz
- A-3_H3K27ac_R1.fastq.gz
- A-3_H3K27ac_R2.fastq.gz
- B-1_input_R1.fastq.gz
- B-1_input_R2.fastq.gz
- B-2_H3K4me3_R1.fastq.gz
- B-2_H3K4me3_R2.fastq.gz
- B-3_H3K27ac_R1.fastq.gz
- B-3_H3K27ac_R2.fastq.gz

## Genomic data
