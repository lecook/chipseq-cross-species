# Mouse ChIP-seq analysis

ChIP-seq aligned reads downloaded from ENCODE for H3K4me3 and H3K27ac for embryonic facial prominances collected from embryos (E10.5 to E15.5). Further details on the samples can be found in the master folder README.

1. [Run snakemake pipeline](#run-snakemake-pipeline)
2. [Filtering aligned reads](#filtering-aligned-reads)
3. [QC with deepTools](#qc-with-deeptools)
4. [Cross-correlation analysis](cross-correlation-analysis)
5. [Peak calling with MACS2](peak-calling-with-MACS2)


## Run snakemake pipeline

Create conda environment:

```
# create conda environment with all packages and dependencies
conda create --name chip bowtie2 samtools picard deeptools multiqc phantompeakqualtools preseq macs2 bedtools fastqc pybedtools

# save environment to a yaml file
conda env export > chip_environment.yml

# activate environment before running pipeline
conda activate chip
```

Run pipeline with the following command:

```
snakemake --snakefile [snakemake file] -j 6 --cluster-config envs/cluster.json --cluster "sbatch -A {cluster.account} -t {cluster.time} -p {cluster.partition} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --mem {cluster.mem}" &
```

H3K27ac:
```
snakemake -j 6 --snakefile Snakefile_H3K27ac --cluster-config configs/cluster.json --cluster "sbatch -A {cluster.account} -t {cluster.time} -p {cluster.partition} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --mem {cluster.mem}" &

```

H3K4me3:
```
snakemake -j 6 --snakefile Snakefile_H3K4me3 --cluster-config configs/cluster.json --cluster "sbatch -A {cluster.account} -t {cluster.time} -p {cluster.partition} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --mem {cluster.mem}" &

```

Cluster configuration file: configs/cluster.json\
Sample configuration file: configs/config.yaml\
Sample text file: configs/SSR.text\
multiQC configuration file: configs/.multiqc_config.yaml\

## Filtering aligned reads

add details here.

## QC with deepTools
### rule deeptools_coverage:
Normalised to the reads per genomic content (normalized to 1x coverage)
Produces a coverage file

The bigWig format is an indexed binary format useful for dense, continuous data that will be displayed in a genome browser as a graph/track, but also is used as input for some of the visualization commands in deepTools.

- `normalizeUsing`: Possible choices: RPKM, CPM, BPM, RPGC. We will use BPM (Bins Per Million), which is similar to TPM in RNA-seq. BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions).
- `binSize`: size of bins in bases
- `smoothLength`: defines a window, larger than the binSize, to average the number of reads over. This helps produce a more continuous plot.
- `centerReads`: reads are centered with respect to the fragment length as specified by extendReads. This option is useful to get a sharper signal around enriched regions.

## Cross-correlation analysis

Information from: https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit

This set of programs operate on mapped Illumina single-end read datasets in tagAlign or BAM format. Because my data is paired-end I need to only use the forward read.

A high-quality ChIP-seq experiment will produce significant clustering of enriched DNA sequence tags/reads at locations bound by the protein of interest; the expectation is that we can observe a bimodal enrichment of reads (sequence tags) on both the forward and the reverse strands.

Cross-correlation analysis is done on a filtered (but not-deduped) and subsampled BAM. There is a special fastq trimming for cross-correlation analysis. Read1 fastq is trimmed to 50bp first using trimfastq.py (last modified 2017/11/08, https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/src/trimfastq.py). And then it is separately mapped as SE. Reads are filtered but duplicates are not removed. Then 15 million reads are randomly sampled and used for cross-correlation analysis.

ENCODE standards:

| col. | abbreviation     | description                                                                                         |
|------|------------------|-----------------------------------------------------------------------------------------------------|
| 1    | Filename         | tagAlign/BAM filename                                                                               |
| 2    | numReads         | effective sequencing depth i.e. total number of mapped reads in input file                          |
| 3    | estFragLen       | comma separated strand cross-correlation peak(s) in decreasing order of correlation.                |
| 4    | corr_estFragLen  | comma separated strand cross-correlation value(s) in decreasing order (COL2 follows the same order) |
| 5    | phantomPeak      | Read length/phantom peak strand shift                                                               |
| 6    | corr_phantomPeak | Correlation value at phantom peak                                                                   |
| 7    | argmin_corr      | strand shift at which cross-correlation is lowest                                                   |
| 8    | min_corr         | minimum value of cross-correlation                                                                  |
| 9    | NSC              | Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8                                 |
| 10   | RSC              | Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8)                 |
| 11   | QualityTag       | Quality tag based on thresholded RSC (codes= -2:veryLow, -1:Low, 0:Medium, 1:High, 2:veryHigh)      |


NSC; NSC>1.1 (higher values indicate more enrichment; 1 = no enrichment)

RSC; RSC>0.8 (0 = no signal; <1 low quality ChIP; >1 high enrichment

Quality tag based on thresholded RSC (codes: -2:veryLow,-1:Low,0:Medium,1:High; 2:veryHigh)


## Peak calling with MACS2

__Input file options__

- `-t`: The IP data file (this is the only REQUIRED parameter for MACS)
- `-c`: The control or mock data file
- `-f`: format of input file; Default is “AUTO” which will allow MACS to decide the format automatically.
- `-f BAMPE`: Here, the fragments are defined by the paired alignments' ends, and there is no modelling or artificial extension. Singleton alignments are ignored. This is the preferred option for using only properly paired alignments.
- `-g`: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.

__Output arguments__

- `--outdir`: MACS2 will save all output files into speficied folder for this option
- `-n`:The prefix string for output files

__Shifting model arguments__

- `-s`: size of sequencing tags. Default, MACS will use the first 10 sequences from your input treatment file to determine it
- `--bw`: The bandwidth which is used to scan the genome ONLY for model building. Can be set to the expected sonication fragment size.
- `--mfold`: upper and lower limit for model building


__Peak calling arguments__

- `-q`: q-value (minimum FDR) cutoff, default 0.05
- `--nolambda`: do not consider the local bias/lambda at peak candidate regions
- `--broad`: broad peak calling

- `--keep-dup all`: PCR duplicates have been removed by another program, such as Picard's MarkDuplicates. So need to specify to keep all duplicates. Default is to remove all potential PCR duplicates. Picard's MarkDuplicates is better.

I've left all the shifting model and peak calling arguments as default


### rule call_peaks_macs2:

- `peaks.narrowPeak`: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
- `peaks.xls`: a tabular file which contains information about called peaks. Additional information includes pileup and fold enrichment
- `summits.bed`: peak summits locations for every peak. To find the motifs at the binding sites, this file is recommended

#  Create consensus peaksets for replicates

Edited version of ENCODE `overlap_peaks.py` - recommended for histone marks.

Overlap with `BEDTools intersect` with the following parameters

- `-a` = bed file 1
- `-b` = bed file 2
- `-wo` = Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported.
- `s1` = length of pooled peak
- `s2` = length of replicate peak
- `$21` = amount that overlaps between the two bed coordinates that intersect
- `FS` = field-separator
- `OFS` = output field-separator

therefore if amount that overlaps between each replicate divided by the length of each peak needs to be greater than 50%.


# Consensus peaks and grouping into putative enhancers and promoters

Find uniquely H3K4me3 sites (i.e. peaks that don't overlap with H3K27ac):

```
bedtools intersect -v -a H3K4me3_E10.5_overlap.narrowPeak -b H3K27ac_E10.5_overlap.narrowPeak > H3K4me3_E10.5_only.narrowPeak
bedtools intersect -v -a H3K4me3_E11.5_overlap.narrowPeak -b H3K27ac_E11.5_overlap.narrowPeak > H3K4me3_E11.5_only.narrowPeak
bedtools intersect -v -a H3K4me3_E12.5_overlap.narrowPeak -b H3K27ac_E12.5_overlap.narrowPeak > H3K4me3_E12.5_only.narrowPeak
bedtools intersect -v -a H3K4me3_E13.5_overlap.narrowPeak -b H3K27ac_E13.5_overlap.narrowPeak > H3K4me3_E13.5_only.narrowPeak
bedtools intersect -v -a H3K4me3_E14.5_overlap.narrowPeak -b H3K27ac_E14.5_overlap.narrowPeak > H3K4me3_E14.5_only.narrowPeak
bedtools intersect -v -a H3K4me3_E15.5_overlap.narrowPeak -b H3K27ac_E15.5_overlap.narrowPeak > H3K4me3_E15.5_only.narrowPeak

```

Find unique H3K27ac sites (i.e. peaks that don't overlap with H3K4me3):
```
bedtools intersect -v -a H3K27ac_E10.5_overlap.narrowPeak -b H3K4me3_E10.5_overlap.narrowPeak > H3K27ac_E10.5_only.narrowPeak
bedtools intersect -v -a H3K27ac_E11.5_overlap.narrowPeak -b H3K4me3_E11.5_overlap.narrowPeak > H3K27ac_E11.5_only.narrowPeak
bedtools intersect -v -a H3K27ac_E12.5_overlap.narrowPeak -b H3K4me3_E12.5_overlap.narrowPeak > H3K27ac_E12.5_only.narrowPeak
bedtools intersect -v -a H3K27ac_E13.5_overlap.narrowPeak -b H3K4me3_E13.5_overlap.narrowPeak > H3K27ac_E13.5_only.narrowPeak
bedtools intersect -v -a H3K27ac_E14.5_overlap.narrowPeak -b H3K4me3_E14.5_overlap.narrowPeak > H3K27ac_E14.5_only.narrowPeak
bedtools intersect -v -a H3K27ac_E15.5_overlap.narrowPeak -b H3K4me3_E15.5_overlap.narrowPeak > H3K27ac_E15.5_only.narrowPeak
```

Find peaks common between H3K27ac & H3K4me3 with a reciprocal overlap of at least 50%.
```
bedtools intersect -f 0.5 -r -a H3K4me3_E10.5_overlap.narrowPeak -b H3K27ac_E10.5_overlap.narrowPeak > H3K4me3_and_H3K27ac_E10.5.narrowPeak
bedtools intersect -f 0.5 -r -a H3K4me3_E11.5_overlap.narrowPeak -b H3K27ac_E11.5_overlap.narrowPeak > H3K4me3_and_H3K27ac_E11.5.narrowPeak
bedtools intersect -f 0.5 -r -a H3K4me3_E12.5_overlap.narrowPeak -b H3K27ac_E12.5_overlap.narrowPeak > H3K4me3_and_H3K27ac_E12.5.narrowPeak
bedtools intersect -f 0.5 -r -a H3K4me3_E13.5_overlap.narrowPeak -b H3K27ac_E13.5_overlap.narrowPeak > H3K4me3_and_H3K27ac_E13.5.narrowPeak
bedtools intersect -f 0.5 -r -a H3K4me3_E14.5_overlap.narrowPeak -b H3K27ac_E14.5_overlap.narrowPeak > H3K4me3_and_H3K27ac_E14.5.narrowPeak
bedtools intersect -f 0.5 -r -a H3K4me3_E15.5_overlap.narrowPeak -b H3K27ac_E15.5_overlap.narrowPeak > H3K4me3_and_H3K27ac_E15.5.narrowPeak
```

Combine peak files

```
cat H3K4me3_E10.5_only.narrowPeak H3K4me3_and_H3K27ac_E10.5.narrowPeak | sort -k 4,4 > E10.5_promoter_peaks.narrowPeak
cat H3K4me3_E11.5_only.narrowPeak H3K4me3_and_H3K27ac_E11.5.narrowPeak | sort -k 4,4 > E11.5_promoter_peaks.narrowPeak
cat H3K4me3_E12.5_only.narrowPeak H3K4me3_and_H3K27ac_E12.5.narrowPeak | sort -k 4,4 > E12.5_promoter_peaks.narrowPeak
cat H3K4me3_E13.5_only.narrowPeak H3K4me3_and_H3K27ac_E13.5.narrowPeak | sort -k 4,4 > E13.5_promoter_peaks.narrowPeak
cat H3K4me3_E14.5_only.narrowPeak H3K4me3_and_H3K27ac_E14.5.narrowPeak | sort -k 4,4 > E14.5_promoter_peaks.narrowPeak
cat H3K4me3_E15.5_only.narrowPeak H3K4me3_and_H3K27ac_E15.5.narrowPeak | sort -k 4,4 > E15.5_promoter_peaks.narrowPeak


cat H3K27ac_E10.5_only.narrowPeak | sort -k 4,4 > E10.5_enhancer_peaks.narrowPeak
cat H3K27ac_E11.5_only.narrowPeak | sort -k 4,4 > E11.5_enhancer_peaks.narrowPeak
cat H3K27ac_E12.5_only.narrowPeak | sort -k 4,4 > E12.5_enhancer_peaks.narrowPeak
cat H3K27ac_E13.5_only.narrowPeak | sort -k 4,4 > E13.5_enhancer_peaks.narrowPeak
cat H3K27ac_E14.5_only.narrowPeak | sort -k 4,4 > E14.5_enhancer_peaks.narrowPeak
cat H3K27ac_E15.5_only.narrowPeak | sort -k 4,4 > E15.5_enhancer_peaks.narrowPeak


```



# Plot DAG

```
snakemake --dag | dot -Tsvg > dag.svg
```


# Annotate peaks
