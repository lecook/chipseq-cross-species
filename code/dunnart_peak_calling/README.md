# Dunnart ChIP-seq analysis
ChIP-seq for H3K4me3 and H3K27ac for craniofacial tissue collection from day of birth dunnart (_Sminthopsis crassicaudata_) pouch young.

1. [Run snakemake pipeline](#run-snakemake-pipeline)
2. [Read alignment](#read-alignment)
3. [Filtering aligned reads](#filtering-aligned-reads)
4. [QC with deepTools](#qc-with-deeptools)
5. [Cross-correlation analysis](cross-correlation-analysis)
6. [Peak calling with MACS2](peak-calling-with-MACS2)

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
snakemake -j 6 --cluster-config envs/cluster.json --cluster "sbatch -A {cluster.account} -t {cluster.time} -p {cluster.partition} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --mem {cluster.mem}" &
```

Cluster configuration file: configs/cluster.json\
Sample configuration file: configs/config.yaml\
Sample text file: configs/SSR.text\
multiQC configuration file: configs/.multiqc_config.yaml\

## Read alignment

<details><summary>_Notes on Bowtie2 Mapping Scores_</summary>
<p>

__Scores: higher = more similar__/
An alignment score quantifies how similar the read sequence is to the reference sequence aligned to. The higher the score, the more similar they are. A score is calculated by subtracting penalties for each difference (mismatch, gap, etc.) and, in local alignment mode, adding bonuses for each match./

The scores can be configured with the --ma (match bonus), --mp (mismatch penalty), --np (penalty for having an N in either the read or the reference), --rdg (affine read gap penalty) and --rfg (affine reference gap penalty) options./

__End-to-end alignment score example__/
A mismatched base at a high-quality position in the read receives a penalty of -6 by default. A length-2 read gap receives a penalty of -11 by default (-5 for the gap open, -3 for the first extension, -3 for the second extension). Thus, in end-to-end alignment mode, if the read is 50 bp long and it matches the reference exactly except for one mismatch at a high-quality position and one length-2 read gap, then the overall score is -(6 + 11) = -17./

The best possible alignment score in end-to-end mode is 0, which happens when there are no differences between the read and the reference./

__Local alignment score example__/
A mismatched base at a high-quality position in the read receives a penalty of -6 by default. A length-2 read gap receives a penalty of -11 by default (-5 for the gap open, -3 for the first extension, -3 for the second extension). A base that matches receives a bonus of +2 be default. Thus, in local alignment mode, if the read is 50 bp long and it matches the reference exactly except for one mismatch at a high-quality position and one length-2 read gap, then the overall score equals the total bonus, 2 * 49, minus the total penalty, 6 + 11, = 81./

The best possible score in local mode equals the match bonus times the length of the read. This happens when there are no differences between the read and the reference./

__Valid alignments meet or exceed the minimum score threshold__/
For an alignment to be considered “valid” (i.e. “good enough”) by Bowtie 2, it must have an alignment score no less than the minimum score threshold. The threshold is configurable and is expressed as a function of the read length. In end-to-end alignment mode, the default minimum score threshold is -0.6 + -0.6 * L, where L is the read length. In local alignment mode, the default minimum score threshold is 20 + 8.0 * ln(L), where L is the read length. This can be configured with the --score-min option. For details on how to set options like --score-min that correspond to functions, see the section on setting function options./

__Mapping quality: higher = more unique__/
The aligner cannot always assign a read to its point of origin with high confidence. For instance, a read that originated inside a repeat element might align equally well to many occurrences of the element throughout the genome, leaving the aligner with no basis for preferring one over the others./

Aligners characterize their degree of confidence in the point of origin by reporting a mapping quality: a non-negative integer Q = -10 log10 p, where p is an estimate of the probability that the alignment does not correspond to the read’s true point of origin. Mapping quality is sometimes abbreviated MAPQ, and is recorded in the SAM MAPQ field./

Mapping quality is related to “uniqueness.” We say an alignment is unique if it has a much higher alignment score than all the other possible alignments. The bigger the gap between the best alignment’s score and the second-best alignment’s score, the more unique the best alignment, and the higher its mapping quality should be./

Accurate mapping qualities are useful for downstream tools like variant callers. For instance, a variant caller might choose to ignore evidence from alignments with mapping quality less than, say, 10. A mapping quality of 10 or less indicates that there is at least a 1 in 10 chance that the read truly originated elsewhere./
</p>
</details>

__Parameters__

- `-X <int>`: Maximum DNA fragment length (default 500bp). If you anticipate that you may have DNA fragments longer than the default value, you should increase this parameter accordingly; other, alignments from such fragment are considered no properly paired). See figure below.
- `--very sensitive`: Bowtie2 has a number of alignment and effort parameters that interact in a complex (and sometimes unexpected) ways. Present collections of these parameters are provided for convenience; the defaults is `--sensitive`, but better alignment results are frequently achieved with this option. Very sensitive is the same as setting -D 20 -R 3 -N 0 -L 20 -i S,1,0.50.
- `-D`: Up to <int> consecutive seed extension attempts can “fail” before Bowtie 2 moves on, using the alignments found so far. A seed extension “fails” if it does not yield a new best or a new second-best alignment. Default: 15.
- `-R`: <int> is the maximum number of times Bowtie 2 will “re-seed” reads with repetitive seeds. When “re-seeding,” Bowtie 2 simply chooses a new set of reads (same length, same number of mismatches allowed) at different offsets and searches for more alignments. A read is considered to have repetitive seeds if the total number of seed hits divided by the number of seeds that aligned at least once is greater than 300. Default: 2.
- `-N`: Sets the number of mismatches to allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) but increases sensitivity. Default: 0.
- `-L`: Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default 22.
- `-i`: Sets a function governing the interval between seed substrings to use during multiseed alignment. Since it’s best to use longer intervals for longer reads, this parameter sets the interval as a function of the read length, rather than a single one-size-fits-all number. __We used S,1,0.5 sets interval function f to f(x) = 1 + 0.5 * sqrt(x)__ which for 150bp read length ours seed length is __18bp__.
- `-x`: genome index name
- `-1`: pair 1
- `-2`: pair 2


__Effective Genome Length__

We can approximate effective genome size for various read lengths using the khmer program and `unique-kmers.py`. This will estimate the number of unique kmers (for a specified length kmer) which can be used to infer the total uniquely mappable genome. (I.e it doesn't include highly repetitive regions). https://khmer.readthedocs.io/en/v2.1.1/user/scripts.html

This was a suggestion of deepTools: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

Install khmer program:
```{bash eval=FALSE}
module load 
module load 
pip2.7 install khmer
```

Run `unique-kmers.py` on dunnart genome for read length of 150bp:

```{bash eval=FALSE}
/usr/local/bin/unique-kmers.py -k 150 Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr.fa

/usr/local/bin/unique-kmers.py -k 150 Sminthopsis_crassicaudata_HiC.fasta

```
```
Estimated number of unique 150-mers in /Users/lauracook/../../Volumes/macOS/genomes/Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr.fasta: 2740338543
Total estimated number of unique 150-mers: 2740338543
```


__Indexing genome file__

Build Index

__Load modules:__

```{bash eval=FALSE}
module load gcc/8.3.0
module load bowtie2/2.3.5.1
```

__Build index__

```{bash eval=FALSE}
bowtie2-build Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr.fa Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr

bowtie2-build Sminthopsis_crassicaudata_HiC.fasta Sminthopsis_crassicaudata_HiC
```

## Filtering aligned reads

### rule filter:

Uses samtools and picards MarkDuplicates to remove low quality reads. This is based on and in accordance with the ENCODE Guidelines and pipeline (https://github.com/ENCODE-DCC/chip-seq-pipeline2)

- `-q 30`: remove reads with MAPQ score below 30
- `-f 2`: keep proper pairs

### rule markDups:

`picard MarkDuplicates`: marks optical and PCR duplicate reads (eg. reads that have exactly the same start coordinate)

### rule dedup:

`-F 1804`: removes unmapped, mate unmapped, not primary alignemnt, read fails platform/vendor quality checks, read is PCR or optical duplicate (marked by picard).

### rule indexBam:

Use `samtools index` to index BAM file for use with deeptools.

### rule mappingStats:

Use `samtools flagstat` to get stats at each step in the filtering so I know exactly how many reads we lose at each filtering step.

### rule downsample_bam:

`Picard DownsampleSam` is used to sample initial experiment for calculating library complexity. The `-P` (probability) parameter is adjusted in such a way that if it is paired end sequencing, around 10M mapped mates are sampled.

### rule preseq:

The preseq package is aimed at predicting and estimating the complexity of a genomic sequencing library, equivalent to predicting and estimating the number of redundant reads from a given sequencing depth and how many will be expected from additional sequencing using an initial sequencing experiment. The estimates can then be used to examine the utility of further sequencing, optimize the sequencing depth, or to screen multiple libraries to avoid low complexity samples.

### rule estimate_lib_complexity:

This follows exactly the code used the in thee ENCODE chip-seq pipeline for computing the non-redundant fraction (NRF), PCR bottlenecking cofficient 1 & 2 (PBC1 and PBC2). The non-redundant fraction is the proportion of reads that have a distinct start coordinate in read pair 1 and a distinct end coordinate in read pair 2 (as opposed to duplicate reads that have exactly the same coordinate) out of all uniquely mapped read pairs.
Library complexity is calculated from a subsample of 10-15M reads from each BAM file.

1. Sort by name
2. Convert to bedPE and obtain fragment coordinates
3. Sort by position and strand
4. Obtain unique count statistics

The output file looks like this:

TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

Library Complexity ChIP-seq Standards:

| PBC1 | PBC2 | Bottlenecking level | NRF | Complexity | Flag colors |
|:-:|:-:|:-:|:-:|:-:|:-:|
| < 0.5 | < 1 | Severe | < 0.5 | Concerning | Orange |
| 0.5 ≤ PBC1 < 0.8 | 1 ≤ PBC2 < 3 | Moderate | 0.5 ≤ NRF < 0.8 | Acceptable | Yellow |
| 0.8 ≤ PBC1 < 0.9 | 3 ≤ PBC2 < 10 | Mild | 0.8 ≤ NRF < 0.9 | Compliant | None |
| ≥ 0.9 | ≥ 10 | None | > 0.9 | Ideal | None |

##  QC with deepTools

### rule deeptools_summary:

computes the read coverages for genomic regions for typically two or more BAM files. The standard output of `multiBamSummary` is a compressed numpy array (`.npz`). It can be directly used to calculate and visualize pairwise correlation values between the read coverages using the tool `plotCorrelation`.

### rule deeptools_correlation:

Pearson or Spearman methods are available to compute correlation coefficients. Results can be saved as multiple scatter plots depicting the pairwise correlations or as a clustered heatmap, where the colors represent the correlation coefficients and the clusters are constructed using complete linkage. Optionally, the values can be saved as tables, too.

`plotCorrelation` computes the overall similarity between two or more files based on read coverage (or other scores) within genomic regions.
This helps to determine whether the different sample types can be separated, i.e., samples of different conditions are expected to be more dissimilar to each other than replicates within the same condition.

### rule deeptools_coverage:

Normalised to the reads per genomic content (normalized to 1x coverage). Produces a coverage file.

The bigWig format is an indexed binary format useful for dense, continuous data that will be displayed in a genome browser as a graph/track, but also is used as input for some of the visualization commands in deepTools.

- `normalizeUsing`: Possible choices: RPKM, CPM, BPM, RPGC. We will use BPM (Bins Per Million), which is similar to TPM in RNA-seq. BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions).
- `binSize`: size of bins in bases
- `smoothLength`: defines a window, larger than the binSize, to average the number of reads over. This helps produce a more continuous plot.
- `centerReads`: reads are centered with respect to the fragment length as specified by extendReads. This option is useful to get a sharper signal around enriched regions.

### rule deeptools_fingerprint:

This tool samples indexed BAM files and plots a profile of cumulative read coverages for each. All reads overlapping a window (bin) of the specified length are counted; these counts are sorted and the cumulative sum is finally plotted.

It determines how well the signal in the ChIP-seq sample can be differentiated from the background distribution of reads in the control sample.

An ideal input with perfect uniform distribution of reads along the genome (i.e. without enrichments in open chromatin etc.) and infinite sequencing coverage should generate a straight diagonal line. A very specific and strong ChIP enrichment will be indicated by a prominent and steep rise of the cumulative sum towards the highest rank. This means that a big chunk of reads from the ChIP sample is located in few bins which corresponds to high, narrow enrichments typically seen for transcription factors.

### rule deeptools_bamPEFragmentSize:

This tool calculates the fragment sizes for read pairs given a BAM file from paired-end sequencing. Several regions are sampled depending on the size of the genome and number of processors to estimate thesummary statistics on the fragment lengths.

##  Cross-correlation analysis

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

### rule trim_read1:

Trim R1 fastq to 50bp using the ENCODE script `scripts/trimfastq.py`.

### rule align_trimmed_read1:

Align trimmed read 1 (not paired) with bowtie2 and use it for filtering step (no deduped).

### rule filter_sort_trimmed_alignment:

Uses samtools and picards MarkDuplicates to remove low quality reads. This is based on and in accordance with the ENCODE Guidelines and pipeline (https://github.com/ENCODE-DCC/chip-seq-pipeline2)

- `-q 30`: remove reads with MAPQ score below 30
- `-f 2`: keep proper pairs

### rule bamtobed_crossC:

Make a tagAlign file for filtered BAM and subsample it for cross-correlation analysis. Use `bedtools bamtobed`

### rule subsample_aligned_reads:

Subsample tagAlign file for 15M reads.

### rule cross_correlation_SSP:

Need to determine the exclusion range for fragment length estimation. (use the fragment length calculated by deepTools `bamPEFragmentSize` above).

Notes from ENCODE pipeline:
1. Use a fixed lowerbound at -500.
2. Upperbound EXCLUSION_RANGE_MAX is
    - TF ChIP-seq:  max(read_len + 10, 50)
    - Histone ChIP-seq:  max(read_len + 10, 100)

__CC_SCORE FILE format__
Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag

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

## Calling peaks after normalising to 10M reads

This is so we can compare to mouse peaks.
Based on https://davemcg.github.io/post/easy-bam-downsampling/ script for doing this. 

```
bash subsample.sh
```

Run script in the directory with the files you want to subsample.

### rule call_peaks_macs2:

- `peaks.narrowPeak`: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
- `peaks.xls`: a tabular file which contains information about called peaks. Additional information includes pileup and fold enrichment
- `summits.bed`: peak summits locations for every peak. To find the motifs at the binding sites, this file is recommended


#  Peak QC


### rule get_narrow_peak_counts_for_multiqc:
Peak counts into a format usable for multiQC. Use Adrian Foucal script `count_peaks.py` (https://gitlab.univ-nantes.fr/foucal-a/full-chipseq/-/blob/master/scripts/count_peaks.py)

### rule bamToBed:

Convert BAM to tagAlign file for calculating FRiP QC metric (Fraction of reads in peaks) using `bedtools bamtobed`.

`-mate1`:	When writing BEDPE (-bedpe) format, always report mate one as the first BEDPE “block”.


### rule frip:

Calculate fraction of reads in peaks using a custom script based on ENCODE pipeline code - `scripts/encode_frip.py`


# Create consensus peaksets for replicates

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

### rule pool_peaks:

Overlap peaks scripts takes a file of pooled peaks (repA merged with repB peaks).

### rule overlap_peaks_H3K4me3:

Overlap replicate peaks for H3K4me3.

### rule overlap_peaks_H3K27ac:

Overlap replicate peaks for H3K27ac.

# Consensus peaks and grouping into putative enhancers and promoters
## Merge peaks with max distance between peaks of 50bp 

```
sort -k1,1 -k2,2n H3K4me3_overlap.narrowPeak > H3K4me3_overlap.sorted.narrowPeak
sort -k1,1 -k2,2n H3K27ac_overlap.narrowPeak > H3K27ac_overlap.sorted.narrowPeak
bedtools merge -i H3K4me3_overlap.sorted.narrowPeak -d 50 > H3K4me3_overlap_merge.narrowPeak
bedtools merge -i H3K27ac_overlap.sorted.narrowPeak -d 50 > H3K27ac_overlap_merge.narrowPeak
```

Find uniquely H3K4me3 sites (i.e. peaks that don't overlap with H3K27ac):

```
bedtools intersect -v -a H3K4me3_overlap.narrowPeak -b H3K27ac_overlap.narrowPeak > H3K4me3_only.narrowPeak
```

Find unique H3K27ac sites (i.e. peaks that don't overlap with H3K4me3):
```
bedtools intersect -v -a H3K27ac_overlap.narrowPeak -b H3K4me3_overlap.narrowPeak > H3K27ac_only.narrowPeak
```

Find peaks common between H3K27ac & H3K4me3 with a reciprocal overlap of at least 50%.
```
bedtools intersect -f 0.5 -r -a H3K4me3_overlap.narrowPeak -b H3K27ac_overlap.narrowPeak > H3K4me3_and_H3K27ac.narrowPeak
```

Combine files

```
cat H3K4me3_only.narrowPeak H3K4me3_and_H3K27ac.narrowPeak | sort -k 4,4 > dunnart_promoter_peaks.narrowPeak
cat H3K27ac_only.narrowPeak | sort -k 4,4 > dunnart_enhancer_peaks.narrowPeak

sort -k1,1 -k2,2n dunnart_enhancers.narrowPeak | bedtools merge -i stdin -d 50 > dunnart_enhancer_mergedPeaks.narrowPeak
```

### Peak Features, Annotation, GO, motif enrichment
See analysis directory.
