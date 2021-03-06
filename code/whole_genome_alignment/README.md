# Cross-species ChIPseq analysis

## Notes on Whole Genome Alignment

__Email from Camille Berthelot:__

We have run into a similar problem with non-model species in a recent project.
Here is the pipeline that we have had good success with, which is very similar to the conservation pipeline used by UCSC:

- build a custom whole-genome alignment between the dunnart and mouse genomes with Lastz
- produce the chains and nets for this alignment (e.g. hierarchical organisation of aligned regions to get a 1-to-1 alignment over the whole genomes)
- use LiftOver to map your dunnart regions of interest to mouse (and reciprocally)

We have benchmarked this pipeline against the one we used in the 20 mammals Cell paper, and there was almost no difference in the recovered regions.

There are a number of online tutorials to build the alignment and nets using Lastz:
http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto
https://www.bioconductor.org/packages/release/bioc/vignettes/CNEr/inst/doc/PairwiseWholeGenomeAlignment.html
https://github.com/hillerlab/GenomeAlignmentTools (more fine-tuned and advanced)

__Methods from Hecker & Hiller 2020 Whole Genome Alignment:__
Going to try this method.

To compute pairwise and multiple genome alignments, we used the human hg38 assembly as the reference (Supplementary Fig. 1 shows the entire workflow). We first built pairwise alignments between human and a query species using lastz and axtChain to compute co-linear alignment chains [82]. To align placental mammals, we used previously determined lastz parameters (K = 2400, L = 3000, Y = 9400, H = 2000, and the lastz default scoring matrix) that have a sufficient sensitivity to capture orthologous exons [16]. To align chimpanzee, bonobo, and gorilla, we changed the lastz parameters (K = 4500 and L = 4500).

After building chains, we applied RepeatFiller (RRID:SCR_017414), a method that performs another round of local alignment, considering unaligning regions ≤20 kb in size that are bounded by co-linear alignment blocks up- and downstream. RepeatFiller removes any repeat masking from the unaligned region and is therefore able to detect novel alignments between repetitive regions. We have previously shown that RepeatFiller detects several megabases of aligning repetitive sequences that would be missed otherwise. After RepeatFiller, we applied chainCleaner with parameters -LRfoldThreshold = 2.5 -doPairs -LRfoldThresholdPairs = 10 -maxPairDistance = 10000 -maxSuspectScore = 100000 -minBrokenChainScore = 75000 to improve alignment specificity. Pairwise alignment chains were converted into alignment nets using a modified version of chainNet that computes real scores of partial nets. Nets were filtered using NetFilterNonNested.perl with parameters -doUCSCSynFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000, which applies the UCSC “syntenic net” score thresholds (minTopScore of 300000 and minSynScore of 200000) and keeps nested nets that align to the same locus (inversions or local translocations; net type “inv” or “syn” according to netClass) if they score ≥5,000. For the Mongolian gerbil, tarsier, Malayan flying lemur, sperm whale, Przewalski's horse, Weddell seal, Malayan pangolin, Chinese pangolin, Hoffmann's two-fingered sloth, and Cape rock hyrax that have genome assemblies with a scaffold N50 ≤1,000,000 and a contig N50 ≤100,000, we just required that nets have a score ≥100,000. For marsupials and platypus, we lowered the score threshold for nets to 10,000 and kept inv or syn nets with scores ≥3,000.

### Some definitions

In chain and net lingo, the __target__ is the reference genome sequence and the __query__ is some other genome sequence. For example, if you are viewing Human-Mouse alignments in the Human genome browser, human is the target and mouse is the query.

A __gapless block__ is a base-for-base alignment between part of the target and part of the query, possibly including mismatching bases. It has the same length in bases on the target and the query. This is the output of the most primitive alignment algorithms.

A __gap__ is a link between two gapless blocks, indicating that the target or the query has sequence that should be skipped over in order to make the best-scoring alignment. In other words, the scoring penalty for skipping over one or more bases is less than the penalty for continuing to align the sequences without skipping.

A __single-sided gap__ is a gap in which sequence in either target or query must be skipped over. A plausible explanation for needing to skip over a base in the target while not skipping a base in the query is that either the target has an inserted base or the query has a deleted base. Many alignment tools produce alignments with single-sided gaps between gapless blocks.

A __double-sided gap__ skips over sequence in both target and query because the sum of penalties for mismatching bases exceeds the penalty for extending a gap across them. This is possible only when the penalty for extending a gap is less than the penalty for creating a new gap and less than the penalty for a mismatch, and when the alignment algorithm is capable of considering double-sided gaps.

A __chain__ is a sequence of non-overlapping gapless blocks, with single- or double-sided gaps between blocks. Within a chain, target and query coords are monotonically non-decreasing (i.e. always increasing or flat). Chains are constructed by the axtChain program which finds pairwise alignments with the same target and query sequence, on the same strand, that can be merged if overlapping and joined into one longer alignment with a higher score under an affine gap-scoring system (progressively decreasing penalties for longer gaps).

* double-sided gaps are a new capability (blastz can't do that) that allow extremely long chains to be constructed.
* not just orthologs, but paralogs too, can result in good chains. but that's useful!
* chains should be symmetrical -- e.g. swap human-mouse -> mouse-human chains, and you should get approx. the same chains as if you chain swapped mouse-human blastz alignments. However, Blastz's dynamic masking is asymmetrical, so in practice those results are not exactly symmetrical. Also, dynamic masking in conjunction with changed chunk sizes can cause differences in results from one run to the next.
* chained blastz alignments are not single-coverage in either target or query unless some subsequent filtering (like netting) is done.
* chain tracks can contain massive pileups when a piece of the target aligns well to many places in the query. Common causes of this include insufficient masking of repeats and high-copy-number genes (or paralogs).

A __net__ is a hierarchical collection of chains, with the highest-scoring non-overlapping chains on top, and their gaps filled in where possible by lower-scoring chains, which in turn may have gaps filled in by lower-level chains and so on.

* I think a chain's qName also helps to determine which level it lands in, i.e. it makes a difference whether a chain's qName is the same as the top-level chain's qName or not, because the levels have meanings associated with them -- see details page.
* a net is single-coverage for target but not for query, unless it has been filtered to be single-coverage on both target and query. By convention we add "rbest" to the net filename in that case.
* because it's single-coverage in the target, it's no longer symmetrical.
* the netter has two outputs, one of which we usually ignore: the target-centric net in query coordinates. The reciprocal best process uses that output: the query-referenced (but target-centric / target single-cov) net is turned back into component chains, and then those are netted to get single coverage in the query too; the two outputs of that netting are reciprocal-best in query and target coords. Reciprocal-best nets are symmetrical again.
* nets do a good job of filtering out massive pileups by collapsing them down to (usually) a single level.
* "LiftOver chains" are actually chains extracted from nets, or chains filtered by the netting process.


### Preparation

Spartan modules

```
module load foss
module load lastz
module load ucsc/21072020
module load perl

```

#### Repeat mask dunnart genome

Create conda environment will all the dependencies:

```
conda create -n wga
conda activate wga
conda config --add channels conda-forge
conda config --add channels biocore
conda config --add channels bioconda
conda install repeatmasker
```

Run RepeatModeler to de novo find repeat regions in the dunnart genome:
```
BuildDatabase -name dunnart -engine ncbi Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr.fasta

nohup RepeatModeler -database dunnart -pa 20 >& repeatmodeler.out

```

Run RepeatMasker to mask repeats in dunnart genome (makes repeats lowercase). Run as an array for scaffolds to make it quicker.
Create commands for array slurm script: `repeatMasker.sh`

```
RepeatMasker -q -xsmall Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr.fasta -default_search_engine hmmer -trf_prgm /home/lecook/.conda/envs/wga/bin/trf -hmmer_dir /home/lecook/.conda/envs/wga/bin/
```

#### Split into scaffolds

Using faSplit from the UCSC Kent Tools

```
faSplit byName Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr.fasta smiCra1/
```


#### Create .2bit and .sizes files

```
faToTwoBit ../../dunnart/genomes/Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr_RM.fasta smiCra1.2bit
```

```
twoBitInfo smiCra1.2bit stdout | sort -k2rn > smiCra1.chrom.sizes
```

#### mm10 target genome
http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit

__mm10.2bit__ - contains the complete mouse/mm10 genome sequence in the 2bit file format.  Repeats from __RepeatMasker__ and __Tandem Repeats Finder__ (with period of 12 or less) are shown in lower case; non-repeating sequence is shown in upper case.  


### lastZ

To align placental mammals, we used previously determined lastz parameters (K = 2400, L = 3000, Y = 9400, H = 2000, and the lastz default scoring matrix) that have a sufficient sensitivity to capture orthologous exons

To align placental mammals, we used the lastz alignment parameters K = 2400, L = 3000, Y = 9400, H = 2000 and the lastz default scoring matrix, correspond- ing to parameter set 2 in Table 1. To align vertebrates, we used K = 2400, L = 3000, Y = 3400, H = 2000 and the HoxD55 scoring matrix. Citation: Increased alignment sensitivity improves the usage of genome alignments for comparative gene annotation. Nucleic Acids Res. 2017;45(14):8369–77.

Create commands for running lastZ for all scaffolds: `lastz.sh`

Repeated with vertebrate alignment parameters and HoxD55 scoring matrix. Retrieve A LOT more aligned sequences. For example scaffold00002 with mammal parameters retrieved 863M of data, while with the new parameters it's 2.5GB.

Run as an array on slurm: `array_wrapper.slurm`


#### Convert maf to axt-format
http://last.cbrc.jp/doc/maf-convert.html

```
module load last/last/1066
maf-convert axt ${tr}.maf > ${tr}.axt
```

### axtChain

We use axtChain (http://www.soe.ucsc.edu/~kent; default parameters) to build co-linear alignment chains.

```
axtChain -linearGap=loose -scoreScheme=../../bin/GenomeAlignmentTools/HoxD55.q mm10_smiCra1.axt mm10.2bit smiCra1.2bit mm10_smiCra1.chain
```


### chainMergeSort
Merge short chains into longer ones, concatenate chains and sort

```
chainMergeSort *.chain > smiCra1_mm10.chain
```

__directory:__ /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/wholeGenomeAlignment/chainMerge/smiCra1_mm10.chain


### Installing Genomic Alignment Tools (Hiller group)

```
git clone https://github.com/hillerlab/GenomeAlignmentTools.git

cd GenomeAlignmentTools/kent/src

make

# Kent binaries
PATH=$PATH:/Users/lauracook/bin/GenomeAlignmentTools/kent/bin;export PATH

export KENTSRC_DIR=/Users/lauracook/bin/GenomeAlignmentTools/kent/src/

cd /Users/lauracook/bin/GenomeAlignmentTools/src
export MACHTYPE=x86_64
make 

PATH=$PATH:/Users/lauracook/bin/GenomeAlignmentTools/src

```
### RepeatFiller

https://github.com/hillerlab/GenomeAlignmentTools

RepeatFiller [5] is a tool to incorporate newly-detected repeat-overlapping alignments into pairwise alignment chains [4]. Its runtime adds little to the computationally more expensive step of generating chains in pairwise whole-genome alignments. RepeatFiller circumvents the problem that considering all repeat-overlapping alignment seeds during whole genome alignment is computationally not feasible. Therefore, RepeatFiller only aligns local genomic regions that are bounded by colinear aligning blocks, as provided in the chains, which makes it feasible to consider all seeds including those that overlap repetitive regions. RepeatFiller application to mammalian genome alignment chains can add between 22 and 84 Mb of previously-undetected alignments that mostly originate from transposable elements [5]. This helps to comprehensively align repetitive regions and improves the annotation of conserved non-coding elements.

```
python3 RepeatFiller.py -c smiCra1_mm10.chain -T2 mm10.2bit -Q2 smiCra1.2bit

python3 RepeatFiller.py -c smiCra1_mm10_patched_sorted.chain -T2 ../../data/genomes/mm10.2bit -Q2 ../../data/genomes/smiCra1.2bit

```
RepeatFiller adds 0.3G data to the alignment - 2 Feb 2021


__directory:__ /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/wholeGenomeAlignment/repeatFiller/


### patchChain

patchChain.perl performs a highly sensitive local pairwise alignment for loci flanked by aligning blocks [1,3]. Given an alignment chain [4], it considers all chains that pass the score and span filters (optional parameters), extracts all the unaligning loci and creates local alignment jobs. After executing these alignment jobs, the newly found and the original local alignments are combined and used to produce a new set of improved chains.

This procedure is recommended for comparisons between species that are separated by >0.75 substitutions per neutral site [1].

```
patchChain.perl smiCra1_mm10.chain /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.2bit /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.2bit /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.chrom.sizes /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.chrom.sizes -chainMinScore 5000 -gapMaxSizeT 500000 -gapMaxSizeQ 500000 -gapMinSizeT 30 -gapMinSizeQ 30 -numJobs 10 -jobDir jobs -jobList jobList -outputDir pslOutput -minEntropy 1.8 -windowSize 30 -minIdentity 60 -lastzParameters "--format=axt K=1500 L=2500 M=0 T=0 W=5 Q=/data/projects/punim0586/lecook/chipseq-pipeline/cross_species/bin/GenomeAlignmentTools/example/HoxD55.q"
```

```
this results in 10 alignment jobs that are located in jobs/ and listed in 'jobList'
now execute these jobs on a compute cluster or run them sequentially by doing 'chmod +x jobList; ./jobList'

# concatenate all new results
find pslOutput -name "*.psl" | xargs -i cat {} > newAlignments.psl

# concatenate all genomewide psl lastz alignments
find psl -name "*.psl" | xargs -i cat {} > genomeWide.lastz.psl

# combine the genome-wide lastz results (the combined psl file that was used to create the input chains) and the newly found psl alignments
cat psl/genomeWide.lastz.psl patchChain/newAlignments.psl > all.psl 

patchChain adds 1.1G to the genome wide alignment file

# use axtChain from the Kent source to compute alignment chains that include the new alignments


```


This create a folder with the jobs and a jobList which can be called in an array slurm script. Then you can run the jobs in parallel.

### axtChain again...
axtChain on psl alignments (patchChain new alignments plus genomeWide lastz alignments)

```
axtChain -psl -linearGap=loose -scoreScheme=../../bin/GenomeAlignmentTools/HoxD55.q all.psl ../../data/genomes/mm10.2bit ../../data/genomes/smiCra1.2bit smiCra1_mm10_patched.chain
```

4.7G Jan 31 18:42 smiCra1_mm10.chain
5.3G Feb  1 17:52 smiCra1_mm10_patched.chain

### chainCleaner

https://github.com/hillerlab/GenomeAlignmentTools

After RepeatFiller, we applied chainCleaner with parameters -LRfoldThreshold = 2.5 -doPairs -LRfoldThresholdPairs = 10 -maxPairDistance = 10000 -maxSuspectScore = 100000 -minBrokenChainScore = 75000 to improve alignment specificity.

chainCleaner improves the specificity in genome alignment chains by detecting and removing local alignments that obscure the evolutionary history of genomic rearrangements [2]. The input is a chain file, ideally after adding alignments found with highly sensitive parameters if distal species are compared. The output is a chain file that contains re-scored and score-sorted chains after removing the local alignments from the parent chains and adding them as individual chains. The resulting output file can be used to get alignment nets by running chainNet [4].


```
chainCleaner smiCra1_mm10_repFil.chain mm10.2bit smiCra1.2bit smiCra1_mm10_repFill_chainCl.chain smiCra1_mm10_repFill_chainCl.bed -tSizes=mm10.chrom.sizes -qSizes=smiCra1.chrom.sizes -LRfoldThreshold=2.5 -doPairs -LRfoldThresholdPairs=10 -maxPairDistance=10000 -maxSuspectScore=100000 -minBrokenChainScore=75000 -linearGap=loose
```

Run chain cleaner on patched chain file
```
chainCleaner smiCra1_mm10_patched_sorted_repFil.chain ../../data/genomes/mm10.2bit ../../data/genomes/smiCra1.2bit smiCra1_mm10_patched_sorted_repFil_chainCl.chain smiCra1_mm10_patched_sorted_repFil_chainCl.bed -tSizes=../../data/genomes/mm10.sizes -qSizes=../../data/genomes/smiCra1.sizes -LRfoldThreshold=2.5 -doPairs -LRfoldThresholdPairs=10 -maxPairDistance=10000 -maxSuspectScore=100000 -minBrokenChainScore=75000 -linearGap=loose
```


`-LRfoldThreshold` = threshold for removing local alignment blocks if the score of the left and right fill of brokenChain. Default is 2.5
`-doPairs` = flag: if set, do test if pairs of chain breaking alignments can be removed
`-LRfoldThresholdPairs` = threshold for removing local alignment blocks if the score of the left and right fill of broken chains (for pairs). Default = 10
`-maxPairDistance` = only consider pairs of chain breaking alignments where the distance between the end of the upstream CBA and the start of the downstream CBA is at most that many bp (default 10000)
`-maxSuspectScore` = threshold for score of suspect subChain. If higher, do not remove suspect.

`-linearGap`=loose

__directory:__ /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/wholeGenomeAlignment/chainCleaner/


For filtering the chains, we need the size of each chromosome
```
faSize /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.fa -detailed > smiCra1.sizes
faSize /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.fa -detailed > mm10.sizes
```

### chainPreNet

```
chainPreNet smiCra1_mm10_patched_sorted_chainCl.chain /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.sizes /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.sizes smiCra1_mm10_patched_sorted_chainCl.chain
```

chainPreNet without chainCleaner

```
chainPreNet smiCra1_mm10_patched_sorted_repFil_chainCl.chain /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.sizes /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.sizes smiCra1_mm10_patched_sorted_repFil_chainCl_preNet.chain
```

### chainNet
Given a set of alignment chains, chainNet produces alignment nets, which is a hierarchical collection of chains or parts of chains that attempt to capture only orthologous alignments [4]. The original chainNet implementation approximates the score of "sub-nets" (nets that come from a part of a chain and fill a gap in a higher-level net) by the fraction of aligning bases. This can lead to a bias in case the aligning blocks of a chain are not equally distributed. We implemented a new parameter "-rescore" in chainNet that computes the real score of each subnet [2].

Make the alignment nets:

```
chainNet smiCra1_mm10_patched_sorted_repFil_chainCl_preNet.chain -minSpace=1 /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.sizes /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.sizes stdout /dev/null | netSyntenic stdin smiCra1_mm10_patched_sorted_repFil_chainCl_preNet.net
```

### netChainSubset
Creates a single chain file using only the chains that also appear in the net.

```
netChainSubset smiCra1_mm10_patched_sorted_repFil_chainCl_preNet.net ../chainCleaner/smiCra1_mm10_patched_sorted_repFil_chainCl.chain smiCra1_mm10_patched_sorted_repFil_chainCl_preNet_chainNet.chain

```
## Join chain fragments
```
chainStitchId smiCra1_mm10_patched_sorted_chainCl_netted.chain smiCra1_mm10_patched_sorted_chainCl_netted_stitched.chain

chainStitchId smiCra1_mm10_patched_sorted_repFil_chainCl_preNet_chainNet.chain smiCra1_mm10_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain

```


From the synteny files (positions), get the sequences and re-create alignments.

### Create MAF files

smiCra1_mm10_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain

```
netToAxt chainNet/smiCra1_mm10_patched_sorted_repFil_chainCl_preNet.net liftOver_chains/smiCra1_mm10_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.2bit /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.2bit stdout | axtSort stdin mm10-smiCra1.axt

axtToMaf netToAxt/mm10-smiCra1.axt /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/mm10.sizes /data/projects/punim0586/lecook/chipseq-pipeline/cross_species/data/genomes/smiCra1.sizes mm10-smiCra1.maf -tPrefix=mm10. -qPrefix=smiCra1.
```

### chainSwap

```
chainSwap mm10_to_smiCra1_patched_sorted_chainCl_netted_stitched.chain stdout | nice chainSort stdin stdout | nice gzip -c > smiCra1_to_mm10_patched_sorted_chainCl_netted_stitched.chain

chainSwap smiCra1_mm10_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain stdout | nice chainSort stdin stdout > mm10_smiCra1_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain
```

### Quality control 

MafFilter.
Using mm10-smiCra1 maf file generated using axtToMaf look at the block lengths and genomic coverage.

```
maffilter param=options-file
```

Parameter options file
```
DATA=mm10-smiCra1 //A user-defined variable, representing the input maf file, without extension
input.file=$(DATA).maf.gz  //Input maf file, gzipped
input.file.compression=gzip
output.log=$(DATA).maffilter.log //Output log file
maf.filter=                                 \
    SequenceStatistics(                     \
        statistics=(                        \
                BlockLength,                \
            AlnScore,                       \
            BlockCounts,                    \
            PairwiseDivergence(             \
                species1=mm10,              \
                species2=smiCra1)),         \
        ref_species=mm10,                   \
        file=divergence.statistics.csv) 
```

Summarise block lengths, block counts, pairwise divergence and calculate exon coverage in the maf file between species in R = `alignment_QC.R`


- Fraction of mouse exons are in the alignment, and what fraction you actually recover by using liftover.

mm10 exons from refseq (downloaded from UCSC)

dunnart exons from lifted devil annotation gtf file converted to a bed file using BEDOPS `gtf2bed` function
```
gff2bed < Scras_dunnart_assem1.0_pb-ont-illsr_flyeassem_red-rd-scfitr2_pil2xwgs2_60chr2.gff > dunnart_gff.bed
less dunnart_gff.bed | grep "exon" > dunnart_exons.bed
```

In R: `alignment_QC.R'

## LiftOver

smiCra1exonsMerged.bed
```
liftOver smiCra1exonsMerged.bed mm10_smiCra1_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain smiCra1exonsMapped.txt smiCra1exonsUnmapped.txt
```
number of dunnart exons = 193182
number of dunnart exons that can be mapped to mouse = 111693

```
liftOver mm10exonsMerged.bed smiCra1_mm10_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain mm10exonsMapped.txt mm10exonsUnmapped.txt
```
number of mouse exons = 193776
number of mouse exons that can be lifted over to dunnart = 160213

### References

[1] Sharma V, Hiller M. Increased alignment sensitivity improves the usage of genome alignments for comparative gene annotation. Nucleic Acids Res., 45(14), 8369–8377, 2017

[2] Suarez H, Langer BE, Ladde P, Hiller M. chainCleaner improves genome alignment specificity and sensitivity. Bioinformatics, 33(11):1596-1603, 2017

[3] Hiller M, Agarwal S, Notwell JH, Parikh R, Guturu H, Wenger AM, Bejerano G. Computational methods to detect conserved non-genic elements in phylogenetically isolated genomes: application to zebrafish. Nucleic Acids Res, 41(15):e151.

[4] Kent WJ, Baertsch R, Hinrichs A, Miller W, Haussler D. Evolution's cauldron: duplication, deletion, and rearrangement in the mouse and human genomes. PNAS, 100(20):11484-9, 2003

[5] Osipova E, Hecker N, Hiller M. RepeatFiller newly identifies megabases of aligning repetitive sequences and improves annotations of conserved non-exonic elements, submitted

[6] MafFilter

[7] Kent tools UCSC