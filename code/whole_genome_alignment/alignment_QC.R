## The scripts looks at some basic QC to assess the 
## mm10-smiCra1 liftOver chains

library(plyr)
library(dplyr)
library(data.table)
library(seqinr)
library(Biostrings)

setwd("/Users/lauracook/OneDrive - The University of Melbourne/PhD (2018-2021)/8-WholeGenomeAlignment/2-UCSC_hoxD55_alignment/1-UCSC_alignment_QC/2-MafFilter/")
# Summarise statistics from MafFilter
mm10 <- read.table("mm10.statistics.csv", header=TRUE, sep="\t") # read in table
smiCra1 <- read.table("smiCra1.statistics.csv", header=TRUE, sep="\t") # read in table
head(mm10) #check file has been imported correctly

mm10ChrCount <- mm10 %>% count(Chr, sort=TRUE)#how many blocks per chromosome
smiCra1ChrCount <- smiCra1 %>% count(Chr, sort=TRUE) #how many blocks per chromosome

#total size (incl gaps) just use mouse as the statistics are the same
size1to10bp <- mm10[which(mm10$BlockLength>=1 & mm10$BlockLength<=10),]
size10to100bp <- mm10[which(mm10$BlockLength>=10 & mm10$BlockLength<=100),]
size100to1kb <- mm10[which(mm10$BlockLength>=100 & mm10$BlockLength<=1000),]
size1kbto10kb <- mm10[which(mm10$BlockLength>=1000 & mm10$BlockLength<=10000),]
size10kbto100kb <- mm10[which(mm10$BlockLength>=10000 & mm10$BlockLength<=100000),]
size100kbto1Mb <- mm10[which(mm10$BlockLength>=100000 & mm10$BlockLength<=1000000),]

count(size1to10bp)
count(size10to100bp)
count(size100to1kb)
count(size1kbto10kb)
count(size10kbto100kb)
count(size100kbto1Mb)

sum(size1to10bp$BlockLength)
sum(size10to100bp$BlockLength)
sum(size100to1kb$BlockLength)
sum(size1kbto10kb$BlockLength)
sum(size10kbto100kb$BlockLength)
sum(size100kbto1Mb$BlockLength)

setwd("/Users/lauracook/OneDrive - The University of Melbourne/PhD (2018-2021)/8-WholeGenomeAlignment/2-UCSC_hoxD55_alignment/1-UCSC_alignment_QC/")

# MOUSE exon coverage
#mm10 refseq ncbi curated exons downloaded from ucsc
system("bedtools sort -i mm10exons.bed > mm10exonsSorted.bed", intern=TRUE) ## sort by coordinates
system("bedtools merge -i mm10exonsSorted.bed > mm10exonsMerged.bed", intern=TRUE) ## bedtools merge where transcripts are the same coordinates or within each other
system("bedtools merge -i mm10_mafAlignment.bed > mm10_mafAlignmentMerged.bed", intern=TRUE) ## also merge the maf alignments from the liftOver
system("bedtools intersect -a mm10exonsMerged.bed -b mm10_mafAlignmentBlocks.bed -wo > mm10exonOverlaps.bed", intern=TRUE) ## intersect exon coordinates with maf alignment file

## calculate coverage 
mm10exon <- read.table("mm10exonsMerged.bed") ## exon file
mm10exon$length <- (mm10exon$V3 - mm10exon$V2) ## exon length in ncbi refseq 
mm10overlaps <- read.table("mm10exonOverlaps.bed") ## overlap file
dim(mm10overlaps) ## check file imported correctly
sum(mm10overlaps$V10) ## 31776094bps of exons found in maf alignment
sum(mm10exon$length) ## 34221920bps of exons in refseq mm10
mm10exonCoverage <- sum(mm10overlaps$V10)/sum(mm10exon$length) ## converge alignment exon bps/total exon bps

## DUNNART exon coverage 
system("bedtools sort -i dunnart_exons.bed > dunnart_exonsSorted.bed", intern = TRUE) ## sort by coordinates
system("bedtools merge -i dunnart_exonsSorted.bed > smiCra1exonsMerged.bed", intern = TRUE) ## bedtools merge where transcripts are the same coordinates or within each other
system("bedtools merge -i smiCra1_mafAlignmentBlocksSorted.bed > smiCra1_mafAlignmentMerged.bed", intern = TRUE) ## also merge the maf alignments from the liftOver
system("bedtools intersect -a smiCra1exonsMerged.bed -b smiCra1_mafAlignmentMerged.bed -wo > smiCra1_exonOverlaps.bed", intern=TRUE) ## intersect exon coordinates with maf alignment file
smiCra1exon <- read.table("smiCra1exonsMerged.bed", sep="\t", header=FALSE) ## exon file
dim(smiCra1exon) ## check file imported correctly
smiCra1exon$length <- (smiCra1exon$V3 - smiCra1exon$V2)  ## exon length
smiCra1overlaps <- read.table("smiCra1_exonOverlaps.bed", sep="\t", header=FALSE) ## exon overlaps with maf alignment file
dim(smiCra1overlaps) ## check imported correctly
sum(smiCra1exon$length) # 63502228
sum(smiCra1overlaps$V7) # 43939135
smiCra1exonCoverage <- sum(smiCra1overlaps$V7)/sum(smiCra1exon$length)

## average sequence percentage mismatches for dunnart and mouse
mm10 <- read.table("mm10.divergence.statistics.csv", header=TRUE, sep="\t")
smiCra1 <- read.table("smiCra1.divergence.statistics.csv", header=TRUE, sep="\t")

mean(mm10$Div.mm10.smiCra1)
mean(smiCra1$Div.smiCra1.mm10)

## Genome coverage
readDNAStringSet("mm10_mafAlignment_nogaps.fasta") #number bp in alignment 743350114
readDNAStringSet("smiCra1_mafAlignment_nogaps.fasta") #number bp in alignment 747935780

mm10genomeSize <- 2652783500
smiCra1genomeSize <-2838290115

mm10coverage <- 743350114/mm10genomeSize
smiCra1coveage <- 747935780/smiCra1genomeSize
