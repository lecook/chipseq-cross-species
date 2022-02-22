# chipseq_cross_species

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

## Dunnart ChIP-seq data
<details><summary>Experimental design</summary>
<p>
  
__H3K4me3__

- Signature of active promoters.
- Closely linked with TSSs.
- Active, and prefers promoters to enhancers

__H3K27ac__

- Signature of active and poised enhancers
- Active, and prefers enhancers to promoters

__Pouch young heads were pooled to generate 2 replicates__


| collection date | PY ID   | weight (g) | head_shape | maternal | paternal | litter | sex | replicate pool |
|-----------------|---------|------------|------------|----------|----------|--------|-----|----------------|
| 17/09/2018      | Fr113.1 | 0.014      | flat head  | Fr113    | My76     | 1      |     | A              |
| 17/09/2018      | Fr113.2 | 0.016      | flat head  | Fr113    | My76     | 1      |     | A              |
| 17/09/2018      | Fr113.3 | 0.0145     | flat head  | Fr113    | My76     | 1      |     | A              |
| 17/09/2018      | Fr113.4 | 0.0162     | flat head  | Fr113    | My76     | 1      |     | A              |
| 14/11/2018      | Fb148.1 | 0.016      | round head | Fb148    | Mg120    | 2      |     | B              |
| 14/11/2018      | Fb148.2 | 0.016      | round head | Fb148    | Mg120    | 2      |     | B              |
| 14/11/2018      | Fb148.3 | 0.016      | round head | Fb148    | Mg120    | 2      |     | B              |
| 14/11/2018      | Fb148.4 | 0.016      | round head | Fb148    | Mg120    | 2      |     | B              |
| 9/01/2019       | Fg123.1 | 0.0155     | flat head  | Fg123    | Mr89     | 3      |     | B              |
| 9/01/2019       | Fg123.2 | 0.0155     | flat head  | Fg123    | Mr89     | 3      |     | B              |
| 9/01/2019       | Fg123.3 | 0.0155     | flat head  | Fg123    | Mr89     | 3      |     | B              |
| 9/01/2019       | Fg123.4 | 0.0155     | flat head  | Fg123    | Mr89     | 3      |     | B              |
| 9/01/2019       | Fg123.5 | 0.0155     | flat head  | Fg123    | Mr89     | 3      |     | B              |
| 9/01/2019       | Fg123.6 | 0.0155     | flat head  | Fg123    | Mr89     | 3      |     | B              |
| 16/01/2019      | Fb148.5 | 0.014      | flat head  | Fb148    | Mg120    | 4      |     | A              |
| 16/01/2019      | Fb148.3 | 0.015      | flat head  | Fb148    | Mg120    | 4      |     | A              |
| 16/01/2019      | Fb148.4 | 0.015      | flat head  | Fb148    | Mg120    | 4      |     | A              |
| 16/01/2019      | Fb148.1 | 0.015      | flat head  | Fb148    | Mg120    | 4      |     | A              |
| 16/01/2019      | Fb148.7 | 0.014      | flat head  | Fb148    | Mg120    | 4      |     | A              |
| 16/01/2019      | Fb148.2 | 0.014      | flat head  | Fb148    | Mg120    | 4      |     | A              |
| 16/01/2019      | Fb148.6 | 0.014      | flat head  | Fb148    | Mg120    | 4      |     | A              |
| 24/09/2019      | Fw263.1 | 0.0164     | round head | Fw263    | Mw270    | 5      |     | A              |
| 24/09/2019      | Fw263.2 | 0.015      | round head | Fw263    | Mw270    | 5      |     | A              |
| 24/09/2019      | Fw263.3 | 0.0157     | round head | Fw263    | Mw270    | 5      |     | A              |
| 24/09/2019      | Fw263.4 | 0.0176     | round head | Fw263    | Mw270    | 5      |     | A              |
| 24/09/2019      | Fw263.5 | 0.0176     | round head | Fw263    | Mw270    | 5      |     | A              |
| 24/09/2019      | Fw263.6 | 0.0146     | round head | Fw263    | Mw270    | 5      |     | A              |
| 24/09/2019      | Fw263.7 | 0.0173     | round head | Fw263    | Mw270    | 5      |     | A              |
| 23/10/2019      | Fb197.1 | 0.0137     | flat head  | Fb197    | My234    | 6      | M   | B              |
| 23/10/2019      | Fb197.2 | 0.0112     | flat head  | Fb197    | My234    | 6      | M   | B              |
| 23/10/2019      | Fb197.3 | 0.0132     | flat head  | Fb197    | My234    | 6      | M   | B              |
| 23/10/2019      | Fb197.4 | 0.0114     | flat head  | Fb197    | My234    | 6      | M   | B              |
| 23/10/2019      | Fb197.5 | 0.0098     | flat head  | Fb197    | My234    | 6      | F   | B              |
| 23/10/2019      | Fb197.6 | 0.0121     | flat head  | Fb197    | My234    | 6      | F   | B              |
| 23/10/2019      | Fb197.7 | 0.0118     | flat head  | Fb197    | My234    | 6      | F   | B              |
| 23/10/2019      | Fb197.8 | 0.0116     | flat head  | Fb197    | My234    | 6      | F   | B              |
| 23/10/2019      | Fb197.9 | 0.0115     | flat head  | Fb197    | My234    | 6      | F   | B              |
| 29/10/2019      | Fb255.1 | 0.0186     | round head | Fb255    | My234    | 7      | F   | B              |
| 29/10/2019      | Fb255.2 | 0.0157     | round head | Fb255    | My234    | 7      | F   | B              |
| 29/10/2019      | Fb255.3 | 0.019      | round head | Fb255    | My234    | 7      | M   | B              |
| 29/10/2019      | Fb255.4 | 0.0187     | round head | Fb255    | My234    | 7      | F   | B              |
| 29/10/2019      | Fw264.1 | 0.0182     | round head | Fw264    | My234    | 8      | F   | A              |
| 29/10/2019      | Fw264.2 | 0.0185     | round head | Fw264    | My234    | 8      | F   | A              |
| 29/10/2019      | Fw264.3 | 0.0168     | round head | Fw264    | My234    | 8      | M   | A              |
| 29/10/2019      | Fw264.4 | 0.0187     | round head | Fw264    | My234    | 8      | M   | A              |
| 19/11/2019      | Fb148.1 | 0.017      | round head | Fb148    | My234    | 9      | F   | B              |
| 19/11/2019      | Fb148.2 | 0.018      | round head | Fb148    | My234    | 9      | M   | B              |
| 19/11/2019      | Fb148.3 | 0.019      | round head | Fb148    | My234    | 9      | F   | B              |
| 19/11/2019      | Fb148.4 | 0.018      | round head | Fb148    | My234    | 9      | M   | B              |
| 27/11/2019      | Fb198.1 | 0.013      | flat head  | Fb198    | My234    | 10     | F   | A              |
| 27/11/2019      | Fb198.2 | 0.016      | flat head  | Fb198    | My234    | 10     | F   | A              |
| 27/11/2019      | Fb198.3 | 0.012      | flat head  | Fb198    | My234    | 10     | F   | A              |
| 27/11/2019      | Fb198.5 | 0.015      | flat head  | Fb198    | My234    | 10     | M   | A              |

</br>

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

### Sequencing 

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

</p>
</details>

## Mouse ChIP-seq data
<details><summary>Experimental design</summary>
<p>
  
balhbalh
  
</p>
</details>

## Data management
```
chipseq-cross-species/
├── analysis
│   ├── dunnart_peaks
│   └── cross_species_peak_comparisons
├── data
│   ├── genomes
│   └── raw_reads
├── docs
├── code
│   ├── dunnart_peak_calling
│   ├── mouse_peak_calling
│   └── whole_genome_alignment
└── output
```
