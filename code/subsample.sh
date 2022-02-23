## Command to subsample bam files

TRA=($(for file in B*_PPq30.sorted.dedup.bam; do echo $file |cut -d "_" -f 1-2;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}

frac=$( samtools idxstats ${tr}_PPq30.sorted.dedup.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=10000000/total; if (frac > 1) {print 1} else {print frac}}' )

samtools view -bs $frac ${tr}_PPq30.sorted.dedup.bam > ../../results_10M/bowtie2/${tr}_PPq30.sorted.dedup.bam

done