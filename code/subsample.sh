## Command to subsample bam files

TRA=($(for file in *_downSampled_q30.sorted.dedup.bam; do echo $file |cut -d "_" -f 1-3;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do
echo ${tr}

samtools index ${tr}_downSampled_q30.sorted.dedup.bam ${tr}_downSampled_q30.sorted.dedup.bai

#frac=$( samtools idxstats ${tr}_PPq30.sorted.dedup.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=10000000/total; if (frac > 1) {print 1} else {print frac}}' )

#samtools view -bs $frac ${tr}_PPq30.sorted.dedup.bam > ${tr}_downSampled_PPq30.sorted.dedup.bam

done
