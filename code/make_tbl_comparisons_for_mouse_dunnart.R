
system("bedtools intersect -a \
output/liftover/dunnart_promoter_50bpsummits_annotation_smiCraTOmm10.bed \
-b output/filtered_peaks/E15_cluster1_peaks.narrowPeak -wo > \
output/liftover/E15.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed")
system("bedtools intersect -a \
output/liftover/dunnart_promoter_50bpsummits_annotation_smiCraTOmm10.bed \
-b output/filtered_peaks/E14_cluster1_peaks.narrowPeak -wo > \
output/liftover/E14.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed")
system("bedtools intersect -a \
output/liftover/dunnart_promoter_50bpsummits_annotation_smiCraTOmm10.bed \
-b output/filtered_peaks/E13_cluster1_peaks.narrowPeak -wo > \
output/liftover/E13.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed")
system("bedtools intersect -a \
output/liftover/dunnart_promoter_50bpsummits_annotation_smiCraTOmm10.bed \
-b output/filtered_peaks/E12_cluster1_peaks.narrowPeak -wo > \
output/liftover/E12.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed")
system("bedtools intersect -a \
output/liftover/dunnart_promoter_50bpsummits_annotation_smiCraTOmm10.bed \
-b output/filtered_peaks/E11_cluster1_peaks.narrowPeak -wo > \
output/liftover/E11.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed")
system("bedtools intersect -a \
output/liftover/dunnart_promoter_50bpsummits_annotation_smiCraTOmm10.bed \
-b output/filtered_peaks/E10_cluster1_peaks.narrowPeak -wo > \
output/liftover/E10.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed")

system("liftOver -bedPlus=4 \
output/liftover/E10.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed \
output/wga/liftOver_chains/mm10TOsmiCra1_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain \
output/liftover/E10.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra.bed \
output/liftover/E10.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra_unmapped.bed")

system("liftOver -bedPlus=4 \
output/liftover/E11.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed \
output/wga/liftOver_chains/mm10TOsmiCra1_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain \
output/liftover/E11.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra.bed \
output/liftover/E11.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra_unmapped.bed")

system("liftOver -bedPlus=4 \
output/liftover/E12.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed \
output/wga/liftOver_chains/mm10TOsmiCra1_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain \
output/liftover/E12.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra.bed \
output/liftover/E12.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra_unmapped.bed")

system("liftOver -bedPlus=4 \
output/liftover/E13.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed \
output/wga/liftOver_chains/mm10TOsmiCra1_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain \
output/liftover/E13.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra.bed \
output/liftover/E13.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra_unmapped.bed")

system("liftOver -bedPlus=4 \
output/liftover/E14.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed \
output/wga/liftOver_chains/mm10TOsmiCra1_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain \
output/liftover/E14.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra.bed \
output/liftover/E14.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra_unmapped.bed")

system("liftOver -bedPlus=4 \
output/liftover/E15.5_mm10intersect_dunnart_promoter_50bpsummits_annotated.bed \
output/wga/liftOver_chains/mm10TOsmiCra1_patched_sorted_repFil_chainCl_preNet_chainNet_stitched.chain \
output/liftover/E15.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra.bed \
output/liftover/E15.5_mm10intersect_dunnart_promoter_50bpsummits_annotated_smiCraTOmm10_mm10TOsmiCra_unmapped.bed")

# dunnartscaff\tstart\tstop\tmousescaff\tstart\tstop\targetgene

