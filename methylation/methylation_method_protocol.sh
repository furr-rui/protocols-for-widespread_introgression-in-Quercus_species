########################################################################################################################
# Documents and scripts were written by: Ruirui Fu & Yuxiang Zhu
# For manuscript: Yuxiang Zhu,Ruirui Fu et al. (2026). "Widespread introgression across a phylogeny of 168 genomes from 67 Quercus species and its association to gene expression regulation"
# email: yxzhu@zju.edu.cn,furuirui@zju.edu.cn
# Jun Chen Lab. 
########################################################################################################################


########################################################################################################################
# Methylation analysis of 19 species of the Quercus 
########################################################################################################################
#############Methylation site extraction
## step1  Read Mapping (Loop processing can be done using scripts methylation_alignment_flexible.pl and flt_bam.pl)
#Creating the Reference Index
biscuit index ref_file
#Aligning Reads to the Reference
biscuit align -@ 40 -R \"\@RG\\tID:$_\\tSM:$_\" $ref_file ./$fq1 ./$fq2 | dupsifter $ref_file | samtools sort -m 4g -@ 30 -o $out_path/$_\'_srt.markdup.bam\' -O BAM -
#Filter BAM files to retain unique properly paired
samtools view -q 20 -f 0x0002 -F 0X0004 -F 0X0008 -b sample_srt.markdup.bam > sample_srt_flt.bam
samtools view -h ./sample_srt_flt.bam | awk '$1~/^@/ || $7=="=" {print}' | samtools view -b - > sample_srt_flt_chr.bam
## step2  Extract DNA methylation 
#Generate Standard VCF Output
samtools index sample_srt_flt_chr.bam
biscuit pileup -@ 20 -o sample.vcf $ref_file sample_srt_flt_chr.bam
#Compresses and indexes the VCF
bgzip -@ 20 sample.vcf
tabix -p vcf sample.vcf.gz
#Extract DNA methylation into BED format
#Also compresses and indexes the BED
biscuit vcf2bed -e -t c sample.vcf.gz sample.methylation_tce_alltype.bed
bgzip sample.methylation_tce_alltype.bed
tabix -p bed sample.methylation_tce_alltype.bed.gz
#############Methylation down analysis
#Script for plotting methylation correlations among individuals
run.spe19_corrletion_plot.R
#Scripts for mapping the different regions of introgression gene and non-introgression in 19 species
methylation_lineplot_spe19_mean.mean.submit.R