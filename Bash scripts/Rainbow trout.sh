#!/usr/bin/env bash
# Extracted analysis commands for: Rainbow trout (RainbowFish)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19

mkdir -p Genome_data/RainbowFish
cd Genome_data/RainbowFish
# Download genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/265/735/GCF_013265735.2_USDA_OmykA_1.1/GCF_013265735.2_USDA_OmykA_1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/034/753/235/GCA_034753235.1_USDA_OmykKC_1.0/GCA_034753235.1_USDA_OmykKC_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/024/535/GCA_036024535.1_ASM3602453v1/GCA_036024535.1_ASM3602453v1_genomic.fna.gz

# Rename to match the pattern used across species
mv -f GCF_013265735.2_USDA_OmykA_1.1_genomic.fna.gz RAINBOWFISH_GCF_013265735.2_USDA_OmykA_1.1_genomic_ref.fna.gz
mv -f GCA_034753235.1_USDA_OmykKC_1.0_genomic.fna.gz N_RAINBOWFISH_GCA_034753235.1_USDA_OmykKC_1.0_genomic.fna.gz
mv -f GCA_036024535.1_ASM3602453v1_genomic.fna.gz IN_RAINBOWFISH_GCA_036024535.1_ASM3602453v1_genomic.fna.gz

gunzip -k RAINBOWFISH_GCF_013265735.2_USDA_OmykA_1.1_genomic_ref.fna.gz
gunzip -k N_RAINBOWFISH_GCA_034753235.1_USDA_OmykKC_1.0_genomic.fna.gz
gunzip -k IN_RAINBOWFISH_GCA_036024535.1_ASM3602453v1_genomic.fna.gz

###### native
minimap2 -ax asm5 -t 20 RAINBOWFISH_GCF_013265735.2_USDA_OmykA_1.1_genomic_ref.fna N_RAINBOWFISH_GCA_034753235.1_USDA_OmykKC_1.0_genomic.fna > N_RainbowFish.sam
samtools view -@ 20 -bS N_RainbowFish.sam | samtools sort -@ 20 -o N_RainbowFish.sorted.bam
samtools index N_RainbowFish.sorted.bam
bcftools mpileup -Ou -f RAINBOWFISH_GCF_013265735.2_USDA_OmykA_1.1_genomic_ref.fna N_RainbowFish.sorted.bam | bcftools call -mv -Ov -o N_RainbowFish.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_RainbowFish.variants.vcf > N_RainbowFish.filtered.variants.vcf

###### invasive
minimap2 -ax asm5 -t 20 RAINBOWFISH_GCF_013265735.2_USDA_OmykA_1.1_genomic_ref.fna IN_RAINBOWFISH_GCA_036024535.1_ASM3602453v1_genomic.fna > IN_RainbowFish.sam
samtools view -@ 20 -bS IN_RainbowFish.sam | samtools sort -@ 20 -o IN_RainbowFish.sorted.bam
samtools index IN_RainbowFish.sorted.bam
bcftools mpileup -Ou -f RAINBOWFISH_GCF_013265735.2_USDA_OmykA_1.1_genomic_ref.fna IN_RainbowFish.sorted.bam | bcftools call -mv -Ov -o IN_RainbowFish.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_RainbowFish.variants.vcf > IN_RainbowFish.filtered.variants.vcf

# Post-processing and diversity metrics
bgzip IN_RainbowFish.filtered.variants.vcf
bgzip N_RainbowFish.filtered.variants.vcf
tabix -p vcf IN_RainbowFish.filtered.variants.vcf.gz
tabix -p vcf N_RainbowFish.filtered.variants.vcf.gz

# Step 1: Nucleotide Diversity (π) per status
vcftools --vcf N_RainbowFish.filtered.variants.vcf.gz --window-pi 10000 --out N_RainbowFish.10000_pi
vcftools --vcf IN_RainbowFish.filtered.variants.vcf.gz --window-pi 10000 --out IN_RainbowFish.10000_pi

# Step 2: Observed Heterozygosity per status
vcftools --vcf N_RainbowFish.filtered.variants.vcf.gz --het --out N_RainbowFish.het
vcftools --vcf IN_RainbowFish.filtered.variants.vcf.gz --het --out IN_RainbowFish.het

#############################
cd ..
