#!/usr/bin/env bash
# Extracted analysis commands for: Old World Bollworm (Old_World_Bollworm)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19

mkdir -p Genome_data/Old_World_Bollworm
cd Genome_data/Old_World_Bollworm

# Download genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/705/265/GCF_030705265.1_ASM3070526v1/GCF_030705265.1_ASM3070526v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/262/555/GCA_026262555.1_ASM2626255v1/GCA_026262555.1_ASM2626255v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/033/439/325/GCA_033439325.1_ASM3343932v1/GCA_033439325.1_ASM3343932v1_genomic.fna.gz

# Rename
mv -f GCF_030705265.1_ASM3070526v1_genomic.fna.gz BOLLWORM_GCF_030705265.1_ASM3070526v1_genomic_ref.fna.gz
mv -f GCA_026262555.1_ASM2626255v1_genomic.fna.gz N_BOLLWORM_GCA_026262555.1_ASM2626255v1_genomic.fna.gz
mv -f GCA_033439325.1_ASM3343932v1_genomic.fna.gz IN_BOLLWORM_GCA_033439325.1_ASM3343932v1_genomic.fna.gz

cd Old_World_Bollworm
gunzip -k BOLLWORM_GCF_030705265.1_ASM3070526v1_genomic_ref.fna.gz
gunzip -k N_BOLLWORM_GCA_026262555.1_ASM2626255v1_genomic.fna.gz
gunzip -k IN_BOLLWORM_GCA_033439325.1_ASM3343932v1_genomic.fna.gz

###### native
minimap2 -ax asm5 -t 20 BOLLWORM_GCF_030705265.1_ASM3070526v1_genomic_ref.fna N_BOLLWORM_GCA_026262555.1_ASM2626255v1_genomic.fna > N_Old_World_Bollworm.sam
samtools view -@ 20 -bS N_Old_World_Bollworm.sam | samtools sort -@ 20 -o N_Old_World_Bollworm.sorted.bam
samtools index N_Old_World_Bollworm.sorted.bam
bcftools mpileup -Ou -f BOLLWORM_GCF_030705265.1_ASM3070526v1_genomic_ref.fna N_Old_World_Bollworm.sorted.bam | bcftools call -mv -Ov -o N_Old_World_Bollworm.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Old_World_Bollworm.variants.vcf > N_Old_World_Bollworm.filtered.variants.vcf

###### invasive
minimap2 -ax asm5 -t 20 BOLLWORM_GCF_030705265.1_ASM3070526v1_genomic_ref.fna IN_BOLLWORM_GCA_033439325.1_ASM3343932v1_genomic.fna > IN_Old_World_Bollworm.sam
samtools view -@ 20 -bS IN_Old_World_Bollworm.sam | samtools sort -@ 20 -o IN_Old_World_Bollworm.sorted.bam
samtools index IN_Old_World_Bollworm.sorted.bam
bcftools mpileup -Ou -f BOLLWORM_GCF_030705265.1_ASM3070526v1_genomic_ref.fna IN_Old_World_Bollworm.sorted.bam | bcftools call -mv -Ov -o IN_Old_World_Bollworm.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Old_World_Bollworm.variants.vcf > IN_Old_World_Bollworm.filtered.variants.vcf

# Post-processing and diversity metrics
bgzip IN_Old_World_Bollworm.filtered.variants.vcf
bgzip N_Old_World_Bollworm.filtered.variants.vcf
tabix -p vcf IN_Old_World_Bollworm.filtered.variants.vcf.gz
tabix -p vcf N_Old_World_Bollworm.filtered.variants.vcf.gz

# Step 1: Nucleotide Diversity (π) per status
vcftools --vcf N_Old_World_Bollworm.filtered.variants.vcf.gz --window-pi 10000 --out N_Old_World_Bollworm.10000_pi
vcftools --vcf IN_Old_World_Bollworm.filtered.variants.vcf.gz --window-pi 10000 --out IN_Old_World_Bollworm.10000_pi

# Step 2: Observed Heterozygosity per status
vcftools --vcf N_Old_World_Bollworm.filtered.variants.vcf.gz --het --out N_Old_World_Bollworm.het
vcftools --vcf IN_Old_World_Bollworm.filtered.variants.vcf.gz --het --out IN_Old_World_Bollworm.het

#############################
cd ..
