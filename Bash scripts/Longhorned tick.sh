#!/usr/bin/env bash
# Extracted analysis commands for: Longhorned tick (Ticks)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19

mkdir -p Genome_data/Ticks
cd Genome_data/Ticks

# Download genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/122/185/GCA_008122185.1_HLAgriLifeRun1/GCA_008122185.1_HLAgriLifeRun1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/414/705/GCA_022414705.1_ASM2241470v1/GCA_022414705.1_ASM2241470v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/339/765/GCA_013339765.2_BIME_HaeL_1.3/GCA_013339765.2_BIME_HaeL_1.3_genomic.fna.gz

# Rename to match pattern used in other species blocks
mv GCA_013339765.2_BIME_HaeL_1.3_genomic.fna.gz TICKS_GCA_013339765.2_BIME_HaeL_1.3_genomic_ref.fna.gz
mv GCA_022414705.1_ASM2241470v1_genomic.fna.gz N_TICKS_GCA_022414705.1_ASM2241470v1_genomic.fna.gz
mv GCA_008122185.1_HLAgriLifeRun1_genomic.fna.gz IN_TICKS_GCA_008122185.1_HLAgriLifeRun1_genomic.fna.gz

# Prepare FASTA (keep originals)
gunzip -k TICKS_GCA_013339765.2_BIME_HaeL_1.3_genomic_ref.fna.gz
gunzip -k N_TICKS_GCA_022414705.1_ASM2241470v1_genomic.fna.gz
gunzip -k IN_TICKS_GCA_008122185.1_HLAgriLifeRun1_genomic.fna.gz

###### native
minimap2 -ax asm5 -t 20 TICKS_GCA_013339765.2_BIME_HaeL_1.3_genomic_ref.fna N_TICKS_GCA_022414705.1_ASM2241470v1_genomic.fna > N_Ticks.sam
samtools view -@ 20 -bS N_Ticks.sam | samtools sort -@ 20 -o N_Ticks.sorted.bam
samtools index N_Ticks.sorted.bam
bcftools mpileup -Ou -f TICKS_GCA_013339765.2_BIME_HaeL_1.3_genomic_ref.fna N_Ticks.sorted.bam | bcftools call -mv -Ov -o N_Ticks.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Ticks.variants.vcf > N_Ticks.filtered.variants.vcf

###### invasive
minimap2 -ax asm5 -t 20 TICKS_GCA_013339765.2_BIME_HaeL_1.3_genomic_ref.fna IN_TICKS_GCA_008122185.1_HLAgriLifeRun1_genomic.fna > IN_Ticks.sam
samtools view -@ 20 -bS IN_Ticks.sam | samtools sort -@ 20 -o IN_Ticks.sorted.bam
samtools index IN_Ticks.sorted.bam
bcftools mpileup -Ou -f TICKS_GCA_013339765.2_BIME_HaeL_1.3_genomic_ref.fna IN_Ticks.sorted.bam | bcftools call -mv -Ov -o IN_Ticks.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Ticks.variants.vcf > IN_Ticks.filtered.variants.vcf

# To calculate nucleotide diversity (π), heterozygosity
bgzip IN_Ticks.filtered.variants.vcf
bgzip N_Ticks.filtered.variants.vcf
tabix -p vcf IN_Ticks.filtered.variants.vcf.gz
tabix -p vcf N_Ticks.filtered.variants.vcf.gz

# Step 1: Calculate Nucleotide Diversity (π) per status
vcftools --vcf N_Ticks.filtered.variants.vcf.gz --window-pi 10000 --out N_Ticks.10000_pi
vcftools --vcf IN_Ticks.filtered.variants.vcf.gz --window-pi 10000 --out IN_Ticks.10000_pi

# Step 2: Calculate Observed Heterozygosity per status
vcftools --vcf N_Ticks.filtered.variants.vcf.gz --het --out N_Ticks.het
vcftools --vcf IN_Ticks.filtered.variants.vcf.gz --het --out IN_Ticks.het

#############################
cd ..
