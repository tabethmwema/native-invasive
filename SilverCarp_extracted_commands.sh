#!/usr/bin/env bash
# Extracted analysis commands for: Silver carp (SilverCarp)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19

mkdir -p Genome_data/SilverCarp
cd Genome_data/SilverCarp

# Download genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/037/950/675/GCA_037950675.1_Hypophthalmichthys_molitrix_1.0/GCA_037950675.1_Hypophthalmichthys_molitrix_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/041/475/455/GCA_041475455.1_Hmo_1.0/GCA_041475455.1_Hmo_1.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/764/525/GCA_004764525.1_HypMol1.0/GCA_004764525.1_HypMol1.0_genomic.fna.gz

# Rename
mv -f GCA_037950675.1_Hypophthalmichthys_molitrix_1.0_genomic.fna.gz SILVERCARP_GCA_037950675.1_Hypophthalmichthys_molitrix_1.0_genomic_ref.fna.gz
mv -f GCA_041475455.1_Hmo_1.0_genomic.fna.gz N_SILVERCARP_GCA_041475455.1_Hmo_1.0_genomic.fna.gz
mv -f GCA_004764525.1_HypMol1.0_genomic.fna.gz IN_SILVERCARP_GCA_004764525.1_HypMol1.0_genomic.fna.gz

cd SilverCarp
gunzip -k SILVERCARP_GCA_037950675.1_Hypophthalmichthys_molitrix_1.0_genomic_ref.fna.gz
gunzip -k N_SILVERCARP_GCA_041475455.1_Hmo_1.0_genomic.fna.gz
gunzip -k IN_SILVERCARP_GCA_004764525.1_HypMol1.0_genomic.fna.gz

###### native
minimap2 -ax asm5 -t 20 SILVERCARP_GCA_037950675.1_Hypophthalmichthys_molitrix_1.0_genomic_ref.fna N_SILVERCARP_GCA_041475455.1_Hmo_1.0_genomic.fna > N_SilverCarp.sam
samtools view -@ 20 -bS N_SilverCarp.sam | samtools sort -@ 20 -o N_SilverCarp.sorted.bam
samtools index N_SilverCarp.sorted.bam
bcftools mpileup -Ou -f SILVERCARP_GCA_037950675.1_Hypophthalmichthys_molitrix_1.0_genomic_ref.fna N_SilverCarp.sorted.bam | bcftools call -mv -Ov -o N_SilverCarp.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_SilverCarp.variants.vcf > N_SilverCarp.filtered.variants.vcf

###### invasive
minimap2 -ax asm5 -t 20 SILVERCARP_GCA_037950675.1_Hypophthalmichthys_molitrix_1.0_genomic_ref.fna IN_SILVERCARP_GCA_004764525.1_HypMol1.0_genomic.fna > IN_SilverCarp.sam
samtools view -@ 20 -bS IN_SilverCarp.sam | samtools sort -@ 20 -o IN_SilverCarp.sorted.bam
samtools index IN_SilverCarp.sorted.bam
bcftools mpileup -Ou -f SILVERCARP_GCA_037950675.1_Hypophthalmichthys_molitrix_1.0_genomic_ref.fna IN_SilverCarp.sorted.bam | bcftools call -mv -Ov -o IN_SilverCarp.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_SilverCarp.variants.vcf > IN_SilverCarp.filtered.variants.vcf

# Post-processing and diversity metrics
bgzip IN_SilverCarp.filtered.variants.vcf
bgzip N_SilverCarp.filtered.variants.vcf
tabix -p vcf IN_SilverCarp.filtered.variants.vcf.gz
tabix -p vcf N_SilverCarp.filtered.variants.vcf.gz

# Step 1: Nucleotide Diversity (π) per status
vcftools --vcf N_SilverCarp.filtered.variants.vcf.gz --window-pi 10000 --out N_SilverCarp.10000_pi
vcftools --vcf IN_SilverCarp.filtered.variants.vcf.gz --window-pi 10000 --out IN_SilverCarp.10000_pi

# Step 2: Observed Heterozygosity per status
vcftools --vcf N_SilverCarp.filtered.variants.vcf.gz --het --out N_SilverCarp.het
vcftools --vcf IN_SilverCarp.filtered.variants.vcf.gz --het --out IN_SilverCarp.het

#############################__
cd ..
