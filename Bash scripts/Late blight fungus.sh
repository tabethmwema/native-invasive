# Source: Raw_script_script_updated_V2
# Extracted analysis commands for: Late blight fungus (Fungus)
# Extracted on: 2025-08-19

mkdir -p Genome_data/Fungus
cd Genome_data/Fungus
gunzip -k FUNGUS_GCF_000142945.1_ASM14294v1_genomic_ref.fna.gz
gunzip -k IN_FUNGUS_GCA_044509515.1_PHIN_VZR-VIR21a_genomic.fna.gz
gunzip -k N_FUNGUS_GCA_026225685.1_Pinf1306_UCR_1_genomic.fna.gz

######  Native Fungus
minimap2 -ax asm5 -t 20 FUNGUS_GCF_000142945.1_ASM14294v1_genomic_ref.fna N_FUNGUS_GCA_026225685.1_Pinf1306_UCR_1_genomic.fna > N_Fungus.sam
samtools view -@ 20 -bS N_Fungus.sam | samtools sort -@ 20 -o N_Fungus.sorted.bam
samtools index N_Fungus.sorted.bam
bcftools mpileup -Ou -f FUNGUS_GCF_000142945.1_ASM14294v1_genomic_ref.fna N_Fungus.sorted.bam | bcftools call -mv -Ov -o N_Fungus_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Fungus_variants.vcf > N_Fungus_filtered_variants.vcf

###### Invasive Fungus
minimap2 -ax asm5 -t 20 FUNGUS_GCF_000142945.1_ASM14294v1_genomic_ref.fna IN_FUNGUS_GCA_044509515.1_PHIN_VZR-VIR21a_genomic.fna > IN_Fungus.sam
samtools view -@ 20 -bS IN_Fungus.sam | samtools sort -@ 20 -o IN_Fungus.sorted.bam
samtools index IN_Fungus.sorted.bam
bcftools mpileup -Ou -f FUNGUS_GCF_000142945.1_ASM14294v1_genomic_ref.fna IN_Fungus.sorted.bam | bcftools call -mv -Ov -o IN_Fungus_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Fungus_variants.vcf > IN_Fungus_filtered_variants.vcf
#############################__

#To calculate nucleotide diversity (π), heterozygosity
bgzip IN_Fungus.filtered.variants.vcf
bgzip N_Fungus.filtered.variants.vcf
tabix -p vcf IN_Fungus.filtered.variants.vcf
tabix -p vcf N_Fungus.filtered.variants.vcf

# Step 1: Calculate Nucleotide Diversity (π) per status
vcftools --vcf N_Fungus.filtered.variants.vcf --window-pi 10000 --out N_Fungus.10000_pi
vcftools --vcf IN_Fungus.filtered.variants.vcf --window-pi 10000 --out IN_Fungus.10000_pi

# Step 2: Calculate Observed Heterozygosity per status
vcftools --vcf N_Fungus.filtered.variants.vcf --het --out N_Fungus.het
vcftools --vcf IN_Fungus.filtered.variants.vcf --het --out IN_Fungus.het

#############################
cd ..
