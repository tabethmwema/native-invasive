# Extracted analysis commands for: Chinese mitten crab (Crab)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19

mkdir -p Genome_data/Crab
cd Genome_data/Crab
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/024/679/095/GCF_024679095.1_ASM2467909v1/GCF_024679095.1_ASM2467909v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/436/485/GCA_013436485.1_ASM1343648v1/GCA_013436485.1_ASM1343648v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/348/025/GCA_023348025.1_ASM2334802v1/GCA_023348025.1_ASM2334802v1_genomic.fna.gz

gunzip -k CRAB_GCF_024679095.1_ASM2467909v1_genomic_ref.fna.gz
gunzip -k N_CRAB_GCA_023348025.1_ASM2334802v1_genomic.fna.gz
gunzip -k IN_CRAB_GCA_013436485.1_ASM1343648v1_genomic.fna.gz

###### Native Crabs
minimap2 -ax asm5 -t 20 CRAB_GCF_024679095.1_ASM2467909v1_genomic_ref.fna N_CRAB_GCA_023348025.1_ASM2334802v1_genomic.fna > N_Crab.sam
samtools view -@ 20 -bS N_Crab.sam | samtools sort -@ 20 -o N_Crab.sorted.bam
samtools index N_Crab.sorted.bam
bcftools mpileup -Ou -f CRAB_GCF_024679095.1_ASM2467909v1_genomic_ref.fna N_Crab.sorted.bam | bcftools call -mv -Ov -o N_Crab_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Crab_variants.vcf > N_Crab_filtered_variants.vcf

###### Invasive Crabs
minimap2 -ax asm5 -t 20 CRAB_GCF_024679095.1_ASM2467909v1_genomic_ref.fna IN_CRAB_GCA_013436485.1_ASM1343648v1_genomic.fna > IN_Crab.sam
samtools view -@ 20 -bS IN_Crab.sam | samtools sort -@ 20 -o IN_Crab.sorted.bam
samtools index IN_Crab.sorted.bam
bcftools mpileup -Ou -f CRAB_GCF_024679095.1_ASM2467909v1_genomic_ref.fna IN_Crab.sorted.bam | bcftools call -mv -Ov -o IN_Crab_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Crab_variants.vcf > IN_Crab_filtered_variants.vcf

# To calculate nucleotide diversity (π), heterozygosity
bgzip IN_Crab_filtered_variants.vcf
bgzip N_Crab_filtered_variants.vcf
tabix -p vcf IN_Crab_filtered_variants.vcf.gz
tabix -p vcf N_Crab_filtered_variants.vcf.gz

# Step 1: Calculate Nucleotide Diversity (π) per status
vcftools --vcf N_Crab.filtered.variants.vcf --window-pi 10000 --out N_Crab.10000_pi
vcftools --vcf IN_Crab.filtered.variants.vcf --window-pi 10000 --out IN_Crab.10000_pi

# Step 2: Calculate Observed Heterozygosity per status
vcftools --vcf N_Crab.filtered.variants.vcf --het --out N_Crab.het
vcftools --vcf IN_Crab.filtered.variants.vcf --het --out IN_Crab.het

#############################
cd ..