# Extracted analysis commands for: Wild boar (Pig)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19

mkdir -p Genome_data/Pig
cd Genome_data/Pig

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/046/630/025/GCA_046630025.1_UWBC_WMS_v1/GCA_046630025.1_UWBC_WMS_v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/869/115/GCA_040869115.1_ASM4086911v1/GCA_040869115.1_ASM4086911v1_genomic.fna.gz

gunzip -k PIG_GCF_000003025.6_Sscrofa11.1_genomic_ref.fna.gz
gunzip -k N_PIG_GCA_040869115.1_ASM4086911v1_genomic.fna.gz
gunzip -k IN_PIG_GCA_046630025.1_UWBC_WMS_v1_genomic.fna.gz

###### Native Pig
minimap2 -ax asm5 -t 20 PIG_GCF_000003025.6_Sscrofa11.1_genomic_ref.fna N_PIG_GCA_040869115.1_ASM4086911v1_genomic.fna > N_Pig.sam
samtools view -@ 20 -bS N_Pig.sam | samtools sort -@ 20 -o N_Pig.sorted.bam
samtools index N_Pig.sorted.bam
bcftools mpileup -Ou -f PIG_GCF_000003025.6_Sscrofa11.1_genomic_ref.fna N_Pig.sorted.bam | bcftools call -mv -Ov -o N_Pig_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Pig_variants.vcf > N_Pig_filtered_variants.vcf

###### Invasive Pig
minimap2 -ax asm5 -t 20 PIG_GCF_000003025.6_Sscrofa11.1_genomic_ref.fna IN_PIG_GCA_046630025.1_UWBC_WMS_v1_genomic.fna > IN_Pig.sam
samtools view -@ 20 -bS IN_Pig.sam | samtools sort -@ 20 -o IN_Pig.sorted.bam
samtools index IN_Pig.sorted.bam
bcftools mpileup -Ou -f PIG_GCF_000003025.6_Sscrofa11.1_genomic_ref.fna IN_Pig.sorted.bam | bcftools call -mv -Ov -o IN_Pig_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Pig_variants.vcf > IN_Pig_filtered_variants.vcf

# To calculate nucleotide diversity (π), heterozygosity
bgzip IN_Pig.filtered.variants.vcf
bgzip N_Pig.filtered.variants.vcf
tabix -p vcf IN_Pig.filtered.variants.vcf.gz
tabix -p vcf N_Pig.filtered.variants.vcf.gz

# Step 1: Nucleotide Diversity (π) per status
vcftools --vcf N_Pig.filtered.variants.vcf.gz --window-pi 10000 --out N_Pig.10000_pi
vcftools --vcf IN_Pig.filtered.variants.vcf.gz --window-pi 10000 --out IN_Pig.10000_pi

# Step 2: Observed Heterozygosity per status
vcftools --vcf N_Pig.filtered.variants.vcf.gz --het --out N_Pig.het
vcftools --vcf IN_Pig.filtered.variants.vcf.gz --het --out IN_Pig.het

#############################
cd ..
