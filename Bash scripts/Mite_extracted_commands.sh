# Extracted analysis commands for: Two-spotted spider mite (Mite)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19S

mkdir -p Genome_data/Mite
cd  Genome_data/Mite
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/239/435/GCF_000239435.1_ASM23943v1/GCF_000239435.1_ASM23943v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/877/765/GCA_036877765.1_tu.genome.baafs_v2/GCA_036877765.1_tu.genome.baafs_v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/522/325/GCA_030522325.1_ASM3052232v1/GCA_030522325.1_ASM3052232v1_genomic.fna.gz
cd Mite
gunzip -k MITE_GCF_000239435.1_ASM23943v1_genomic_ref.fna.gz
gunzip -k N_MITE_GCA_030522325.1_ASM3052232v1_genomic.fna.gz
gunzip -k IN_MITE_GCA_036877765.1_tu.genome.baafs_v2_genomic.fna.gz
######  Native Mite
minimap2 -ax asm5 -t 20 MITE_GCF_000239435.1_ASM23943v1_genomic_ref.fna N_MITE_GCA_030522325.1_ASM3052232v1_genomic.fna > N_Mite.sam
samtools view -@ 20 -bS N_Mite.sam | samtools sort -@ 20 -o N_Mite.sorted.bam
samtools index N_Mite.sorted.bam
bcftools mpileup -Ou -f MITE_GCF_000239435.1_ASM23943v1_genomic_ref.fna N_Mite.sorted.bam | bcftools call -mv -Ov -o N_Mite_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Mite_variants.vcf > N_Mite_filtered_variants.vcf
###### Invasive Mite
minimap2 -ax asm5 -t 20 MITE_GCF_000239435.1_ASM23943v1_genomic_ref.fna IN_MITE_GCA_036877765.1_tu.genome.baafs_v2_genomic.fna > IN_Mite.sam
samtools view -@ 20 -bS IN_Mite.sam | samtools sort -@ 20 -o IN_Mite.sorted.bam
samtools index IN_Mite.sorted.bam
bcftools mpileup -Ou -f MITE_GCF_000239435.1_ASM23943v1_genomic_ref.fna IN_Mite.sorted.bam | bcftools call -mv -Ov -o IN_Mite_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Mite_variants.vcf > IN_Mite_filtered_variants.vcf

# To calculate nucleotide diversity (π), heterozygosity
bgzip IN_Mite.filtered.variants.vcf
bgzip N_Mite.filtered.variants.vcf
tabix -p vcf IN_Mite.filtered.variants.vcf.gz
tabix -p vcf N_Mite.filtered.variants.vcf.gz

# Step 1: Calculate Nucleotide Diversity (π) per status
vcftools --vcf N_Mite.filtered.variants.vcf --window-pi 10000 --out N_Mite.10000_pi
vcftools --vcf IN_Mite.filtered.variants.vcf --window-pi 10000 --out IN_Mite.10000_pi

# Step 2: Calculate Observed Heterozygosity per status
vcftools --vcf N_Mite.filtered.variants.vcf --het --out N_Mite.het
vcftools --vcf IN_Mite.filtered.variants.vcf --het --out IN_Mite.het

cd ..
