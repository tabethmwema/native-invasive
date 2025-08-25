# Extracted analysis commands for: Goldfish (Fish)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19

mkdir -p Genome_data/Fish
cd Genome_data/Fish
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/368/295/GCF_003368295.1_ASM336829v1/GCF_003368295.1_ASM336829v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/037/022/325/GCA_037022325.1_ASM3702232v1/GCA_037022325.1_ASM3702232v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/720/715/GCA_019720715.2_ASM1972071v2/GCA_019720715.2_ASM1972071v2_genomic.fna.gz
cd Fish
gunzip -k FISH_GCF_003368295.1_ASM336829v1_genomic_ref.fna.gz
gunzip -k N_FISH_GCA_019720715.2_ASM1972071v2_genomic.fna.gz
gunzip -k IN_FISH_GCA_037022325.1_ASM3702232v1_genomic.fna.gz

###### Native Fish
minimap2 -ax asm5 -t 20 FISH_GCF_003368295.1_ASM336829v1_genomic_ref.fna N_FISH_GCA_019720715.2_ASM1972071v2_genomic.fna > N_Fish.sam
samtools view -@ 20 -bS N_Fish.sam | samtools sort -@ 20 -o N_Fish.sorted.bam
samtools index N_Fish.sorted.bam
bcftools mpileup -Ou -f FISH_GCF_003368295.1_ASM336829v1_genomic_ref.fna N_Fish.sorted.bam | bcftools call -mv -Ov -o N_Fish_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Fish_variants.vcf > N_Fish_filtered_variants.vcf

###### Invasive Fish
minimap2 -ax asm5 -t 20 FISH_GCF_003368295.1_ASM336829v1_genomic_ref.fna IN_FISH_GCA_037022325.1_ASM3702232v1_genomic.fna > IN_Fish.sam
samtools view -@ 20 -bS IN_Fish.sam | samtools sort -@ 20 -o IN_Fish.sorted.bam
samtools index IN_Fish.sorted.bam
bcftools mpileup -Ou -f FISH_GCF_003368295.1_ASM336829v1_genomic_ref.fna IN_Fish.sorted.bam | bcftools call -mv -Ov -o IN_Fish_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Fish_variants.vcf > IN_Fish_filtered_variants.vcf

# To calculate nucleotide diversity (π), heterozygosity
bgzip IN_Fish_filtered_variants.vcf
bgzip N_Fish_filtered_variants.vcf
tabix -p vcf IN_Fish_filtered_variants.vcf.gz
tabix -p vcf N_Fish_filtered_variants.vcf.gz

# Step 1: Calculate Nucleotide Diversity (π) per status
vcftools --vcf N_Fish_filtered_variants.vcf --window-pi 10000 --out N_Fish.10000_pi
vcftools --vcf IN_Fish_filtered_variants.vcf --window-pi 10000 --out IN_Fish.10000_pi

# Step 2: Calculate Observed Heterozygosity per status
vcftools --vcf N_Fish_filtered_variants.vcf --het --out N_Fish.het
vcftools --vcf IN_Fish_filtered_variants.vcf --het --out IN_Fish.het

#############################
cd ..
