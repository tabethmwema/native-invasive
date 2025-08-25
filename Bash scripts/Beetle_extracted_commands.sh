# Extracted analysis commands for: Colorado potato beetle (Beetle)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19

mkdir -p Genome_data/Beetle
cd Genome_data/Beetle
###### Download genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/532/075/GCA_025532075.1_Ldec_3_OR/GCA_025532075.1_Ldec_3_OR_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/712/935/GCA_024712935.1_ASM2471293v1/GCA_024712935.1_ASM2471293v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/500/325/GCF_000500325.1_Ldec_2.0/GCF_000500325.1_Ldec_2.0_genomic.fna.gz

###### Decompress but keep originals
gunzip -k BEETLE_GCF_000500325.1_Ldec_2.0_genomic_ref.fna.gz
gunzip -k IN_BEETLE_GCA_024712935.1_ASM2471293v1_genomic.fna.gz
gunzip -k N_BEETLE_GCA_025532075.1_Ldec_3_OR_genomic.fna.gz

###### Align native & invasive genomes to the reference
###### native
minimap2 -ax asm5 -t 20  BEETLE_GCF_000500325.1_Ldec_2.0_genomic_ref.fna N_BEETLE_GCA_025532075.1_Ldec_3_OR_genomic.fna > N_Beetle.sam
samtools view -@ 20 -bS N_Beetle.sam | samtools sort -@ 20 -o N_Beetle.sorted.bam
samtools index N_Beetle.sorted.bam
bcftools mpileup -Ou -f BEETLE_GCF_000500325.1_Ldec_2.0_genomic_ref.fna N_Beetle.sorted.bam | bcftools call -mv -Ov -o N_Beetle.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Beetle.variants.vcf > N_Beetle.filtered.variants.vcf

###### invasive
minimap2 -ax asm5 -t 20  BEETLE_GCF_000500325.1_Ldec_2.0_genomic_ref.fna IN_BEETLE_GCA_024712935.1_ASM2471293v1_genomic.fna > IN_Beetle.sam
samtools view -@ 20 -bS IN_Beetle.sam | samtools sort -@ 20 -o IN_Beetle.sorted.bam
samtools index IN_Beetle.sorted.bam
bcftools mpileup -Ou -f BEETLE_GCF_000500325.1_Ldec_2.0_genomic_ref.fna IN_Beetle.sorted.bam | bcftools call -mv -Ov -o IN_Beetle.variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Beetle.variants.vcf > IN_Beetle.filtered.variants.vcf

###### Calculate nucleotide diversity (π), heterozygosity
bgzip IN_Beetle.filtered.variants.vcf
bgzip N_Beetle.filtered.variants.vcf
tabix -p vcf IN_Beetle.filtered.variants.vcf.gz
tabix -p vcf N_Beetle.filtered.variants.vcf.gz

# Step 1: Calculate Nucleotide Diversity (π) per status
vcftools --vcf N_Beetle.filtered.variants.vcf --window-pi 10000 --out N_Beetle.10000_pi
vcftools --vcf IN_Beetle.filtered.variants.vcf --window-pi 10000 --out IN_Beetle.10000_pi

# Step 2: Calculate Observed Heterozygosity per status
vcftools --vcf N_Beetle.filtered.variants.vcf --het --out N_Beetle.het
vcftools --vcf IN_Beetle.filtered.variants.vcf --het --out IN_Beetle.het
#############################
cd ..
