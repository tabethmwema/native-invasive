# Extracted analysis commands for: Asian tiger mosquito (Mosq)
# Source: Raw_script_script_updated_V2
# Extracted on: 2025-08-19S

mkdir -p Genome_data/Mosq
cd Genome_data/Mosq
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/444/175/GCA_001444175.2_A.albopictus_v1.1/GCA_001444175.2_A.albopictus_v1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/035/046/485/GCF_035046485.1_AalbF5/GCF_035046485.1_AalbF5_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/516/635/GCA_006516635.1_Aalbo_alts_of_PRJNA530512/GCA_006516635.1_Aalbo_alts_of_PRJNA530512_genomic.fna.gz

##### running Native
minimap2 -ax asm5 -t 20 Mosquitoes_GCF_035046485.1_ref.fna N_Mosquitoes_GCA_001444175.2.fna > N_Mosquitoes.sam
samtools view -@ 20 -bS N_Mosquitoes.sam | samtools sort -@ 20 -o N_Mosquitoes.sorted.bam
samtools index N_Mosquitoes.sorted.bam
bcftools mpileup -Ou -f Mosquitoes_GCF_035046485.1_ref.fna N_Mosquitoes.sorted.bam | bcftools call -mv -Ov -o N_Mosquitoes_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' N_Mosquitoes_variants.vcf > N_Mosquitoes_filtered_variants.vcf

###### running Invasiv
minimap2 -ax asm5 -t 20 Mosquitoes_GCF_035046485.1_ref.fna In_Mosquitoes_GCA_006516635.1.fna > IN_Mosquitoes.sam
samtools view -bS IN_Mosquitoes.sam | samtools sort -o IN_Mosquitoes.sorted.bam
samtools index IN_Mosquitoes.sorted.bam
bcftools mpileup -Ou -f Mosquitoes_GCF_035046485.1_ref.fna IN_Mosquitoes.sorted.bam | bcftools call -mv -Ov -o IN_Mosquitoes_variants.vcf
bcftools filter -s LOWQUAL -e '%QUAL<20' IN_Mosquitoes_variants.vcf > IN_Mosquitoes_filtered_variants.vcf

# To calculate nucleotide diversity (π), heterozygosity
bgzip IN_Mosquitoes_filtered_variants.vcf
bgzip N_Mosquitoes_filtered_variants.vcf
tabix -p vcf IN_Mosquitoes_filtered_variants.vcf.gz
tabix -p vcf N_Mosquitoes_filtered_variants.vcf.gz

# Step 1: Calculate Nucleotide Diversity (π) per status
vcftools --vcf N_Mosquitoes_filtered_variants.vcf.gz --window-pi 10000 --out N_Mosquitoes.10000_pi
vcftools --vcf IN_Mosquitoes_filtered_variants.vcf.gz --window-pi 10000 --out IN_Mosquitoes.10000_pi

# Step 2: Calculate Observed Heterozygosity per status
vcftools --vcf N_Mosquitoes_filtered_variants.vcf.gz --het --out N_Mosquitoes.het
vcftools --vcf IN_Mosquitoes_filtered_variants.vcf.gz --het --out IN_Mosquitoes.het

#############################
cd ..
