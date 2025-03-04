#!/bin/bash

#SBATCH --job-name=Mummer2VCF                      #Job name
#SBATCH --output=logs/Mummer2VCF.log
#SBATCH --time=5-00                                   #Runtime in D-HH:MM
#SBATCH --partition=jrw0107_std,general,nova           # Change if needed
#SBATCH --ntasks=5                                     # Number of CPU cores
#SBATCH --mem=100G                                     # Memory as in GB
#SBATCH --mail-type=END,FAIL                           # Email user when job finishes $
#SBATCH --mail-user=tzm0087@auburn.edu

# Load dependencies
module load bcftools/1.17
module load samtools

source /home/tzm0087/miniforge3/bin/activate Python

# Define directories
ALIGN_DIR="/home/tzm0087/Genomes/Alignment_Results2"
VCF_DIR="/home/tzm0087/Genomes/VCFs"
REF_DIR="/home/tzm0087/Genomes"

# Ensure output directory exists
mkdir -p "$VCF_DIR"

# List of species
species_list=("Ticks" "Crab" "Drosophila" "Mosquitoes" "Fish" "Birds" "Ants" "Fungus" "Mites" "Pig" "Snail" "Worms" "Beetle")

# Run conversion for each species
for species in "${species_list[@]}"; do
    echo "Processing $species..."

    # Define input SNP files
    NATIVE_SNPS="$ALIGN_DIR/${species}_native_alignment.snps"
    INVASIVE_SNPS="$ALIGN_DIR/${species}_invasive_alignment.snps"

    # Determine reference genome
    REFERENCE_GENOME="$REF_DIR/${species}/Reference/Ref_genome.fna"

    # Define output VCF files
    NATIVE_VCF="$VCF_DIR/${species}_native.vcf"
    INVASIVE_VCF="$VCF_DIR/${species}_invasive.vcf"

    # Generate VCF for native population using the native reference genome
    echo "Generating VCF for Native $species..."
    /home/tzm0087/Genomes/all2vcf/all2vcf mummer \
        --snps "$NATIVE_SNPS" \
        --reference "$REFERENCE_GENOME" \
        --input-header --output-header > "$NATIVE_VCF"

    # Compress and index native VCF
    # gzip -c "$NATIVE_VCF" > "${NATIVE_VCF}.gz"
    # tabix -p vcf "${NATIVE_VCF}.gz"

    # Generate VCF for invasive population using the native reference genome
    echo "Generating VCF for Invasive $species..."
    /home/tzm0087/Genomes/all2vcf/all2vcf mummer \
        --snps "$INVASIVE_SNPS" \
        --reference "$REFERENCE_GENOME" \
        --input-header --output-header > "$INVASIVE_VCF"

    # Compress and index invasive VCF
    # gzip -c "$INVASIVE_VCF" > "${INVASIVE_VCF}.gz"
    # tabix -p vcf "${INVASIVE_VCF}.gz"

    echo "Completed: ${NATIVE_VCF}.gz and ${INVASIVE_VCF}.gz"
done

