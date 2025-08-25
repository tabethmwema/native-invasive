#!/bin/bash

#SBATCH --job-name=ROH_species_strict
#SBATCH --output=logs/roh_strict_%j.log
#SBATCH --error=logs/roh_strict_%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=jrw0107_std,general,nova
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tzm0087@auburn.edu

# Load PLINK
module load plink/1.9

# Set directories
VCF_DIR="VCFs"
OUT_DIR="ROH_output_snp30"
mkdir -p ${OUT_DIR} logs

# Species to process (must match the prefixes in .vcf filenames)
species_list=(
  "Mosq" "Ticks" "Fish" "Fungus" "Mite"
  "Crab" "Pig" "Beetle" "Rainbow_fish"
  "Old_World_Bollworm" "Silver_carp_fish"
)

# Loop through Native (N) and Invasive (IN)
for status in IN N; do
  for species in "${species_list[@]}"; do
    vcf_file="${VCF_DIR}/${status}_${species}.filtered.variants.vcf"
    out_prefix="${OUT_DIR}/${status}_${species}"

    if [[ -f "$vcf_file" ]]; then
      echo ">>> Starting ROH for ${status}_${species}"

      # Step 1: Convert to PLINK binary format
      plink \
        --vcf "$vcf_file" \
        --allow-extra-chr \
        --make-bed \
        --out "$out_prefix"

      # Step 2: Run ROH detection
      plink \
        --bfile "$out_prefix" \
        --homozyg \
        --homozyg-density 50 \
        --homozyg-gap 1000 \
        --homozyg-kb 100 \
        --homozyg-snp 100 \
        --homozyg-window-het 1 \
        --homozyg-window-missing 5 \
        --homozyg-window-snp 50 \
        --homozyg-window-threshold 0.05 \
        --allow-extra-chr \
        --out "${out_prefix}_roh"

      echo "✔ Done with ${status}_${species}"
    else
      echo "⚠ VCF not found or corrupt: ${vcf_file}"
    fi
  done
done

