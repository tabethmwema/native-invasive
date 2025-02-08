#!/bin/bash

module load mummer 2>/dev/null || echo "Warning: Could not load MUMmer. Make sure it is installed."

BASE_DIR="/home/tzm0087/Genomes"
OUTPUT_DIR="/home/tzm0087/Genomes/Alignment_Results"

mkdir -p "$OUTPUT_DIR"

species_list=("Ants" "Birds" "Fish" "Frogs" "Mosquitoes" "Prawn" "Ticks")

for species in "${species_list[@]}"; do
    echo "Submitting MUMmer for $species..."

REFERENCE="${BASE_DIR}/${species}/Native/native_genome.fasta"
QUERY="${BASE_DIR}/${species}/Invasive/invasive_genome.fasta"
OUT_PREFIX="${OUTPUT_DIR}/${species}_alignment"

# Check if output files already exist (skip if they do)
if [[ -f "${OUT_PREFIX}.delta" && -f "${OUT_PREFIX}.coords" && -f "${OUT_PREFIX}.snps" ]]; then
    echo "Skipping $species: Alignment already completed."
    continue
fi

sbatch  <<- EOF
#!/bin/bash

#SBATCH --job-name=Mummer_$species                     #Job name
#SBATCH --output=logs/mummer_alignment_$species.log
#SBATCH --error=logs/mummer_alignment_$species.err
#SBATCH --time=5-00                                    #Runtime in D-HH:MM
#SBATCH --partition=general,jrw0107_std,nova           # Change if needed
#SBATCH --ntasks=5                                     # Number of CPU cores
#SBATCH --mem=500G                                     # Memory as in GB
#SBATCH --mail-type=END,FAIL                           # Email user when job finishes or fails
#SBATCH --mail-user=tzm0087@auburn.edu                 # Email to which notifications will be sent

echo "Running MUMmer for $species..."
echo "Alignment started at $(date)"

echo ref ${REFERENCE}


if [[ -f "$REFERENCE" && -f "$QUERY" ]]; then
    echo "Found genomes for $species, starting alignment..."
    nucmer --threads 10 --prefix="$OUT_PREFIX" "$REFERENCE" "$QUERY"

    if [[ $? -eq 0 ]]; then
        echo "Processing alignment results for $species..."
        delta-filter -r -q "${OUT_PREFIX}.delta" > "${OUT_PREFIX}_filtered.delta"
        show-coords -rcl "${OUT_PREFIX}_filtered.delta" > "${OUT_PREFIX}.coords"
        show-snps -Clr "${OUT_PREFIX}_filtered.delta" > "${OUT_PREFIX}.snps"
        echo "Alignment completed for $species!"
    else
        echo "Error: NUCmer alignment failed for $species!"
    fi
    echo "Alignment process finished at $(date)"
else
    echo "Missing genome files for $species. Skipping alignment."
fi
EOF
done

echo "All MUMmer alignment jobs submitted!"
