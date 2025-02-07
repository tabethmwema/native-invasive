#!/bin/bash
#SBATCH --job-name=mummer_alignment                    #Job name
#SBATCH --output=logs/mummer_alignment.log
#SBATCH --error=logs/mummer_alignment.err
#SBATCH --time=5-00                                    #Runtime in D-HH:MM
#SBATCH --partition=general,jrw0107_std,nova                            # Change if needed
#SBATCH --ntasks=30 #Number of CPU cores
#SBATCH --mem=150G                                     # Memory as in GB
#SBATCH --mail-type=END,FAIL                           # Email user when job finishes or fails
#SBATCH --mail-user=tzm0087@auburn.edu                 # Email to which notifications will be sent

module load mummer 2>/dev/null || echo "Warning: Could not load MUMmer. Make sure it is installed."

BASE_DIR="$HOME/Genomes"
OUTPUT_DIR="$HOME/Genomes/Alignment_Results"

mkdir -p "$OUTPUT_DIR"
echo "Alignment started at $(date)"

species_list=("Ants" "Birds" "Fish" "Frogs" "Mosquitoes" "Prawn" "Ticks")

for species in "${species_list[@]}"; do
    echo "Running MUMmer for $species..."
    REFERENCE="$BASE_DIR/$species/Native/native_genome.fasta"
    QUERY="$BASE_DIR/$species/Invasive/invasive_genome.fasta"
    OUT_PREFIX="$OUTPUT_DIR/${species}_alignment"

    if [[ -f "$REFERENCE" && -f "$QUERY" ]]; then
        echo "Found genomes for $species, starting alignment..."
        nucmer --prefix="$OUT_PREFIX" "$REFERENCE" "$QUERY"

        if [[ $? -eq 0 ]]; then
            echo "Processing alignment results for $species..."
            delta-filter -r -q "${OUT_PREFIX}.delta" > "${OUT_PREFIX}_filtered.delta"
            show-coords -rcl "${OUT_PREFIX}_filtered.delta" > "${OUT_PREFIX}.coords"
            show-snps -Clr "${OUT_PREFIX}_filtered.delta" > "${OUT_PREFIX}.snps"
            echo "Alignment completed for $species!"
        else
            echo "Error: NUCmer alignment failed for $species!"
        fi
    else
        echo "Missing genome files for $species. Skipping alignment."
    fi
done

echo "Alignment process finished at $(date)"

