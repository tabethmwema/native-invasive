#!/bin/bash

module load mummer 2>/dev/null || echo "Warning: Could not load MUMmer. Make sure it is installed."

BASE_DIR="/home/tzm0087/Genomes"
OUTPUT_DIR="/home/tzm0087/Genomes/Alignment_Results2"

mkdir -p "$OUTPUT_DIR"

species_list=("Ticks" "Crab" "Mosquitoes" "Fish" "Birds" "Ants" "Fungus" "Mites" "Pig" "Snail" "Worms" "Beetle")

for species in "${species_list[@]}"; do
    echo "Submitting MUMmer for $species ..."

    REFERENCE="${BASE_DIR}/${species}/Reference/Ref_genome.fna"

    for category in native invasive; do

    echo "Submitting MUMmer for $species $category ..."

    QUERY="${BASE_DIR}/${species}/${category}/${category}_genome.fna"
    OUT_PREFIX="${OUTPUT_DIR}/${species}_${category}_alignment"

    # Check if output files already exist (skip if they do)
    if [[ -f "${OUT_PREFIX}.delta" && -f "${OUT_PREFIX}_filtered.delta" && -f "${OUT_PREFIX}.coords" && -f "${OUT_PREFIX}.snps" ]]; then
        echo "Skipping $species: Alignment already completed."
        continue
    fi

sbatch  <<- EOF
#!/bin/bash

#SBATCH --job-name=Mummer_${species}_$category                    #Job name
#SBATCH --output=logs/mummer_alignment_${species}_$category.log
#SBATCH --error=logs/mummer_alignment_${species}_$category.err
#SBATCH --time=10-00                                            #Runtime in D-HH:MM
#SBATCH --partition=jrw0107_std,general,nova                    # Change if needed
#SBATCH --ntasks=20                                             # Number of CPU cores
#SBATCH --mem=100G                                              # Memory as in GB
#SBATCH --mail-type=END,FAIL                                    # Email user when job finishes or fails
#SBATCH --mail-user=tzm0087@auburn.edu                          # Email to which notifications will be sent

echo "Running MUMmer for $species $category ..."
echo "Alignment started at $(date)"

echo ref ${REFERENCE}

nucmer --threads 20 --prefix="${OUT_PREFIX}" "$REFERENCE" "$QUERY"
    
if [[ $? -eq 0 ]]; then
     echo "Processing alignment results for $species $category ..."
	    delta-filter -r -q "${OUT_PREFIX}.delta" > "${OUT_PREFIX}_filtered.delta"
        show-coords -rcl "${OUT_PREFIX}_filtered.delta" > "${OUT_PREFIX}.coords"
        show-snps -T "${OUT_PREFIX}_filtered.delta" > "${OUT_PREFIX}.snps"

	echo "Alignment completed for $species!"
else
        echo "Error: NUCmer alignment failed for $species!"
    fi
    echo "Alignment process finished at $(date)"
EOF
    done
done
echo "All MUMmer alignment jobs submitted!"

