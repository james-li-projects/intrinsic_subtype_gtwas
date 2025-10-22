#!/bin/bash

# Check if the required argument (output subdirectory) is provided
if [ "$#" -ne 1 ]; then
    echo "Approach for running GWAS is missing. Please include approach as a command line argument."
    echo "Usage: $0 <approach>"
    exit 1
fi

# Define directories
approach="$1"  # Get output subdirectory from command line argument
DATA_DIR="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_${approach}"
OUTPUT_DIR="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/${approach}"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define the subtypes to analyze
SUBTYPES=("HRNEG_HER2NEG" "HRPOS_HER2NEG" "HRNEG_HER2POS" "HRPOS_HER2POS")

# Loop through each subtype and perform meta-analysis
for SUBTYPE in "${SUBTYPES[@]}"; do
    METAL_SCRIPT="${OUTPUT_DIR}/metal_analysis_${SUBTYPE}.txt"
    OUTPUT_PREFIX="${OUTPUT_DIR}/METAANALYSIS_${SUBTYPE}"

    # Find all relevant input files for this subtype
    file_list=$(ls ${DATA_DIR}/${SUBTYPE}.*.sumstats 2>/dev/null)

    # Skip if no files are found for this subtype
    if [ -z "$file_list" ]; then
        echo "No input files found for ${SUBTYPE}. Skipping..."
        continue
    fi

    # Create the METAL script for this subtype
    echo "# METAL meta-analysis script for ${SUBTYPE}" > "$METAL_SCRIPT"
    echo "SCHEME STDERR" >> "$METAL_SCRIPT"
    echo "MARKER ID" >> "$METAL_SCRIPT"
    echo "ALLELE EffectAllele BaselineAllele" >> "$METAL_SCRIPT"
    echo "EFFECT BETA" >> "$METAL_SCRIPT"
    echo "STDERR SE" >> "$METAL_SCRIPT"
    echo "PVALUE P" >> "$METAL_SCRIPT"

    # Loop through all found files and add them to the METAL script
    for file in $file_list; do
        echo "PROCESS $file" >> "$METAL_SCRIPT"
    done

    # Specify output files (store in the output directory)
    echo "OUTFILE ${OUTPUT_PREFIX}_ .tbl" >> "$METAL_SCRIPT"
    echo "ANALYZE HETEROGENEITY" >> "$METAL_SCRIPT"
    echo "QUIT" >> "$METAL_SCRIPT"

    # Run METAL with the generated script and store logs
    metal "$METAL_SCRIPT" > "${OUTPUT_DIR}/metal_log_${SUBTYPE}.txt" 2>&1

    echo "Meta-analysis for ${SUBTYPE} completed. Results stored in: $OUTPUT_DIR"
done

echo "All meta-analyses completed."
