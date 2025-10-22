#!/bin/bash
#SBATCH --job-name=RUN_TOP_GWAS
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/RUN_TOP_GWAS_%j.log
#SBATCH --array=1-2943%300

module load gcc/12.1.0
module load miniconda3
source activate r_env


# Convert the SLURM array task ID to a 5-digit zero-padded string
formatted_index_str=$(printf "%05d" $SLURM_ARRAY_TASK_ID)

# Run the R script with the formatted index
/gpfs/data/huo-lab/BCAC/james.li/conda_env/r_env/bin/Rscript /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/code/rscripts/RUN_TOP_GWAS.R partition_${formatted_index_str}.txt
