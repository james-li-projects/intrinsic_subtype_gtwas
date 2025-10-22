#!/bin/bash
#SBATCH --job-name=COLLATING_TOP_PLINK2_SUMSTATS
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/COLLATING_TOP_PLINK2_SUMSTATS_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/COLLATING_TOP_PLINK2_SUMSTATS_%A.err
#SBATCH --partition=tier2q

# importing arguments 
input_sumstats=${ARGS1}

# loading modules
module load gcc/12.1.0
module load miniconda3/23.1.0
source activate r_env

# running Rscript 
Rscript /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/code/rscripts/COLLATING_TOP_PLINK2_SUMSTATS.R ${input_sumstats}
