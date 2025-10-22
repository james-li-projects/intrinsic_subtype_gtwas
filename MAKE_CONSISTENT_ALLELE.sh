#!/bin/bash
#SBATCH --job-name=MAKE_CONSISTENT_ALLELE
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/MAKE_CONSISTENT_ALLELE_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/MAKE_CONSISTENT_ALLELE_%A.err
#SBATCH --partition=tier2q

# importing arguments 
input_directory=${ARGS1}
input_sumstats=${ARGS2}

# loading modules
module load gcc/12.1.0
module load miniconda3/23.1.0
source activate r_env

# running Rscript 
Rscript /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/code/rscripts/MAKE_CONSISTENT_ALLELE.R ${input_directory} ${input_sumstats}
