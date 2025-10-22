#!/bin/bash
#SBATCH --job-name=convert_pgen_traw
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/convert_pgen_traw_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/convert_pgen_traw_%A.err
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load plink/2.0

i=${ARGS1}

# identifying variants that have been pruned out 
plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS --chr ${i} --export Av --memory 100000 --threads 1 --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/genetic_data/GWAS_chr${i}
