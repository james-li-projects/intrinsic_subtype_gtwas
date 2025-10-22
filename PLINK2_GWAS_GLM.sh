#!/bin/bash
#SBATCH --job-name=PLINK2_GWAS_GLM
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PLINK2_GWAS_GLM_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PLINK2_GWAS_GLM_%A.err
#SBATCH --partition=tier2q

# importing arguments 
subtype=${ARGS1}
dataset=${ARGS2}

# generating a string for 10 PCs
PC_index_list=`seq 1 10`
PC_list=`for i in $PC_index_list; do echo PC${i}; done`

# loading modules
module load gcc/12.1.0
module load plink/2.0

# specifying paths for input and output
in_path=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile
pheno_cov_path=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pheno_cov
out_path=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/plink2_gwas

# running GWAS using plink2 [logistic regression] 
plink2 --pfile ${in_path}/GWAS \
    --glm hide-covar \
    --pheno ${pheno_cov_path}/${subtype}.${dataset}.pheno_cov \
    --pheno-name Status \
    --covar ${pheno_cov_path}/${subtype}.${dataset}.pheno_cov \
    --covar-name Age ${PC_list} \
    --covar-variance-standardize \
    --threads 1 \
    --memory 100000 \
    --out ${out_path}/${subtype}.${dataset}
