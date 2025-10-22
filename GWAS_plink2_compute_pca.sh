#!/bin/bash
#SBATCH --job-name=GWAS_plink2_compute_pca
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=256gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/GWAS_plink2_compute_pca.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/GWAS_plink2_compute_pca.err
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load plink/2.0

# identifying variants that have been pruned out 
plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS \
    --indep-pairwise 500kb 0.2 \
    --threads 1 \
    --memory 256000 \
    --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/PCA/pruned

# obtaining a pfile that has been pruned
plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS \
    --extract /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/PCA/pruned.prune.in \
    --threads 1 \
    --memory 256000 \
    --make-pgen \
    --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/PCA/PRUNED
    

# generating PCs
plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/PCA/PRUNED \
    --pca 50 \
    --threads 1 \
    --memory 256000 \
    --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/PCA/AABCG_PCA_50PCs
