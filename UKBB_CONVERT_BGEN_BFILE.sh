#!/bin/bash
#SBATCH --job-name=UKBB_CONVERT_BGEN_BFILE
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/UKBB_CONVERT_BGEN_BFILE_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/logs/UKBB_CONVERT_BGEN_BFILE_%A.err
#SBATCH --partition=tier2q

# importing arguments 
i=${ARGS1}

# setting TMPDIR
export TMPDIR=/scratch/jll1/tmp

# loading modules
module load gcc/12.1.0
module load plink/2.0
module load bcftools/1.20
module load miniconda3/23.1.0
source activate r_env

# list of paths
pathA=/gpfs/data/huo-lab/UKbiobank/genotype
pathB=/gpfs/data/huo-lab/UKbiobank/BreastCa
pathC=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/xancestry_ld_ref/eur

# converting bgen file to bfile
plink2 --sample $pathA/ukb49564_imp_chr17_v3_s487297.sample --bgen $pathA/ukb_imp_chr${i}_v3.bgen ref-first --extract $pathB/SNPlist_ukb_mfi_v3_chr${i}.txt.gz --pheno $pathC/breastCa_pheno.txt --memory 100000 --make-bed --out $pathC/raw_chr${i}

# setting variant IDs
plink2 -bfile $pathC/raw_chr${i} --set-all-var-ids @:#:\$r:\$a --make-bed --out $pathC/set_missing_id_chr${i}

# removing duplicated variant IDs just in case
plink2 -bfile $pathC/set_missing_id_chr${i} --rm-dup --make-bed --out $pathC/set_missing_id_chr_rm_dup${i}
