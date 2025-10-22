# setting TMPDIR
export TMPDIR=/scratch/jll1/tmp

# loading modules
module load gcc/12.1.0
module load plink/2.0
module load bcftools/1.20

# chromosome index
# i=22

# list of paths
pathA=/gpfs/data/huo-lab/UKbiobank/genotype
pathB=/gpfs/data/huo-lab/UKbiobank/BreastCa
pathC=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/xancestry_ld_ref/eur

# restricting phenotype file to only 1 phenotype
cat $pathB/breastCa_pheno.txt | cut -f1-3 > $pathC/breastCa_pheno.txt

# converting bgen file to bed file
for i in `seq 1 22`
do
  echo "Converting bgen to bed for Chromosome:" $i
  plink2 --sample $pathA/ukb49564_imp_chr17_v3_s487297.sample --bgen $pathA/ukb_imp_chr${i}_v3.bgen ref-first --extract $pathB/SNPlist_ukb_mfi_v3_chr${i}.txt.gz --pheno $pathC/breastCa_pheno.txt --make-bed --out $pathC/raw_chr${i}	
done

# come back to combine bed files

# converting bfile to vcf
plink2 -bfile $pathC/raw_chr${i} --keep $pathC/keep_sample.list --recode vcf --output-chr MT --out $pathC/raw_chr${i}

# renaming chromosomes
rename_chr_file=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/xancestry_ld_ref/chr_dict/chr_name_conv.txt
input_vcf_file=$pathC/raw_chr${i}.vcf
bcftools annotate --rename-chrs ${rename_chr_file} ${input_vcf_file} -Oz -o $pathC/rename_chr${i}.vcf.gz

# running CrossMap
cd $pathC
fasta_file=/gpfs/data/pierce-lab/james.li/hg38/hg38.fa
input_chain=/gpfs/data/pierce-lab/james.li/liftOver/hg19ToHg38.over.chain.gz
CrossMap vcf ${input_chain} \
$pathC/rename_chr${i}.vcf.gz \
${fasta_file} \
$pathC/CrossMap_chr${i}.vcf

# converting CrossMap output to a bim file
plink2 --vcf $pathC/CrossMap_chr${i}.vcf --allow-extra-chr --chr $i just-acgt --sort-vars --make-bed --out $pathC/CrossMap_chr${i}






####################################################
# renaming variants for this particular chromosome #
####################################################
R
library(data.table)
library(dplyr)

chr_bim <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/xancestry_ld_ref/eur/CrossMap_chr",i))




# converting to plink bfile
plink2 --vcf /scratch/jll1/imQTL/genetic_data/CrossMap_Bangladesh_HRC_hg37_imputed_r2gtp3_indels_mafgtp005_woAFdiff_wosexchr.vcf --allow-extra-chr --chr 1-22 --snps-only just-acgt --sort-vars --make-pgen --out /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw
plink2 -pfile /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw --chr 1-22 --make-bed --out /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw

#####################################################
# creating variant ID conversion dictionaries based A1 and A2 alleles designated by plink2
conversion_bangladesh <- bangladesh_bim %>% mutate(ID=paste(V1,V4,V6,V5,sep=":")) %>% mutate(newID=paste0("chr",ID)) %>% select(V2,newID)
gtex_bim <- fread("/scratch/jll1/imQTL/genetic_data/processed_gtex.bim")
conversion_gtex <- gtex_bim %>% mutate(ID=paste(V1,V4,V6,V5,sep=":")) %>% mutate(newID=paste0("chr",ID)) %>% select(V2,newID)

# writing out only Bangladesh dictionary since GTEx IDs are all identical between old and new IDs 
write.table(conversion_bangladesh,file="/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/HEALS_ID_CONVERT_DICT_A1A2.txt",quote=F,row.names=F,col.names=F,sep="\t")

# updating variant names for HEALS genotyping data
system("module load plink/2.0; plink2 -bfile /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw --update-name /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/HEALS_ID_CONVERT_DICT_A1A2.txt --make-bed --out /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw_updatename")
# filtering this file for genotype missingness and HWE
system("module load plink/2.0; plink2 -bfile /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw_updatename --geno 0.05 --hwe 1e-6 --make-bed --out /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw_updatename_geno0.05_hwe1e6")
# further filtering to exclude 4 duplicated variants (total 8 variants removed)
system("module load plink/2.0; plink2 -bfile /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw_updatename_geno0.05_hwe1e6 --rm-dup exclude-mismatch --make-bed --out /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data")
# adding chr prefix and filtering MAF 0.01
system("module load plink/2.0; plink2 -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data --maf 0.01 --output-chr chrM --keep-allele-order --make-pgen --out /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix")