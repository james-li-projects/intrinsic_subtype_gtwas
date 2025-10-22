###########################
# importing libraries
library(readxl)
library(dplyr)
library(TOP)
library(data.table)
###########################

###########################
# reading in lead variants specific to AFR
lead_variant_path="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/"
HRNEG_HER2NEG <- fread(paste0(lead_variant_path,"/HRNEG_HER2NEG.txt"))
HRPOS_HER2NEG <- fread(paste0(lead_variant_path,"/HRPOS_HER2NEG.txt")) 
HRPOS_HER2POS <- fread(paste0(lead_variant_path,"/HRPOS_HER2POS.txt")) 
lead_variant_list <- rbind(HRNEG_HER2NEG,HRPOS_HER2NEG,HRPOS_HER2POS) %>% unique()
write.table(lead_variant_list,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/by_study_AF/variant.list",quote=F,row.names=F,col.names=F)
###########################

###########################
# importing covariate data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
# converting breast cancer status to a binary variable
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,`Dataset...3`,Age_GWAS,Status,ER,PR,HER2,AFR_pro)
pheno_data$Status <- 2-pheno_data$Status
colnames(pheno_data)[1] <- "Sample_Name"
colnames(pheno_data)[2] <- "Dataset"

# importing data of principal components
eigenvec <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/PCA/AABCG_PCA_50PCs.eigenvec",header=T)
colnames(eigenvec)[1] <- "Sample_Name"

# joining covariate and PC data.frames by sample IID
combined_cov <- inner_join(pheno_data,eigenvec, by = c("Sample_Name")) 
# recoding the dataset variable 
combined_cov$dataset_recode <- NA
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_01_WGS"] <- "WGS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_02_WGS2"] <- "WGS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_03_MEGA_VANDY"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_04_MEGA_RP"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_05_MEGA_USC"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_06_AMBER"] <- "AMBER" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_07_ROOT"] <- "ROOT" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_08_AABC"] <- "AABC" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_09_GBHS"] <- "GBHS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_10_BCAC_OncoArray"] <- "BCAC_OncoArray" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_11_BCAC_iCOGS"] <- "BCAC_iCOGS" 

# write out sample lists 
platform_list = unique(combined_cov$dataset_recode)

for (current_platform in platform_list) {
  print(paste("Computing allele frequencies for:",current_platform))
  write.table(combined_cov %>% filter(dataset_recode==current_platform) %>% filter(Status==0) %>% select(Sample_Name),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/by_study_AF/",current_platform,".list"),quote=F,row.names=F,col.names=F,sep="\t")
  system(paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS --extract /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/by_study_AF/variant.list --keep /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/by_study_AF/",current_platform,".list --freq --memory 100000 --threads 1 --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/by_study_AF/",current_platform))
}

library(dplyr)
library(readr)
library(purrr)
library(stringr)

# Set the directory path
afreq_dir <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/by_study_AF/"

# Get list of all *.afreq files
afreq_files <- list.files(path = afreq_dir, pattern = "\\.afreq$", full.names = TRUE)

# Function to read and process each afreq file
read_and_rename <- function(file) {
  # Extract name prefix before .afreq
  name_prefix <- str_remove(basename(file), "\\.afreq$")
  
  # Read in the file
  df <- read_tsv(file, show_col_types = FALSE)
  
  # Rename ALT_FREQS column
  df <- df %>%
    select(-OBS_CT) %>%
    rename(!!name_prefix := ALT_FREQS)
  
  return(df)
}

# Read and process all files into a list
afreq_list <- lapply(afreq_files, read_and_rename)

# Reduce list into single dataframe joined by the specified keys
afreq_merged <- reduce(afreq_list, full_join, by = c("#CHROM", "ID", "REF", "ALT"))


# write out output table
write.table(afreq_merged,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/by_study_AF/merged_AF.tsv",quote=F,row.names=F,col.names=T,sep="\t")
