# Load required libraries
library(data.table)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)

########################################
# importing list of lead variants
HRPOS_HER2POS_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRPOS_HER2POS.txt"))$ID
HRPOS_HER2NEG_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRPOS_HER2NEG.txt"))$ID
HRNEG_HER2NEG_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRNEG_HER2NEG.txt"))$ID

########################################
# extracting heterogeneity estimates between ancestries
HRPOS_HER2POS_het <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/XANCESTRY_META/output/XANCESTRY_METAANALYSIS_HRPOS_HER2POS_1.tbl") %>% filter(MarkerName %in% HRPOS_HER2POS_lead_variant_list) %>% rename(ID=MarkerName) %>% select(ID,HetPVal) %>% arrange(ID)
HRPOS_HER2NEG_het <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/XANCESTRY_META/output/XANCESTRY_METAANALYSIS_HRPOS_HER2NEG_1.tbl") %>% filter(MarkerName %in% HRPOS_HER2NEG_lead_variant_list) %>% rename(ID=MarkerName) %>% select(ID,HetPVal) %>% arrange(ID)
HRNEG_HER2NEG_het <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/XANCESTRY_META/output/XANCESTRY_METAANALYSIS_HRNEG_HER2NEG_1.tbl") %>% filter(MarkerName %in% HRNEG_HER2NEG_lead_variant_list) %>% rename(ID=MarkerName) %>% select(ID,HetPVal) %>% arrange(ID)

########################################
# assembling DF with all lead variants and their heterogeneity estimates and writing out this DF as a tsv
heterogeneity_ancestry_df <- rbind(
  HRPOS_HER2POS_het,
  HRPOS_HER2NEG_het,
  HRNEG_HER2NEG_het
)
write.table(heterogeneity_ancestry_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/ancestry_heterogeneity/heterogeneity_ancestry.tsv",quote=F,row.names=F,col.names=T,sep="\t")
