# initializing packages
library(data.table)
library(dplyr)
library(tidyr)

# initializing subtype
for (subtype in c(
  "HRPOS_HER2NEG",
  "HRNEG_HER2NEG",
  "HRPOS_HER2POS",
  "HRNEG_HER2POS"
)) {
  
  print(paste("Processing lead variants for:",subtype))
  
  # importing summary statistic file 
  afr_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_sumstats/afr/",subtype,".tsv"))
  eur_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_sumstats/eur/",subtype,".tsv"))
  
  # importing lists of lead variants
  afr_file <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry/", subtype, "_afr.tsv")
  eur_file <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry/", subtype, "_eur.tsv")
  xancestry_file <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry/", subtype, "_xancestry.tsv")
  
  afr_lead_variants <- if (file.exists(afr_file)) (fread(afr_file))$ID else c()
  eur_lead_variants <- if (file.exists(eur_file)) (fread(eur_file))$ID else c()
  xancestry_lead_variants <- if (file.exists(xancestry_file)) (fread(xancestry_file))$ID else c()
  
  # collating these lists into ancestry-specific twas lead variants
  twas_afr_lead_variants <- unique(c(afr_lead_variants,xancestry_lead_variants))
  twas_eur_lead_variants <- unique(c(eur_lead_variants,xancestry_lead_variants))
  
  # obtaining lead variant DFs with association estimates for each ancestry for COJO TWAS
  twas_afr_lead_variant_df <- afr_sumstats %>% filter(ID%in%twas_afr_lead_variants)
  twas_eur_lead_variant_df <- eur_sumstats %>% filter(ID%in%twas_eur_lead_variants)
  
  # writing out these lead variant DFs
  write.table(twas_afr_lead_variant_df,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_lead_variants/afr/",subtype,".tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  write.table(twas_eur_lead_variant_df,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_lead_variants/eur/",subtype,".tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  
}

