# initializing packages
library(data.table)
library(dplyr)
library(tidyr)

# importing EAF and N annotations for each ancestry's GWAS meta-analysis studies
afr_eaf_n<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/afr_eaf_n.txt") %>% select(-N)
eur_eaf_n<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/eur_eaf_n.txt") %>% select(-N)

# initializing subtype
for (subtype in c(
  "HRPOS_HER2NEG",
  "HRNEG_HER2NEG",
  "HRPOS_HER2POS",
  "HRNEG_HER2POS"
)) {
  
  print(paste("Processing summary statistics for:",subtype))
  
  # importing summary statistics
  eur_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_meta_sumstats/",subtype,".tsv")) 
  print(nrow(eur_sumstats))
  eur_sumstats <- eur_sumstats %>% inner_join(eur_eaf_n,by=c("ID"))
  print(nrow(eur_sumstats))
  
  afr_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats/",subtype,".tsv")) 
  print(nrow(afr_sumstats))
  afr_sumstats <- afr_sumstats %>% inner_join(afr_eaf_n,by=c("ID")) %>% select(-Direction)
  print(nrow(afr_sumstats))
  
  # computing effective sample sizes
  eur_sumstats <- eur_sumstats %>% 
    mutate(N=(-((1/91477)-2*SE^2*EAF*(1-EAF))^-1) + 91477) %>%
    mutate(N_case=(-((1/91477)-2*SE^2*EAF*(1-EAF))^-1)) %>%
    mutate(N_control=91477)
  afr_sumstats <- afr_sumstats %>% 
    mutate(N=(-((1/18800)-2*SE^2*EAF*(1-EAF))^-1) + 18800) %>%  
    mutate(N_case=(-((1/18800)-2*SE^2*EAF*(1-EAF))^-1)) %>%  
    mutate(N_control=18800)
  
  # writing out summary statistics
  write.table(eur_sumstats,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_sumstats/eur/",subtype,".tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  write.table(afr_sumstats,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_sumstats/afr/",subtype,".tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  
}

