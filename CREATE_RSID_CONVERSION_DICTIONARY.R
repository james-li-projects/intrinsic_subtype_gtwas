library(data.table)
library(dplyr)

# importing entire rsid dictionary
info<-readRDS("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/rsid_dictionary/SNP_GRCh37_38_match_update.rds")

# assembling a more direct rsid dictionary to convert with chr:pos:ref:alt
variant_dict <- info %>% mutate(ID=paste(chr,pos38,allele2_38,allele1_38,sep=":")) %>% select(ID,rsid)

# obtaining all variant IDs in this study
all_ID <- data.frame()
for (subtype in c(
  "HRPOS_HER2NEG",
  "HRNEG_HER2NEG",
  "HRPOS_HER2POS",
  "HRNEG_HER2POS"
)) {
  raw_eur_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_meta_sumstats/",subtype,".tsv")) %>% select(ID)
  raw_afr_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats/",subtype,".tsv")) %>% select(ID)
  raw_xancestry_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_xancestry_meta_sumstats/",subtype,".tsv")) %>% select(ID)
  all_ID=rbind(
    all_ID,
    raw_eur_sumstats,
    raw_afr_sumstats,
    raw_xancestry_sumstats
  ) %>% unique()
}

# generating final rsid conversion dictionary
final_variant_dict <- variant_dict %>% filter(ID %in% all_ID$ID)

# writing out this ID to RSID conversion dictionary
write.table(final_variant_dict,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/rsid_dictionary/ID_RSID_dict.tsv",quote=F,row.names=F,col.names=T,sep="\t")
