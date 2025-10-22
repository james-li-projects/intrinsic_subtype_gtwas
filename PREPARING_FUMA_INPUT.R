library(dplyr)
library(data.table)

######################
# PARSING FUMA INPUT #
######################
# importing rsid conversion dictionary
ID_RSID_dict <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/rsid_dictionary/ID_RSID_dict.tsv")

# parsing each subtype
for (subtype in c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRNEG_HER2NEG")) {
  # compiling cross-ancestry sumstats for FUMA
  ancestry="xancestry"
  FUMA_INPUT <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/",ancestry,"/",subtype,".tsv")) %>% inner_join(ID_RSID_dict) %>% rename(SNP=rsid,BP=POS,A1=EffectAllele,A2=BaselineAllele) %>% mutate(OR=exp(BETA),Beta=BETA) %>% select(SNP,CHR,BP,A1,A2,P,OR,Beta,SE)
  # writing out final FUMA input
  write.table(FUMA_INPUT, file = gzfile(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/FUMA_INPUT/FUMA_INPUT_",subtype,".tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

