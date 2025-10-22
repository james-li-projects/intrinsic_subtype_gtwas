library(dplyr)
library(data.table)

###############################################
# REMOVING PROBLEMATIC VARIANTS FROM SUMSTATS #
###############################################
# identifying problematic HRPOS_HER2NEG variants
subtype="HRPOS_HER2NEG"
ancestry="afr"
tmp_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_",ancestry,"_meta_sumstats/",subtype,".tsv")) 
problem_variant_HRPOS_HER2NEG <- (tmp_sumstats %>% filter(CHR==3,POS>197512284-500000,POS<197512284+500000) %>% filter(P<1e-4) %>% arrange(P))$ID

# identifying problematic HRNEG_HER2NEG variants
subtype="HRNEG_HER2NEG"
ancestry="afr"
tmp_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_",ancestry,"_meta_sumstats/",subtype,".tsv")) 
problem_variant_HRNEG_HER2NEG <- (tmp_sumstats %>% filter(CHR==13,POS>110604015-500000,POS<110604015+500000) %>% filter(P<1e-4) %>% arrange(P))$ID

# collating problematic variants
all_problem_variant <- data.frame(unique(c(problem_variant_HRPOS_HER2NEG,problem_variant_HRNEG_HER2NEG)))
write.table(all_problem_variant,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/problematic_variant_list/problematic_variant_list.txt",quote=F,row.names=F,col.names=F)

# filtering out variants
for (ancestry in c("afr","xancestry")) {
  for (subtype in c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRNEG_HER2NEG")) {
    print(paste("Removing problematic variants from:",ancestry,subtype))
    system(paste0("grep -v -f /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/problematic_variant_list/problematic_variant_list.txt /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_",ancestry,"_meta_sumstats/",subtype,".tsv > /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/",ancestry,"/",subtype,".tsv"))
  }
}

# copying european variants just for completeness
for (ancestry in c("eur")) {
  for (subtype in c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRNEG_HER2NEG")) {
    print(paste("Copying sumstats:",ancestry,subtype))
    system(paste0("cp /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_",ancestry,"_meta_sumstats/",subtype,".tsv /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/",ancestry,"/",subtype,".tsv"))
  }
}
