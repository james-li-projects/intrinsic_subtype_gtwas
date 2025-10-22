library(data.table)
library(dplyr)
library(tidyr)
library(topr)

previous_BC_gwas_hits <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/gwas_catalog/gwas-association-downloaded_2025-04-09-MONDO_0007254-withChildTraits.tsv") %>% mutate(P=as.numeric(`P-VALUE`)) %>% filter(P<5e-8) %>% select(CHR_ID,CHR_POS,SNPS) %>% rename(CHROM=CHR_ID,POS=CHR_POS,ID=SNPS) %>% mutate(CHROM=paste0("chr",CHROM),REF="A",ALT="C",BETA=2,SE=0.02,P=1e-200,class="Previous")

novel_lead_df <- data.frame()
for (subtype in c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2NEG")) {
  current_lead_variant_df <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/parsed_output/",subtype,"_FINAL_LEAD_VARIANT_DF.tsv")) %>% select(ID) %>% separate(ID,into=c("CHROM","POS","REF","ALT"),sep="\\:",remove=F) %>% mutate(CHROM=paste0("chr",CHROM),BETA=1,SE=0.02,P=4.9e-8,class="Novel")
  combined_df <- rbind(
    previous_BC_gwas_hits,
    current_lead_variant_df
  )
  tmp_novel_lead_df <- (get_lead_snps(combined_df,thresh = 5e-08,region_size = 4e+06) %>% filter(class=="Novel")) %>% select(ID) %>% mutate(SUBTYPE=subtype)
  novel_lead_df <- rbind(
    novel_lead_df,
    tmp_novel_lead_df
  )
}

# printing out completely novel variants
print(novel_lead_df)

# tmp_novel_lead_df %>% filter(POS==10503552)
