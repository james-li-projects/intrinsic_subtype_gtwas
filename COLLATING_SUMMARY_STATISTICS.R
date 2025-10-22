library(data.table)
library(dplyr)


MEGA<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/top_combined_sumstats/MEGA.sumstats") %>% filter(variant_id!="variant_id")

HRNEG_HER2NEG <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/plink2_gwas/HRNEG_HER2NEG.Status.glm.logistic.hybrid") 
HRNEG_HER2NEG <- HRNEG_HER2NEG %>% mutate(BETA=log(OR))

identical(MEGA$variant_id,HRNEG_HER2NEG$ID)


# summary(lm(as.numeric(MEGA$HRNEG_HER2NEG_BETA)~HRNEG_HER2NEG$BETA))


hi<-data.frame(MEGA$variant_id,as.numeric(MEGA$HRNEG_HER2NEG_BETA),as.numeric(MEGA$HRNEG_HER2NEG_P),HRNEG_HER2NEG$BETA,HRNEG_HER2NEG$P)
colnames(hi) <- c("variant_id","TOP_BETA","TOP_P","PLINK2_BETA","PLINK2_P")
hi_sig <- hi %>% filter(TOP_P<5e-8)
hi_sig

summary(lm(hi_sig$TOP_BETA~hi_sig$PLINK2_BETA))

#hi_sig <- hi %>% filter(PLINK2_P<5e-8)
#hi_sig
