library(data.table)
library(dplyr)
library(tidyr)

# set seed 
set.seed(2025)

# determining number of controls to sample from UKBB based on the EUR-to-AFR ratio of sample sizes utilized in generating summary statistics
#eur_n = 247173
#afr_n = 36191
#afr_control_n = 18800
#eur_control_n = round((eur_n/afr_n)*afr_control_n/10)

# setting the number of controls to subsample at 5000
eur_control_n = 5000
print(paste("Number of UKBB controls to subsample:",eur_control_n))

# importing UKBB phenotype file
pheno_UKBB <- fread("/gpfs/data/huo-lab/UKbiobank/BreastCa/breastCa_pheno.txt") %>% filter(breast==1)

# randomly sample unique rows
subsampled_df <- (pheno_UKBB[sample(nrow(pheno_UKBB), eur_control_n, replace = FALSE), ]) %>% select(FID,IID)

# write out this sample list to keep in subsequent UKBB processing 
write.table(subsampled_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/xancestry_ld_ref/eur/keep_sample.list",quote=F,row.names=F,col.names=F)
