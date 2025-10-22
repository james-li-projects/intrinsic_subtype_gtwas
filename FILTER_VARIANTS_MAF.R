library(data.table)
library(dplyr)
library(stringr)

# importing sumstats name
args <- commandArgs(trailingOnly = TRUE)
dataset=args[1]
sumstats_name=args[2]
print(paste("PROCESSING SUMSTATS FOR:",sumstats_name))

# extracting current_subtype name
current_subtype=unlist(str_split(string=sumstats_name,pattern="\\."))[1]

# importing allele frequencies for the given current_subtype
afreq <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/afreq/",current_subtype,".afreq")) %>% mutate(MAF=ifelse(ALT_FREQS>0.5,1-ALT_FREQS,ALT_FREQS)) %>% select(ID,MAF)

# assembling df with current_subtype sample sizes
sample_size_df <- data.frame(
  subtype=c("HRNEG_HER2NEG","HRNEG_HER2POS","HRPOS_HER2NEG","HRPOS_HER2POS"),
  N=c(2760,795,5882,1296)
)

# identifying list of variants for the given current_subtype to allow results for
allowed_variant_MAF_threshold <- 30/(sample_size_df %>% filter(subtype==current_subtype))$N
print(paste("MAF Threshold based on MAC of 30:",round(allowed_variant_MAF_threshold,4)))
allowed_variant_list <- (afreq %>% filter(MAF>=allowed_variant_MAF_threshold))$ID

# importing TOP sumstats, filtered by whether they qualify to be included in analysis based on MAF
tmp_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/unfiltered_sumstats_",dataset,"/",sumstats_name)) %>% filter(ID%in%allowed_variant_list)

# writing out summary statistics
write.table(tmp_sumstats,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_",dataset,"/",sumstats_name),quote=F,row.names=F,col.names=T,sep="\t")
