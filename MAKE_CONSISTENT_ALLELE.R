# initializing libraries
library(data.table)
library(tidyverse)
library(stringr)

###########################################
# importing arguments
args <- commandArgs(trailingOnly = TRUE)
input_directory <- args[1]
input_sumstats <- args[2]

###########################################
# parsing summary statistics depending on whether the file was generated using TOP or plink2
if (input_directory == "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/plink2_gwas") {
  #input_directory="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/plink2_gwas"
  #input_sumstats="HRPOS_HER2NEG.WGS.Status.glm.logistic.hybrid"
  
  # importing summary statistics
  print("Importing summary statistics generated using plink2")
  current_sumstats <- fread(paste0(input_directory,"/",input_sumstats))
  
  # standardizing alleles based on alternate allele
  print("Standardizing alleles")
  current_sumstats <- current_sumstats %>% separate(ID,into=c("V1","V2","BaselineAllele","EffectAllele"),remove=F,sep="\\:")
  current_sumstats <- current_sumstats %>% mutate(OR=as.numeric(OR),`LOG(OR)_SE`=as.numeric(`LOG(OR)_SE`))
  parsed_current_sumstats <- current_sumstats %>% mutate(BETA=ifelse(EffectAllele==A1,log(OR),-log(OR))) %>% rename(SE=`LOG(OR)_SE`)
  parsed_current_sumstats <- parsed_current_sumstats %>% mutate(
    BETA=ifelse(BETA==-Inf,NA,BETA),
    SE=ifelse(SE==-Inf,NA,SE),
    P=ifelse(P==-Inf,NA,P),
    BETA=ifelse(BETA==Inf,NA,BETA),
    SE=ifelse(SE==Inf,NA,SE),
    P=ifelse(P==Inf,NA,P)
  )
  
  # writing consistent allele summary statistics
  print("Writing out sumstats with consistent alleles")
  output_filename <- gsub(".Status.glm.logistic.hybrid",".sumstats",paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/unfiltered_sumstats_plink2/",input_sumstats))
  write.table(parsed_current_sumstats %>% select(ID,EffectAllele,BaselineAllele,BETA,SE,P),file=output_filename,quote=F,row.names=F,col.names=T,sep="\t")
} else if (input_directory=="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/top_combined_sumstats") {
  #input_directory="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/top_combined_sumstats"
  #input_sumstats="WGS.sumstats"
  
  # importing summary statistics
  print("Importing summary statistics generated using TOP")
  current_sumstats <- fread(paste0(input_directory,"/",input_sumstats))
  
  # standardizing alleles based on alternate allele
  print("Standardizing alleles")
  current_sumstats <- current_sumstats %>% separate(variant_id,into=c("V1","V2","V3","V4"),remove=F,sep="\\:")
  parsed_current_sumstats <- current_sumstats %>%
    mutate(across(contains("_BETA"), ~ ifelse(V4 != EffectAllele, -., .)))
  
  # making data.frames for each subtype 
  for (subtype in c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRNEG_HER2NEG")) {
    current_subtype_df <- parsed_current_sumstats %>%
      select(variant_id, EffectAllele, BaselineAllele, contains(paste0(subtype,"_BETA")), contains(paste0(subtype,"_SE")), contains(paste0(subtype,"_P")))
    colnames(current_subtype_df) <- c("ID","EffectAllele","BaselineAllele","BETA","SE","P")
    output_filename <- gsub("\\.BCAC\\.",".BCAC_OncoArray.",paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/unfiltered_sumstats_TOP/",subtype,".",input_sumstats))
    write.table(current_subtype_df %>% select(ID,EffectAllele,BaselineAllele,BETA,SE,P),file=output_filename,quote=F,row.names=F,col.names=T,sep="\t")
  }
}
