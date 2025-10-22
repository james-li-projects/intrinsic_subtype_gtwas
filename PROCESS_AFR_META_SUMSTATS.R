library(data.table)
library(dplyr)
library(tidyr)

# Define the subtypes to iterate over
subtypes <- c("HRPOS_HER2NEG", "HRPOS_HER2POS", "HRNEG_HER2POS", "HRNEG_HER2NEG")

# parsing sumstats for each subtype
for (subtype in subtypes) {
  print(paste("Arranging sumstats for:",subtype)) 
  
  # importing AFR sumstats and making alleles consistent with IDs
  tmp_afr_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/COLLATED/METAANALYSIS_",subtype,"_1.tbl"))
  tmp_afr_sumstats <- tmp_afr_sumstats %>% separate(MarkerName,into=c("chr","pos","a2","a1"),remove=F,sep="\\:")
  tmp_afr_sumstats <- tmp_afr_sumstats %>% mutate(
    chr=as.numeric(chr),
    pos=as.numeric(pos)
  )
  
  # making uppercase alleles
  parsed_tmp_afr_sumstats <- tmp_afr_sumstats %>% 
    mutate(Allele1=toupper(Allele1),
           Allele2=toupper(Allele2))
  
  # altering direction and effect columns if A1 and Allele1 do not match
  parsed_tmp_afr_sumstats <- parsed_tmp_afr_sumstats %>%
    mutate(
      Direction = if_else(Allele1 != a1, chartr('+-', '-+', Direction), Direction),
      BETA=ifelse(Allele1==a1,Effect,-Effect)
    ) %>% mutate(
      EffectAllele=a1,
      BaselineAllele=a2
    ) %>% mutate(ID=MarkerName,SE=StdErr,P=`P-value`) %>% select(ID,chr,pos,EffectAllele,BaselineAllele,BETA,SE,P,Direction) %>% rename(CHR=chr,POS=pos)
  
  # filtering out any variants with NA values for chromosomes after converting to numeric 
  parsed_tmp_afr_sumstats <- parsed_tmp_afr_sumstats %>% filter(!is.na(CHR)) %>% filter(!is.na(POS)) %>% mutate(CHR=as.integer(CHR),POS=as.integer(POS))
  
  # writing out afr meta sumstats 
  write.table(parsed_tmp_afr_sumstats,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats/",subtype,".tsv"),quote=F,row.names=F,col.names=T,sep="\t")
}
