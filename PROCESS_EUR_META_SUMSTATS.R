library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

# import EUR summary statistics and filter out variants with invalid alleles
EUR_SUMSTATS <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/EUR_SUMSTATS/BCAC_4subtypes_sum_data.txt") %>% 
  mutate(col4 = as.integer(seq(nrow(.)))) %>% 
  filter(Effect.Meta!="") %>% 
  filter(Baseline.Meta!="") %>%
  filter(grepl("^[ATCG]+$", Effect.Meta) & grepl("^[ATCG]+$", Baseline.Meta)) %>%
  mutate(col4 = as.integer(col4)) %>% 
  mutate(chr_hg19=chr.iCOGs,pos_hg19=Position.iCOGs) %>% 
  select(-matches("iCOGs|Onco"))

# preparing DF to liftOver
tmp_liftOver <- EUR_SUMSTATS %>% mutate(col1 = paste0(rep("chr",nrow(EUR_SUMSTATS)), chr_hg19)) %>% mutate(col2 = pos_hg19) %>% mutate(col2 = as.integer(col2)) %>% mutate(col3 = as.integer(col2 + 1)) %>% select(col1,col2,col3,col4) 

# writing temporary table out for liftOver (and operate from the specified tmp directory)
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/tmp")
library("withr")
withr::with_options(
  list(scipen = 999),
  write.table(tmp_liftOver, file="tmp_liftOver", quote=F, row.names=F, col.names=F, sep="\t")
)

# performing liftOver
system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver tmp_liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped liftOver_unmapped")

# reading in liftOver mapped results
liftOver_mapped <- fread("liftOver_mapped",header=F,sep="\t")
colnames(liftOver_mapped) <- c("col1","col2","col3","col4")

# joining with original table 
liftOver_EUR_SUMSTATS <- inner_join(EUR_SUMSTATS, liftOver_mapped, by = c("col4"))

# parsing these liftOver EUR sumstats
liftOver_EUR_SUMSTATS <- liftOver_EUR_SUMSTATS %>%
  rename(EffectAllele1=Effect.Meta,BaselineAllele1=Baseline.Meta) %>%
  mutate(
    EffectAllele2=BaselineAllele1,
    BaselineAllele2=EffectAllele1,
    chr_hg38=gsub("chr","",col1),
    pos_hg38=gsub("chr","",col2),
    ID1 = paste(chr_hg38,pos_hg38,BaselineAllele1,EffectAllele1,sep=":"),
    ID2 = paste(chr_hg38,pos_hg38,BaselineAllele2,EffectAllele2,sep=":")) 

# Rename columns based on specified substitutions
colnames(liftOver_EUR_SUMSTATS) <- gsub("Luminal_A", "HRPOS_HER2NEG", colnames(liftOver_EUR_SUMSTATS))
colnames(liftOver_EUR_SUMSTATS) <- gsub("Luminal_B", "HRPOS_HER2POS", colnames(liftOver_EUR_SUMSTATS))
colnames(liftOver_EUR_SUMSTATS) <- gsub("HER2_Enriched", "HRNEG_HER2POS", colnames(liftOver_EUR_SUMSTATS))
colnames(liftOver_EUR_SUMSTATS) <- gsub("Triple_Neg", "HRNEG_HER2NEG", colnames(liftOver_EUR_SUMSTATS))
colnames(liftOver_EUR_SUMSTATS) <- gsub("se_meta", "SE", colnames(liftOver_EUR_SUMSTATS))
colnames(liftOver_EUR_SUMSTATS) <- gsub("log_or_meta", "BETA", colnames(liftOver_EUR_SUMSTATS))

# compute p-values
liftOver_EUR_SUMSTATS <- liftOver_EUR_SUMSTATS %>% 
  mutate(
    HRPOS_HER2NEG_P=2*pnorm(-abs(HRPOS_HER2NEG_BETA/HRPOS_HER2NEG_SE), lower.tail = T),
    HRPOS_HER2POS_P=2*pnorm(-abs(HRPOS_HER2POS_BETA/HRPOS_HER2POS_SE), lower.tail = T),
    HRNEG_HER2POS_P=2*pnorm(-abs(HRNEG_HER2POS_BETA/HRNEG_HER2POS_SE), lower.tail = T),
    HRNEG_HER2NEG_P=2*pnorm(-abs(HRNEG_HER2NEG_BETA/HRNEG_HER2NEG_SE), lower.tail = T)
  )

# Define the subtypes to iterate over
subtypes <- c("HRPOS_HER2NEG", "HRPOS_HER2POS", "HRNEG_HER2POS", "HRNEG_HER2NEG")

# importing AFR GWAS pvar file to see the format of aligned variants
pvar <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS.pvar")

# parsing sumstats for each subtype
for (subtype in subtypes) {
  print(paste("Arranging sumstats for:",subtype)) 
  
  # Select and rename relevant columns
  df_subset <- liftOver_EUR_SUMSTATS %>%
    select(ID1, 
           ID2, 
           chr_hg38,
           pos_hg38,
           EffectAllele1, BaselineAllele1, 
           EffectAllele2, BaselineAllele2,
           all_of(paste0(subtype, "_BETA")),
           all_of(paste0(subtype, "_SE")),
           all_of(paste0(subtype, "_P"))) %>%
    rename(BETA = all_of(paste0(subtype, "_BETA")),
           SE = all_of(paste0(subtype, "_SE")),
           P = all_of(paste0(subtype, "_P"))) %>% 
    rename(BETA1=BETA) %>%
    mutate(BETA2=-BETA1) %>% rename(CHR=chr_hg38,POS=pos_hg38)
  
  # parsing out variants shared with AFR in orientation #1 
  df_subset_1 <- df_subset %>% select(
    ID1,
    CHR,
    POS,
    EffectAllele1,
    BaselineAllele1,
    BETA1,
    SE,
    P
  ) %>% filter(ID1 %in% pvar$ID) %>% rename(
    ID=ID1,
    EffectAllele=EffectAllele1,
    BaselineAllele=BaselineAllele1,
    BETA=BETA1,
    SE=SE,
    P=P
  )
  
  # parsing out variants shared with AFR in orientation #2
  df_subset_2 <- df_subset %>% select(
    ID1,ID2,
    CHR,
    POS,
    EffectAllele2,
    BaselineAllele2,
    BETA2,
    SE,
    P
  ) %>% filter(ID2 %in% pvar$ID) %>% filter(!(ID1%in%df_subset_1$ID)) %>% select(-ID1) %>% rename(
    ID=ID2,
    EffectAllele=EffectAllele2,
    BaselineAllele=BaselineAllele2,
    BETA=BETA2,
    SE=SE,
    P=P
  )
  
  # parsing out variants shared with AFR in orientation #2 and identifying orientation #1 IDs for these variants
  df_subset_2_ID1_list <- (df_subset %>% select(ID1,ID2) %>% filter(ID2 %in% pvar$ID))$ID1
  
  # parsing out variants unique to EUR dataset
  df_subset_unique_EUR <- df_subset %>% select(
    ID1,
    CHR,
    POS,
    EffectAllele1,
    BaselineAllele1,
    BETA1,
    SE,
    P
  ) %>% filter(!(ID1 %in% c(pvar$ID,df_subset_2_ID1_list))) %>% rename(
    ID=ID1,
    EffectAllele=EffectAllele1,
    BaselineAllele=BaselineAllele1,
    BETA=BETA1,
    SE=SE,
    P=P
  )
  
  # combining all data subsets with proper allele alignment
  df_subset_combined <- rbind(
    df_subset_1,
    df_subset_2,
    df_subset_unique_EUR
  )
  
  # removing records with duplicated IDs
  unique_ids <- df_subset_combined$ID[!df_subset_combined$ID %in% df_subset_combined$ID[duplicated(df_subset_combined$ID) | duplicated(df_subset_combined$ID, fromLast = TRUE)]]
  df_subset_combined <- df_subset_combined[df_subset_combined$ID %in% unique_ids, ]

  # printing out head just to see 
  print(head(df_subset_combined))
  
  # removing variants with NA values for chromosomes
  df_subset_combined <- df_subset_combined %>% mutate(CHR=as.numeric(CHR),POS=as.integer(POS)) %>% filter(CHR%in%c(1:22)) %>% mutate(CHR=as.integer(CHR))
  
  # Define output file path
  output_file <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_meta_sumstats/", subtype, ".tsv")
  
  # Save as a TSV file
  write.table(df_subset_combined, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

