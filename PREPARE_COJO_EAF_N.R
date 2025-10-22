library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

####################################
############# EUROPEAN #############
####################################
# import EUR summary statistics and filter out variants with invalid alleles
EUR_SUMSTATS <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/EUR_SUMSTATS/BCAC_4subtypes_sum_data.txt") %>% 
  mutate(col4 = as.integer(seq(nrow(.)))) %>% 
  filter(Effect.Meta!="") %>% 
  filter(Baseline.Meta!="") %>%
  filter(grepl("^[ATCG]+$", Effect.Meta) & grepl("^[ATCG]+$", Baseline.Meta)) %>%
  mutate(col4 = as.integer(col4)) %>% 
  mutate(chr_hg19=chr.iCOGs,pos_hg19=Position.iCOGs) 

# computing allele frequency and N
EUR_SUMSTATS <- EUR_SUMSTATS %>% mutate(
  N=NumCalled.Onco+NumCalled.iCOGs,
  EAF=
  (EAFcontrols.Onco*NumCalled.Onco+EAFcontrols.iCOGs*NumCalled.iCOGs) /
    (NumCalled.Onco+NumCalled.iCOGs)
)
EUR_SUMSTATS <- EUR_SUMSTATS %>% 
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
    ID2 = paste(chr_hg38,pos_hg38,BaselineAllele2,EffectAllele2,sep=":")) %>% mutate(EAF1=EAF,EAF2=1-EAF)

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
subtype="HRPOS_HER2NEG"

# importing AFR GWAS pvar file to see the format of aligned variants
pvar <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS.pvar")

# parsing sumstats for each subtype
print(paste("Arranging sumstats for:",subtype)) 

# Select and rename relevant columns
df_subset <- liftOver_EUR_SUMSTATS %>%
  select(ID1, 
         ID2, 
         chr_hg38,
         pos_hg38,
         EffectAllele1, BaselineAllele1, 
         EffectAllele2, BaselineAllele2,
         EAF1, EAF2, N,
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
  EAF1,
  N,
  BETA1,
  SE,
  P
) %>% filter(ID1 %in% pvar$ID) %>% rename(
  ID=ID1,
  EffectAllele=EffectAllele1,
  BaselineAllele=BaselineAllele1,
  EAF=EAF1,
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
  EAF2,
  N,
  BETA2,
  SE,
  P
) %>% filter(ID2 %in% pvar$ID) %>% filter(!(ID1%in%df_subset_1$ID)) %>% select(-ID1) %>% rename(
  ID=ID2,
  EffectAllele=EffectAllele2,
  BaselineAllele=BaselineAllele2,
  EAF=EAF2,
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
  EAF1,
  N,
  BETA1,
  SE,
  P
) %>% filter(!(ID1 %in% c(pvar$ID,df_subset_2_ID1_list))) %>% rename(
  ID=ID1,
  EffectAllele=EffectAllele1,
  BaselineAllele=BaselineAllele1,
  EAF=EAF1,
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

# Save N and EAF as a TSV file
write.table(df_subset_combined %>% select(ID,EAF,N), file = "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/eur_eaf_n.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


###################################
############# AFRICAN #############
###################################
# importing libraries
library(readxl)
library(dplyr)
library(data.table)

# importing covariate data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
# converting breast cancer status to a binary variable
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,`Dataset...3`,Age_GWAS,Status,ER,PR,HER2,AFR_pro)
pheno_data$Status <- 2-pheno_data$Status
colnames(pheno_data)[1] <- "Sample_Name"
colnames(pheno_data)[2] <- "Dataset"

# importing data of principal components
eigenvec <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/PCA/AABCG_PCA_50PCs.eigenvec",header=T)
colnames(eigenvec)[1] <- "Sample_Name"

# joining covariate and PC data.frames by sample IID
combined_cov <- inner_join(pheno_data,eigenvec, by = c("Sample_Name")) 
# recoding the dataset variable 
combined_cov$dataset_recode <- NA
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_01_WGS"] <- "WGS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_02_WGS2"] <- "WGS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_03_MEGA_VANDY"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_04_MEGA_RP"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_05_MEGA_USC"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_06_AMBER"] <- "AMBER" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_07_ROOT"] <- "ROOT" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_08_AABC"] <- "AABC" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_09_GBHS"] <- "GBHS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_10_BCAC_OncoArray"] <- "BCAC_OncoArray" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_11_BCAC_iCOGS"] <- "BCAC_iCOGS" 

# selecting final covariates for analysis
combined_cov <- combined_cov %>% select(Sample_Name,Status,ER,PR,HER2,dataset_recode,Age_GWAS,paste0(rep("PC",10),c(1:10))) %>% rename(Platform=dataset_recode,Age=Age_GWAS,case.control=Status)

# recoding missing receptor status in cases and controls
combined_cov <- combined_cov %>% mutate(ER=ifelse(ER==8,NA,ER))
combined_cov <- combined_cov %>% mutate(PR=ifelse(PR==8,NA,PR))
combined_cov <- combined_cov %>% mutate(HER2=ifelse(HER2==8,NA,HER2))
combined_cov <- combined_cov %>% mutate(ER=ifelse(ER==9,888,ER))
combined_cov <- combined_cov %>% mutate(PR=ifelse(PR==9,888,PR))
combined_cov <- combined_cov %>% mutate(HER2=ifelse(HER2==9,888,HER2))

# making receptor negative a 0 value
combined_cov <- combined_cov %>% mutate(ER=ifelse(ER==2,0,ER))
combined_cov <- combined_cov %>% mutate(PR=ifelse(PR==2,0,PR))
combined_cov <- combined_cov %>% mutate(HER2=ifelse(HER2==2,0,HER2))

# initializing subtype column
combined_cov$subtype <- NA

#for first subtype HR+_HER2-
combined_cov$subtype[which((combined_cov$ER==1|combined_cov$PR==1)&combined_cov$HER2==0)] <- "HRPOS_HER2NEG"
# for second subtype HR+_HER2+
combined_cov$subtype[which((combined_cov$ER==1|combined_cov$PR==1)&combined_cov$HER2==1)] <- "HRPOS_HER2POS"
# for third subtype HR-_HER2+
combined_cov$subtype[which(combined_cov$ER==0&combined_cov$PR==0&combined_cov$HER2==1)] <- "HRNEG_HER2POS"
# for third subtype HR-_HER2-
combined_cov$subtype[which(combined_cov$ER==0&combined_cov$PR==0&combined_cov$HER2==0)] <- "HRNEG_HER2NEG"

# Identifying controls in the subtype column #
combined_cov <- combined_cov %>% mutate(subtype=ifelse(case.control==0,"control",subtype))

# Finalizing covariate column selection #
combined_cov <- combined_cov %>% select(Sample_Name,case.control,Platform,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,subtype)
final_combined_cov <- combined_cov %>% mutate(case.control=case.control + 1) %>% rename(Status=case.control,`#IID`=Sample_Name) 

# filtering DF for just controls
control_cov_df <- combined_cov %>% filter(subtype=="control") %>% mutate(V1=0) %>% select(V1,Sample_Name)

# subsetting plink file for controls
write.table(control_cov_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls.list",quote=F,row.names=F,col.names=F,sep="\t")
system("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS --keep /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls.list --make-pgen --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls")

# computing allele frequencies
system(paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls --freq --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls"))
afr.afreq <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls.afreq")
afr.afreq[, ID_ALT := tstrsplit(ID, "\\:", keep = 4)]
afr.afreq[, EAF := ifelse(ALT == ID_ALT, ALT_FREQS, 1 - ALT_FREQS)]
afr.afreq <- afr.afreq %>% select(ID,EAF)

# computing sample counts with non-missing data
system(paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls --missing variant-only vcols=nmissdosage --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls"))
afr.vmiss <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/aabcg_controls.vmiss") %>% rename(ID=`#ID`) %>% mutate(N=18800-MISSING_DOSAGE_CT) %>% select(-MISSING_DOSAGE_CT)

# Save N and EAF as a TSV file
afr_eaf_n <- inner_join(afr.afreq,afr.vmiss,by=c("ID"))
write.table(afr_eaf_n, file = "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/afr_eaf_n.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
