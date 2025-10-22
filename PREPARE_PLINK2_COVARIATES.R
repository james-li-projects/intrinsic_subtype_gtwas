###########################
# importing libraries
library(readxl)
library(dplyr)
library(TOP)
library(data.table)
###########################


###########################
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


####################################################
# Assigning intrinsic subtype indices to each case #
####################################################
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

##############################################
# Identifying controls in the subtype column #
##############################################
combined_cov <- combined_cov %>% mutate(subtype=ifelse(case.control==0,"control",subtype))

#########################################
# Finalizing covariate column selection #
#########################################
combined_cov <- combined_cov %>% select(Sample_Name,case.control,Platform,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,subtype)
final_combined_cov <- combined_cov %>% mutate(case.control=case.control + 1) %>% rename(Status=case.control,`#IID`=Sample_Name) 

######################################
# Writing out subset pheno_cov files #
######################################
# extracting covariates for a given subtype
HRPOS_HER2NEG_df <- final_combined_cov %>% filter(subtype %in% c("HRPOS_HER2NEG","control")) %>% select(-subtype)
HRPOS_HER2POS_df <- final_combined_cov %>% filter(subtype %in% c("HRPOS_HER2POS","control")) %>% select(-subtype)
HRNEG_HER2POS_df <- final_combined_cov %>% filter(subtype %in% c("HRNEG_HER2POS","control")) %>% select(-subtype)
HRNEG_HER2NEG_df <- final_combined_cov %>% filter(subtype %in% c("HRNEG_HER2NEG","control")) %>% select(-subtype)

# further extracting covariates for a given dataset and writing out these subtype-specific covariate files
for (current_Platform in unique(HRPOS_HER2NEG_df$Platform)) {
  print(paste("Writing out covariate files for platform:",current_Platform))
  write.table(HRPOS_HER2NEG_df %>% filter(Platform==current_Platform) %>% select(-Platform),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pheno_cov/HRPOS_HER2NEG.",current_Platform,".pheno_cov"),quote=F,row.names=F,col.names=T)
  write.table(HRPOS_HER2POS_df %>% filter(Platform==current_Platform) %>% select(-Platform),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pheno_cov/HRPOS_HER2POS.",current_Platform,".pheno_cov"),quote=F,row.names=F,col.names=T)
  write.table(HRNEG_HER2POS_df %>% filter(Platform==current_Platform) %>% select(-Platform),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pheno_cov/HRNEG_HER2POS.",current_Platform,".pheno_cov"),quote=F,row.names=F,col.names=T)
  write.table(HRNEG_HER2NEG_df %>% filter(Platform==current_Platform) %>% select(-Platform),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pheno_cov/HRNEG_HER2NEG.",current_Platform,".pheno_cov"),quote=F,row.names=F,col.names=T)
}

#############################################
# Writing out sample lists for each subtype #
#############################################
HRPOS_HER2NEG_sample.list <- HRPOS_HER2NEG_df %>% mutate(V1=0) %>% select(V1,`#IID`)
write.table(HRPOS_HER2NEG_sample.list,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/sample_list/HRPOS_HER2NEG_sample.list",quote=F,row.names=F,col.names=F,sep="\t")

HRPOS_HER2POS_sample.list <- HRPOS_HER2POS_df %>% mutate(V1=0) %>% select(V1,`#IID`)
write.table(HRPOS_HER2POS_sample.list,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/sample_list/HRPOS_HER2POS_sample.list",quote=F,row.names=F,col.names=F,sep="\t")

HRNEG_HER2POS_sample.list <- HRNEG_HER2POS_df %>% mutate(V1=0) %>% select(V1,`#IID`)
write.table(HRNEG_HER2POS_sample.list,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/sample_list/HRNEG_HER2POS_sample.list",quote=F,row.names=F,col.names=F,sep="\t")

HRNEG_HER2NEG_sample.list <- HRNEG_HER2NEG_df %>% mutate(V1=0) %>% select(V1,`#IID`)
write.table(HRNEG_HER2NEG_sample.list,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/sample_list/HRNEG_HER2NEG_sample.list",quote=F,row.names=F,col.names=F,sep="\t")

