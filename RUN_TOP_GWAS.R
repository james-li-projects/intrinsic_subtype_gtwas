###########################
# importing libraries
library(readxl)
library(dplyr)
library(TOP)
library(data.table)
###########################

# set seed
set.seed(1)

# importing partition name
args <- commandArgs(trailingOnly = TRUE)
partition_name <- args[1]
partition_prefix <- sub("\\.txt$", "", partition_name)

# extracting traw for partition
system(paste0("plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS --extract /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/partition/variant_lists/",partition_prefix,".txt --export Av --memory 100000 --threads 1 --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/partition/traw/",partition_prefix))

###########################
# reading in traw/dosage file
input_traw<-fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/partition/traw/",partition_prefix,".traw"))

# convert dosage columns to numeric and subtract each value from 2 since the COUNTED allele is actually the reference allele in the input traw, while we want the counted allele to be the alternate allele 
cols_to_update <- 7:ncol(input_traw)
input_traw[, (cols_to_update) := lapply(.SD, function(x) 2 - as.numeric(x)), .SDcols = cols_to_update]

# accordingly, swap the values in COUNTED and ALT columns
input_traw[, c("COUNTED", "ALT") := .(ALT, COUNTED)]

# adjust the formatting of this traw file
traw<-t(input_traw)
rownames(traw) <- gsub("0_","",rownames(traw))
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

# outputting final covariates df
covariates <- combined_cov %>% select(Sample_Name,Status,ER,PR,HER2,dataset_recode,Age_GWAS,paste0(rep("PC",10),c(1:10))) %>% rename(Platform=dataset_recode,Age=Age_GWAS,case.control=Status)

# recoding missing receptor status in cases and controls
covariates <- covariates %>% mutate(ER=ifelse(ER==8,NA,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==8,NA,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==8,NA,HER2))
covariates <- covariates %>% mutate(ER=ifelse(ER==9,888,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==9,888,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==9,888,HER2))
# making receptor negative a 0 value
covariates <- covariates %>% mutate(ER=ifelse(ER==2,0,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==2,0,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==2,0,HER2))
###########################


###########################
# perform intrinsic subtype GWAS analyses for a single SNP at a time
###########################
# print time
a1<-Sys.time()
print(paste("Starting analysis at:",a1))

# identifying platform list to iterate through for generating summary statistics
platform_list <- sort(unique(covariates$Platform))

# initializing current platform to compute summary statistics
for (platform in platform_list) {
  print(paste0("Calculating summary statistics for ",partition_prefix," - ",platform))
  
  # initializing data.frame for output summary statistics
  final_sumstats <- data.frame() 
  
  for (i in 1:ncol(traw)) {
    print(paste0("Evaluating variant: ",i))
    
    # identifying variant name
    current_variant_name <- toString(traw["SNP",i])
    
    # identifying effect allele name
    current_variant_effect_allele <- toString(traw["COUNTED",i])
    
    # identifying non-effect allele name (ALT because we swapped the COUNTED/ALT columns above)
    current_variant_noneffect_allele <- toString(traw["ALT",i])
    
    # identifying CHR and POS
    current_variant_CHR <- as.numeric(trimws(traw["CHR",i]))
    current_variant_POS <- as.numeric(trimws(traw["POS",i]))
    
    # assembling temporary data frame with a SNP and covariates
    data <- cbind(covariates,as.numeric(traw[covariates$Sample_Name,i])) %>% rename(SNP=`as.numeric(traw[covariates$Sample_Name, i])`) %>% select(case.control, ER, PR, HER2, SNP, Age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, Platform) 
    data <- data %>% filter(Platform==platform) %>% select(-Platform)
    
    # filter out individuals with NA dosage values for the current variant
    data<-data%>%filter(!is.na(SNP))
    
    # continue to evaluate suitability for TOP if some individuals have non-NA variant data
    if (nrow(data) > 0) {
      
      # identifying the number of individuals with each subtype
      data$subtype<-NA
      #for first subtype HR+_HER2-
      data$subtype[which((data$ER==1|data$PR==1)
                         &data$HER2==0)] <- "subtype1"
      # for second subtype HR+_HER2+
      data$subtype[which((data$ER==1|data$PR==1)
                         &data$HER2==1)] <- "subtype2"
      # for third subtype HR-_HER2+
      data$subtype[which(data$ER==0&data$PR==0
                         &data$HER2==1)] <- "subtype3"
      # for third subtype HR-_HER2-
      data$subtype[which(data$ER==0&data$PR==0
                         &data$HER2==0)] <- "subtype4"
      
      ###########################
      # proceed only if there are >10 individuals (cases and controls) for a single subtype with non-NA variant values
      if ((sum(table(data$subtype, useNA = "no") > 10) >= 1) & (nrow(data %>% filter(case.control==0)) > 10)) {
        
        ###########################
        # initializing available subtype list
        subtypes_to_consider <- c(
          "subtype1",
          "subtype2",
          "subtype3",
          "subtype4")
        z.design.cols_to_consider <- c(
          1:4
        )
        
        ###########################
        # identifying problematic subtypes for removal: 1) subtypes with no variability in their SNP dosage values, 2) limited variability with two or less non-extreme (0 or 2) dosages. TRUE status for any of the conditions will indicate removal of a subtype. 
        condition_vector<-c()
        for (tmp_subtype_i in 1:4) {
          if (paste0("subtype",tmp_subtype_i) %in% data$subtype) {
            tmp_condition=(var((data %>% filter(subtype==paste0("subtype",tmp_subtype_i)))$SNP)==0) | (is.na(var((data %>% filter(subtype==paste0("subtype",tmp_subtype_i)))$SNP))) | (sum((data.frame(table((data %>% filter(subtype==paste0("subtype",tmp_subtype_i)))$SNP)) %>% filter(!(Var1 %in% c(0,2))))$Freq) < 4) | (sum((data.frame(table((data %>% filter(subtype==paste0("subtype",tmp_subtype_i)))$SNP)) %>% filter(Var1 == 0))$Freq) < 4 & sum((data.frame(table((data %>% filter(subtype==paste0("subtype",tmp_subtype_i)))$SNP)) %>% filter(Var1 == 2))$Freq) < 4) | nrow((data %>% filter(subtype==paste0("subtype",tmp_subtype_i)))) <= 10 | var((data %>% filter(subtype==paste0("subtype",tmp_subtype_i)))$SNP) == 0
            condition_vector <- c(condition_vector, tmp_condition)
          } else {
            tmp_condition=TRUE
            condition_vector <- c(condition_vector, tmp_condition)
          }
        }

        ###########################
        # removing subtypes that were filtered out as being problematic
        if (condition_vector[1]==TRUE) {
          subtypes_to_consider[1]<-NA
          z.design.cols_to_consider[1]<-NA
        }
        if (condition_vector[2]==TRUE) {
          subtypes_to_consider[2]<-NA
          z.design.cols_to_consider[2]<-NA
        }
        if (condition_vector[3]==TRUE) {
          subtypes_to_consider[3]<-NA
          z.design.cols_to_consider[3]<-NA
        }
        if (condition_vector[4]==TRUE) {
          subtypes_to_consider[4]<-NA
          z.design.cols_to_consider[4]<-NA
        }
        subtypes_to_consider <- subtypes_to_consider[!is.na(subtypes_to_consider)]
        z.design.cols_to_consider <- z.design.cols_to_consider[!is.na(z.design.cols_to_consider)]
        
        ###########################
        # proceed with running TOP if at least two subtypes remain after filtering
        if ((length(subtypes_to_consider) >= 2)) {

          # removing subtypes with no SNP variability
          data_update = data %>% 
            filter(subtype%in%c(subtypes_to_consider,NA))
          
          # obtaining sample sizes for each subtype
          n_subtype1=nrow(data_update %>% filter(subtype=="subtype1"))
          n_subtype2=nrow(data_update %>% filter(subtype=="subtype2"))
          n_subtype3=nrow(data_update %>% filter(subtype=="subtype3"))
          n_subtype4=nrow(data_update %>% filter(subtype=="subtype4"))
          n_case=n_subtype1+n_subtype2+n_subtype3+n_subtype4
          n_control=nrow(data_update %>% filter(case.control==0))
          
          # computing minor allele frequency
          control.EAF_snp <- mean((data_update %>% filter(case.control==0))$SNP) / 2
          case.EAF_snp <- mean((data_update %>% filter(case.control==1))$SNP) / 2

          # removing the temporary "subtype" column
          data_update$subtype <- NULL
          
          ###########################
          # prepping inputs for TOP
          y<-data.matrix(data_update %>% select(case.control, ER, PR, HER2))
          z.standard <- GenerateZstandard(y)
          M <- nrow(z.standard) 
          z.design <- matrix(0,M,4)
          colnames(z.design) <- c(
            "HR+_HER2-",
            "HR+_HER2+",
            "HR-_HER2+",
            "HR-_HER2-"
          )
          
          # constructing design matrix for subtypes
          #for first subtype HR+_HER2-
          idx.1 <- which((z.standard[,1]==1|z.standard[,2]==1)
                         &z.standard[,3]==0)
          z.design[idx.1,1] <- 1
          # for second subtype HR+_HER2+
          idx.2 <- which((z.standard[,1]==1|z.standard[,2]==1)
                         &z.standard[,3]==1)
          z.design[idx.2,2] <- 1
          # for third subtype HR-_HER2+
          idx.3 <- which(z.standard[,1]==0&z.standard[,2]==0
                         &z.standard[,3]==1)
          z.design[idx.3,3] <- 1
          # for third subtype HR-_HER2-
          idx.4 <- which(z.standard[,1]==0&z.standard[,2]==0
                         &z.standard[,3]==0)
          z.design[idx.4,4] <- 1
          
          # one SNP and one Principal components (PC1) are the covariates
          SNP <- data_update[,5,drop=F]
          COV_MAT <- data_update[,6:ncol(data_update),drop=F]
          
          ###########################
          # check if TOP will converge for the current variant
          TOP_CONVERGE_STATUS <- tryCatch({
            # Main code block
            model.3 <- EMmvpolySelfDesign(
              y,
              x.self.design = data.matrix(SNP),
              z.design = cbind(z.design[,c(z.design.cols_to_consider)]),
              baselineonly = COV_MAT,
              missingTumorIndicator = 888
            )
            Sys.time()  # Capture system time if the code runs successfully
            
            # If successful, print a confirmation message
            print("CONVERGED")
          }, error = function(e) {
            # If there is an error, print "failed"
            print("FAILED")
          })
          
          # run obtain TOP results if TOP succesfully converged, otherwise output NA results
          if (TOP_CONVERGE_STATUS == "CONVERGED") {
            ###########################
            # running the self design TOP regression
            Sys.time()
            model.3 <- EMmvpolySelfDesign(
              y,
              x.self.design = data.matrix(SNP),
              z.design = cbind(z.design[,c(z.design.cols_to_consider)]),
              baselineonly=COV_MAT,
              missingTumorIndicator = 888)
            Sys.time()
            
            ###########################
            # parsing a data.frame for the output 
            raw_output <- data.frame(model.3[[4]]) %>% select(-OddsRatio,-`OddsRatio.95.CI.low.`,-`OddsRatio.95.CI.high.`,-Covariate)
            # computing parameters to extract logOR and SE
            M <- nrow(z.standard)
            C = ncol(z.design[,c(z.design.cols_to_consider)])
            # extracting logOR and SE
            raw_output$BETA<-model.3[[1]][(M+1):(M+C)]
            raw_output$SE<-sqrt(diag(model.3[[2]][(M+1):(M+C),(M+1):(M+C)]))
            # rearranging columns
            raw_output <- raw_output %>% select(SecondStageEffect,BETA,SE,Pvalue)
            
            # Replace occurrences in the SecondStageEffect column
            raw_output$SecondStageEffect <- gsub("HR\\+", "HRPOS", raw_output$SecondStageEffect)
            raw_output$SecondStageEffect <- gsub("HR\\-", "HRNEG", raw_output$SecondStageEffect)
            raw_output$SecondStageEffect <- gsub("HER2\\+", "HER2POS", raw_output$SecondStageEffect)
            raw_output$SecondStageEffect <- gsub("HER2\\-", "HER2NEG", raw_output$SecondStageEffect)
            
            # parsing the output for each subtype
            if (nrow(raw_output %>% filter(SecondStageEffect=="HRPOS_HER2NEG"))>0) {
              HRPOS_HER2NEG_vector <- as.numeric(raw_output %>% filter(SecondStageEffect=="HRPOS_HER2NEG") %>% select(-SecondStageEffect))
            } else {
              HRPOS_HER2NEG_vector <- c(rep(NA,3))
            }
            
            if (nrow(raw_output %>% filter(SecondStageEffect=="HRPOS_HER2POS"))>0) {
              HRPOS_HER2POS_vector <- as.numeric(raw_output %>% filter(SecondStageEffect=="HRPOS_HER2POS") %>% select(-SecondStageEffect))
            } else {
              HRPOS_HER2POS_vector <- c(rep(NA,3))
            }
            
            if (nrow(raw_output %>% filter(SecondStageEffect=="HRNEG_HER2POS"))>0) {
              HRNEG_HER2POS_vector <- as.numeric(raw_output %>% filter(SecondStageEffect=="HRNEG_HER2POS") %>% select(-SecondStageEffect))
            } else {
              HRNEG_HER2POS_vector <- c(rep(NA,3))
            }
            
            if (nrow(raw_output %>% filter(SecondStageEffect=="HRNEG_HER2NEG"))>0) {
              HRNEG_HER2NEG_vector <- as.numeric(raw_output %>% filter(SecondStageEffect=="HRNEG_HER2NEG") %>% select(-SecondStageEffect))
            } else {
              HRNEG_HER2NEG_vector <- c(rep(NA,3))
            }
            
            # storing output for each subtype in a dataframe
            transposed_df <- data.frame(t(c(
              HRPOS_HER2NEG_vector,
              HRPOS_HER2POS_vector,
              HRNEG_HER2POS_vector,
              HRNEG_HER2NEG_vector)))
            colnames(transposed_df) <- c(
              "HRPOS_HER2NEG_BETA",
              "HRPOS_HER2NEG_SE",
              "HRPOS_HER2NEG_P",
              "HRPOS_HER2POS_BETA",
              "HRPOS_HER2POS_SE",
              "HRPOS_HER2POS_P",
              "HRNEG_HER2POS_BETA",
              "HRNEG_HER2POS_SE",
              "HRNEG_HER2POS_P",
              "HRNEG_HER2NEG_BETA",
              "HRNEG_HER2NEG_SE",
              "HRNEG_HER2NEG_P"
            )
            
            # additionally including the variant_id and extract columns for overall global tests of association and heterogeneity
            output_df <- cbind(
              variant_id=current_variant_name, 
              CHR=current_variant_CHR,
              POS=current_variant_POS,
              EffectAllele=current_variant_effect_allele, 
              BaselineAllele=current_variant_noneffect_allele,
              control.EAF=control.EAF_snp,
              case.EAF=case.EAF_snp,
              CONTROL_N=n_control,
              CASE_N=n_case,
              HRPOS_HER2NEG_N=n_subtype1,
              HRPOS_HER2POS_N=n_subtype2,
              HRNEG_HER2POS_N=n_subtype3,
              HRNEG_HER2NEG_N=n_subtype4,
              transposed_df)
            rownames(output_df) <- NULL
            
            # storing results in a final results data.frame 
            final_sumstats <- rbind(final_sumstats,output_df)
          } else { 
            # assembling row with all NA values for a given variant (output_df)
            transposed_df <- data.frame(t(c(rep(NA,12))))
            colnames(transposed_df) <- c(
              "HRPOS_HER2NEG_BETA",
              "HRPOS_HER2NEG_SE",
              "HRPOS_HER2NEG_P",
              "HRPOS_HER2POS_BETA",
              "HRPOS_HER2POS_SE",
              "HRPOS_HER2POS_P",
              "HRNEG_HER2POS_BETA",
              "HRNEG_HER2POS_SE",
              "HRNEG_HER2POS_P",
              "HRNEG_HER2NEG_BETA",
              "HRNEG_HER2NEG_SE",
              "HRNEG_HER2NEG_P"
            )
            output_df <- cbind(
              variant_id=current_variant_name, 
              CHR=current_variant_CHR,
              POS=current_variant_POS,
              EffectAllele=current_variant_effect_allele, 
              BaselineAllele=current_variant_noneffect_allele,
              control.EAF=NA,
              case.EAF=NA,
              CONTROL_N=NA,
              CASE_N=NA,
              HRPOS_HER2NEG_N=NA,
              HRPOS_HER2POS_N=NA,
              HRNEG_HER2POS_N=NA,
              HRNEG_HER2NEG_N=NA,
              transposed_df)
            rownames(output_df) <- NULL
            
            # storing results in a final results data.frame 
            final_sumstats <- rbind(final_sumstats,output_df)
          }

        } else if ((length(subtypes_to_consider) < 2)) {
          
          print("Insufficient number of subtypes to run TOP")
          
          # assembling row with all NA values for a given variant (output_df)
          transposed_df <- data.frame(t(c(rep(NA,12))))
          colnames(transposed_df) <- c(
            "HRPOS_HER2NEG_BETA",
            "HRPOS_HER2NEG_SE",
            "HRPOS_HER2NEG_P",
            "HRPOS_HER2POS_BETA",
            "HRPOS_HER2POS_SE",
            "HRPOS_HER2POS_P",
            "HRNEG_HER2POS_BETA",
            "HRNEG_HER2POS_SE",
            "HRNEG_HER2POS_P",
            "HRNEG_HER2NEG_BETA",
            "HRNEG_HER2NEG_SE",
            "HRNEG_HER2NEG_P"
          )
          output_df <- cbind(
            variant_id=current_variant_name, 
            CHR=current_variant_CHR,
            POS=current_variant_POS,
            EffectAllele=current_variant_effect_allele, 
            BaselineAllele=current_variant_noneffect_allele,
            control.EAF=NA,
            case.EAF=NA,
            CONTROL_N=NA,
            CASE_N=NA,
            HRPOS_HER2NEG_N=NA,
            HRPOS_HER2POS_N=NA,
            HRNEG_HER2POS_N=NA,
            HRNEG_HER2NEG_N=NA,
            transposed_df)
          rownames(output_df) <- NULL
          
          # storing results in a final results data.frame 
          final_sumstats <- rbind(final_sumstats,output_df)
        }
        
      } else {
        print("No single subtype has at least 10 cases and 10 controls with non-NA values")
        # assembling row with all NA values for a given variant (output_df)
        transposed_df <- data.frame(t(c(rep(NA,12))))
        colnames(transposed_df) <- c(
          "HRPOS_HER2NEG_BETA",
          "HRPOS_HER2NEG_SE",
          "HRPOS_HER2NEG_P",
          "HRPOS_HER2POS_BETA",
          "HRPOS_HER2POS_SE",
          "HRPOS_HER2POS_P",
          "HRNEG_HER2POS_BETA",
          "HRNEG_HER2POS_SE",
          "HRNEG_HER2POS_P",
          "HRNEG_HER2NEG_BETA",
          "HRNEG_HER2NEG_SE",
          "HRNEG_HER2NEG_P"
        )
        output_df <- cbind(
          variant_id=current_variant_name, 
          CHR=current_variant_CHR,
          POS=current_variant_POS,
          EffectAllele=current_variant_effect_allele, 
          BaselineAllele=current_variant_noneffect_allele,
          control.EAF=NA,
          case.EAF=NA,
          CONTROL_N=NA,
          CASE_N=NA,
          HRPOS_HER2NEG_N=NA,
          HRPOS_HER2POS_N=NA,
          HRNEG_HER2POS_N=NA,
          HRNEG_HER2NEG_N=NA,
          transposed_df)
        rownames(output_df) <- NULL
        
        # storing results in a final results data.frame 
        final_sumstats <- rbind(final_sumstats,output_df)
      }
      
    } else {
      print("Variant has all missing values")
      # assembling row with all NA values for a given variant (output_df)
      transposed_df <- data.frame(t(c(rep(NA,12))))
      colnames(transposed_df) <- c(
        "HRPOS_HER2NEG_BETA",
        "HRPOS_HER2NEG_SE",
        "HRPOS_HER2NEG_P",
        "HRPOS_HER2POS_BETA",
        "HRPOS_HER2POS_SE",
        "HRPOS_HER2POS_P",
        "HRNEG_HER2POS_BETA",
        "HRNEG_HER2POS_SE",
        "HRNEG_HER2POS_P",
        "HRNEG_HER2NEG_BETA",
        "HRNEG_HER2NEG_SE",
        "HRNEG_HER2NEG_P"
      )
      output_df <- cbind(
        variant_id=current_variant_name, 
        CHR=current_variant_CHR,
        POS=current_variant_POS,
        EffectAllele=current_variant_effect_allele, 
        BaselineAllele=current_variant_noneffect_allele,
        control.EAF=NA,
        case.EAF=NA,
        CONTROL_N=NA,
        CASE_N=NA,
        HRPOS_HER2NEG_N=NA,
        HRPOS_HER2POS_N=NA,
        HRNEG_HER2POS_N=NA,
        HRNEG_HER2NEG_N=NA,
        transposed_df)
      rownames(output_df) <- NULL
      
      # storing results in a final results data.frame 
      final_sumstats <- rbind(final_sumstats,output_df)
    }
  }
  
  # writing out summary statistics for current partition and platform  
  write.table(final_sumstats,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/partition/output_sumstats/",platform,"_",partition_prefix,".sumstats"),row.names=F,col.names=T,quote=F,sep="\t")
  print(paste0("Finished writing out summary statistics for partition/platform: ",partition_prefix,"/",platform))
} 

# printing out end time
a2<-Sys.time()
print(paste("Finished all analyses at:",a2))
print(paste("Total runtime:",a2-a1))

# printing success message if the last SNP (5000) for the last platform (GBHS) was reached
if (i == 5000 & platform == "WGS") {
  print("ANALYSES STATUS: SUCCESSFULLY COMPLETED")
} else {
  print("ANALYSES STATUS: INCOMPLETE")
}

