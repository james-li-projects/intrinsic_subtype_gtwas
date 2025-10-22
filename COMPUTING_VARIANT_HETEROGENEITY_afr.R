###########################
# importing libraries
library(readxl)
library(dplyr)
library(TOP)
library(data.table)
###########################

# set seed
set.seed(1)


# defining fixed effects meta-analysis function
fixed_effect_meta_analysis <- function(beta_list, cov_list) {
  # Ensure we have the same number of studies
  K <- length(beta_list)
  if (K != length(cov_list)) {
    stop("Number of beta vectors and covariance matrices must be the same!")
  }
  
  # Compute the sum of inverse covariance matrices
  sum_inv_cov <- Reduce("+", lapply(cov_list, solve))
  
  # Compute the pooled estimate
  pooled_estimate <- solve(sum_inv_cov) %*% Reduce("+", mapply(function(beta, cov) {
    solve(cov) %*% beta
  }, beta_list, cov_list, SIMPLIFY = FALSE))
  
  # Compute covariance of pooled estimate
  pooled_cov <- solve(sum_inv_cov)
  
  # Compute standard errors (square root of diagonal elements)
  se <- sqrt(diag(pooled_cov))
  
  # Compute 95% confidence intervals
  ci_lower <- pooled_estimate - 1.96 * se
  ci_upper <- pooled_estimate + 1.96 * se
  ci <- cbind(ci_lower, ci_upper)
  
  # Return results
  return(list(Estimate = pooled_estimate, SE = se, CI = ci, Covariance = pooled_cov))
}

# reading in lead variants specific to AFR
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/")
HRNEG_HER2NEG <- fread("HRNEG_HER2NEG.txt") %>% filter(Ancestry=="afr") %>% select(ID)
HRPOS_HER2NEG <- fread("HRPOS_HER2NEG.txt") %>% filter(Ancestry=="afr") %>% select(ID)
HRPOS_HER2POS <- fread("HRPOS_HER2POS.txt") %>% filter(Ancestry=="afr") %>% select(ID)
lead_variant_list <- rbind(HRNEG_HER2NEG,HRPOS_HER2NEG,HRPOS_HER2POS) %>% unique()
write.table(lead_variant_list,file="/scratch/jll1/XANCESTRY_GWAS_TWAS/lead_variant/variant.list",quote=F,row.names=F,col.names=F)


# extracting traw for partition
system(paste0("plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS --extract /scratch/jll1/XANCESTRY_GWAS_TWAS/lead_variant/variant.list --export Av --memory 100000 --threads 1 --out /scratch/jll1/XANCESTRY_GWAS_TWAS/lead_variant/dosage"))

###########################
# reading in traw/dosage file
input_traw<-fread(paste0("/scratch/jll1/XANCESTRY_GWAS_TWAS/lead_variant/dosage.traw"))

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

variant_pval_df <- data.frame()

# identifying platform list to iterate through for generating summary statistics
platform_list <- sort(unique(covariates$Platform))

# initializing list to store covariance matrices
covariance_matrix_list <- vector(mode="list")

for (i in 1:ncol(traw)) {
  print(paste0("Evaluating variant: ",i))
  
  # Initialize lists to store results for each platform
  beta_list <- list()
  cov_list <- list()
  
  # initializing current platform to compute summary statistics
  for (platform in platform_list) {
    print(paste0("Obtaining beta and covariance matrix for platform: ",platform))
    
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
      # initializing available subtype list
      subtypes_to_consider <- c(
        "subtype1",
        "subtype2",
        "subtype3",
        "subtype4")
      z.design.cols_to_consider <- c(
        1:4
      )
      
      
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
        
        # Store results in lists
        M <- nrow(z.standard)
        C = ncol(z.design)
        #this is the log-odds ratio for five subtypes
        beta_list[[platform]] <- model.3[[1]][(M+1):(M+C)]
        #this is covariance matrix for the log odds ratio of SNPs
        cov_list[[platform]] <- model.3[[2]][(M+1):(M+C),(M+1):(M+C)]                             
      } else {
      }
    }
    
  }
  
  
  # performing fixed effects meta-analysis on these beta_list and cov_list
  result <- fixed_effect_meta_analysis(beta_list, cov_list)
  
  # initializing meta-analyzed beta estimate and covariance matrix
  beta_estimates <- result$Estimate
  covariance_matrix <- result$Covariance
  covariance_matrix_list[[i]] <- covariance_matrix
  
  # Define the contrast matrix for testing differences (3x4)
  contrast_matrix <- matrix(c(
    1, -1,  0,  0,
    1,  0, -1,  0,
    1,  0,  0, -1
  ), nrow = 3, byrow = TRUE)
  
  # Compute the transformed beta vector under the contrasts
  transformed_beta <- contrast_matrix %*% beta_estimates
  
  # Compute the covariance matrix of the transformed betas
  transformed_covariance <- contrast_matrix %*% covariance_matrix %*% t(contrast_matrix)
  
  # Calculate the chi-square test statistic
  test_statistic <- t(transformed_beta) %*% solve(transformed_covariance) %*% transformed_beta
  
  # Degrees of freedom (number of contrasts)
  df <- nrow(contrast_matrix)
  
  # Calculate the heterogeneity test p-value
  p_value <- pchisq(test_statistic, df = df, lower.tail = FALSE)
  
  # storing results
  variant_pval_df <- rbind(
    variant_pval_df,
    data.frame(
      variant = current_variant_name,
      p_value = p_value,
      stringsAsFactors = FALSE
    ))
}

# printing variant heterogeneity data.frame and correlation matrices
correlation_matrix_list <- lapply(covariance_matrix_list, cov2cor)
print(variant_pval_df)
print(correlation_matrix_list)

# save het p-value output
write.table(variant_pval_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/subtype_heterogeneity/het_p_afr.tsv",quote=F,row.names=F,col.names=T,sep="\t")
