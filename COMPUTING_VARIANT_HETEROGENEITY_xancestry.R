# Load required libraries
library(data.table)
library(dplyr)

# importing list of lead variants
HRPOS_HER2POS_xancestry_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRPOS_HER2POS.txt") %>% filter(Ancestry=="xancestry"))$ID
HRPOS_HER2NEG_xancestry_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRPOS_HER2NEG.txt") %>% filter(Ancestry=="xancestry"))$ID
HRNEG_HER2NEG_xancestry_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRNEG_HER2NEG.txt") %>% filter(Ancestry=="xancestry"))$ID
xancestry_lead_variant_list <- c(
  HRPOS_HER2POS_xancestry_lead_variant_list,
  HRPOS_HER2NEG_xancestry_lead_variant_list,
  HRNEG_HER2NEG_xancestry_lead_variant_list
)

# Read and process each subtype's summary statistics 1
HRPOS_HER2NEG <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_xancestry_meta_sumstats/HRPOS_HER2NEG.tsv") %>%
  select(ID, EffectAllele, BaselineAllele, BETA, SE) %>%
  rename(HRPOS_HER2NEG_BETA = BETA, HRPOS_HER2NEG_SE = SE) %>%
  filter(ID %in% xancestry_lead_variant_list)
HRPOS_HER2POS <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_xancestry_meta_sumstats/HRPOS_HER2POS.tsv") %>%
  select(ID, EffectAllele, BaselineAllele, BETA, SE) %>%
  rename(HRPOS_HER2POS_BETA = BETA, HRPOS_HER2POS_SE = SE) %>%
  filter(ID %in% xancestry_lead_variant_list)
HRNEG_HER2NEG <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_xancestry_meta_sumstats/HRNEG_HER2NEG.tsv") %>%
  select(ID, EffectAllele, BaselineAllele, BETA, SE) %>%
  rename(HRNEG_HER2NEG_BETA = BETA, HRNEG_HER2NEG_SE = SE) %>%
  filter(ID %in% xancestry_lead_variant_list)
HRNEG_HER2POS <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_xancestry_meta_sumstats/HRNEG_HER2POS.tsv") %>%
  select(ID, EffectAllele, BaselineAllele, BETA, SE) %>%
  rename(HRNEG_HER2POS_BETA = BETA, HRNEG_HER2POS_SE = SE) %>%
  filter(ID %in% xancestry_lead_variant_list)

# Join all dataframes by ID, EffectAllele, and BaselineAllele
joined_df <- HRPOS_HER2NEG %>%
  full_join(HRPOS_HER2POS, by = c("ID", "EffectAllele", "BaselineAllele")) %>%
  full_join(HRNEG_HER2NEG, by = c("ID", "EffectAllele", "BaselineAllele")) %>%
  full_join(HRNEG_HER2POS, by = c("ID", "EffectAllele", "BaselineAllele"))

# initializing list to store covariance matrices
covariance_matrix_list <- vector(mode="list")
# initializing df to store heterogeneity test results
variant_pval_df <- data.frame()

# iterating for each variant
for (i in 1:nrow(joined_df)) {
  
  beta_estimates <- c(
    joined_df$HRPOS_HER2NEG_BETA[i],
    joined_df$HRPOS_HER2POS_BETA[i],
    joined_df$HRNEG_HER2POS_BETA[i],
    joined_df$HRNEG_HER2NEG_BETA[i]
  )
  
  se_estimates <- c(
    joined_df$HRPOS_HER2NEG_SE[i],
    joined_df$HRPOS_HER2POS_SE[i],
    joined_df$HRNEG_HER2POS_SE[i],
    joined_df$HRNEG_HER2NEG_SE[i]
  )
  
  covariance_matrix <- diag(se_estimates^2)
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
  
  # Calculate the p-value
  p_value <- pchisq(test_statistic, df = df, lower.tail = FALSE)
  
  # storing results
  variant_pval_df <- rbind(
    variant_pval_df,
    data.frame(
      variant = joined_df$ID[i],
      p_value = p_value,
      stringsAsFactors = FALSE
    ))
  
  # Print results
  print(joined_df$ID[i])
  cat("Chi-square Test Statistic:", test_statistic, "\n")
  cat("Degrees of Freedom:", df, "\n")
  cat("P-value:", p_value, "\n")
}

# printing variant heterogeneity data.frame and correlation matrices
correlation_matrix_list <- lapply(covariance_matrix_list, cov2cor)
print(variant_pval_df)
print(correlation_matrix_list)

# save het p-value output
write.table(variant_pval_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/subtype_heterogeneity/het_p_xancestry.tsv",quote=F,row.names=F,col.names=T,sep="\t")
