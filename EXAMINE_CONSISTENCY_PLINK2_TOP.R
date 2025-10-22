# examining consistency between TOP and plink2 results
library(data.table)
library(dplyr)
library(tidyverse)
library(tidyr)

# setting directory to get list of summary statistics files
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_plink2/")
all_sumstats_list = list.files()

# initializing data.frames to store results
all_TOP_plink2_compare_effect <- data.frame()
all_TOP_plink2_compare_NA <- data.frame()
all_TOP_plink2_compare_numsig <- data.frame()

# iterating across all sumstats files
for (current_sumstats_file in all_sumstats_list) {
  # importing sumstats for a subtype and platform 
  # current_sumstats_file <- "HRNEG_HER2NEG.MEGA.sumstats"
  # current_sumstats_file <- "HRPOS_HER2NEG.MEGA.sumstats"
  print(paste("Analyzing summary statistics for:",current_sumstats_file))
  
  # importing summary statistics
  sumstats_TOP <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_TOP/",current_sumstats_file)) %>% rename(BETA_TOP=BETA,SE_TOP=SE,P_TOP=P)
  sumstats_plink2 <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_plink2/",current_sumstats_file)) %>% rename(BETA_plink2=BETA,SE_plink2=SE,P_plink2=P)
  
  # joining summary statistics between TOP and plink2 regressions
  joined_sumstats <- inner_join(sumstats_TOP,sumstats_plink2,by=c("ID","EffectAllele","BaselineAllele"))
  
  #########################################
  # IDENTIFYING NUMBER OF NON-NA VARIANTS # 
  #########################################
  # creating a collated effect size column to see the gain in using just standard plink2 logistic regression 
  joined_sumstats <- joined_sumstats %>% mutate(BETA_COLLATED = ifelse(!is.na(BETA_TOP),BETA_TOP,BETA_plink2))
  tmp_TOP_plink2_compare_NA <- cbind(subtype_dataset=current_sumstats_file,num_nonNA_TOP=sum(!is.na(joined_sumstats$BETA_TOP)),num_nonNA_plink2=sum(!is.na(joined_sumstats$BETA_plink2)),num_nonNA_COLLATED=sum(!is.na(joined_sumstats$BETA_COLLATED)))
  # rbinding these results with all non-NA counts for all datasets/subtypes 
  all_TOP_plink2_compare_NA <- rbind(tmp_TOP_plink2_compare_NA,all_TOP_plink2_compare_NA)
  
  ##########################################################
  # IDENTIFYING NUMBER OF GENOME-WIDE SIGNIFICANT VARIANTS #
  ##########################################################
  tmp_TOP_plink2_compare_numsig <- cbind(
    subtype_dataset=current_sumstats_file,
    num_sig_TOP=sum(joined_sumstats$P_TOP < 5e-8, na.rm = TRUE),
    num_sig_plink2=sum(joined_sumstats$P_plink2 < 5e-8, na.rm = TRUE))
  all_TOP_plink2_compare_numsig <- rbind(
    all_TOP_plink2_compare_numsig,
    tmp_TOP_plink2_compare_numsig
  )
  
  ################
  # ALL VARIANTS #
  ################
  # identify correlation between the effect sizes for each dataset
  lm_model <- lm(BETA_TOP ~ BETA_plink2, data = joined_sumstats)
  summary_lm <- summary(lm_model)
  
  # Extract R-squared, p-value, and effect size estimate
  r_squared <- summary_lm$r.squared
  p_value <- summary_lm$coefficients[2,4]  # P-value for BETA_plink2
  effect_size <- summary_lm$coefficients[2,1]  # Effect size estimate
  
  # storing results
  tmp_TOP_plink2_comparison <- cbind(subtype_dataset=current_sumstats_file,effect_size=effect_size,p=p_value,r2=r_squared,variant_subset="All Variants")
  # binding results to a cumulative results df
  all_TOP_plink2_compare_effect <- rbind(
    tmp_TOP_plink2_comparison,
    all_TOP_plink2_compare_effect
  )
  
  ################
  # TOP (p<5e-8) #
  ################
  if (nrow(joined_sumstats %>% filter(P_TOP<5e-8))>2) {
    tmp_joined_sumstats_subset <- joined_sumstats %>% filter(P_TOP<5e-8)
    if (sum(!is.na(tmp_joined_sumstats_subset$BETA_plink2))>2) {

    # identify correlation between the effect sizes for each dataset
    lm_model <- lm(BETA_TOP ~ BETA_plink2, data = joined_sumstats %>% filter(P_TOP<5e-8))
    summary_lm <- summary(lm_model)
    
    # Extract R-squared, p-value, and effect size estimate
    r_squared <- summary_lm$r.squared
    p_value <- summary_lm$coefficients[2,4]  # P-value for BETA_plink2
    effect_size <- summary_lm$coefficients[2,1]  # Effect size estimate
    
    # storing results
    tmp_TOP_plink2_comparison <- cbind(subtype_dataset=current_sumstats_file,effect_size=effect_size,p=p_value,r2=r_squared,variant_subset="TOP (p<5e-8)")
    # binding results to a cumulative results df
    all_TOP_plink2_compare_effect <- rbind(
      tmp_TOP_plink2_comparison,
      all_TOP_plink2_compare_effect
    )
    }
    
  }
  
  if (nrow(joined_sumstats %>% filter(P_plink2<5e-8))>2) {
    tmp_joined_sumstats_subset <- joined_sumstats %>% filter(P_plink2<5e-8)
    if (sum(!is.na(tmp_joined_sumstats_subset$BETA_TOP))>2) {
      
    ###################
    # plink2 (p<5e-8) #
    ###################
    # identify correlation between the effect sizes for each dataset
    lm_model <- lm(BETA_TOP ~ BETA_plink2, data = joined_sumstats %>% filter(P_plink2<5e-8))
    summary_lm <- summary(lm_model)
    
    # Extract R-squared, p-value, and effect size estimate
    r_squared <- summary_lm$r.squared
    p_value <- summary_lm$coefficients[2,4]  # P-value for BETA_plink2
    effect_size <- summary_lm$coefficients[2,1]  # Effect size estimate
    
    # storing results
    tmp_TOP_plink2_comparison <- cbind(subtype_dataset=current_sumstats_file,effect_size=effect_size,p=p_value,r2=r_squared,variant_subset="plink2 (p<5e-8)")
    # binding results to a cumulative results df
    all_TOP_plink2_compare_effect <- rbind(
      tmp_TOP_plink2_comparison,
      all_TOP_plink2_compare_effect
    )
    }
  }
}

save(all_TOP_plink2_compare_effect,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/all_TOP_plink2_compare_effect.RData")
save(all_TOP_plink2_compare_NA,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/all_TOP_plink2_compare_NA.RData")
save(all_TOP_plink2_compare_numsig,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/all_TOP_plink2_compare_numsig.RData")

# plotting QC plots
load("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/all_TOP_plink2_compare_effect.RData")
load("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/all_TOP_plink2_compare_NA.RData")

all_TOP_plink2_compare_effect <- all_TOP_plink2_compare_effect %>% mutate(subtype_dataset=gsub("\\.sumstats","",subtype_dataset)) %>% separate(subtype_dataset,into=c("subtype","dataset"),sep="\\.")
all_TOP_plink2_compare_effect <- all_TOP_plink2_compare_effect %>% mutate(
  effect_size=as.numeric(effect_size),
  p=as.numeric(p),
  r2=as.numeric(r2),
)

all_TOP_plink2_compare_NA <- all_TOP_plink2_compare_NA %>% mutate(subtype_dataset=gsub("\\.sumstats","",subtype_dataset)) %>% separate(subtype_dataset,into=c("subtype","dataset"),sep="\\.")
all_TOP_plink2_compare_NA <- all_TOP_plink2_compare_NA %>% mutate(
  num_nonNA_TOP=as.numeric(num_nonNA_TOP),
  num_nonNA_plink2=as.numeric(num_nonNA_plink2),
  num_nonNA_COLLATED=as.numeric(num_nonNA_COLLATED),
)

all_TOP_plink2_compare_numsig <- all_TOP_plink2_compare_numsig %>% mutate(subtype_dataset=gsub("\\.sumstats","",subtype_dataset)) %>% separate(subtype_dataset,into=c("subtype","dataset"),sep="\\.")
all_TOP_plink2_compare_numsig <- all_TOP_plink2_compare_numsig %>% mutate(
  num_sig_TOP=as.numeric(num_sig_TOP),
  num_sig_plink2=as.numeric(num_sig_plink2)
)




# QC plot for effect sizes
library(ggplot2)
# Create scatterplot
p <- ggplot(all_TOP_plink2_compare_effect, aes(x = dataset, y = effect_size, color = variant_subset)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_wrap(~ subtype) +
  theme_classic() +
  labs(title = "Effect Size Comparison",
       x = "Dataset",
       y = "Effect Size",
       color = "Variant Subset") +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", size = 1) 

# Save as high-quality PNG
ggsave("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/effect_compare.png",
       plot = p, width = 10, height = 6, dpi = 300)


# QC plot for number of variants with non-NA values
library(ggplot2)
library(reshape2)

# Reshape data into long format for plotting
all_TOP_plink2_compare_NA_long <- melt(
  all_TOP_plink2_compare_NA,
  id.vars = c("subtype", "dataset"),
  measure.vars = c("num_nonNA_TOP", "num_nonNA_plink2", "num_nonNA_COLLATED"),
  variable.name = "category",
  value.name = "num_nonNA"
)

# Create scatterplot
p <- ggplot(all_TOP_plink2_compare_NA_long, aes(x = dataset, y = num_nonNA, color = category)) +
  geom_point(size = 3, alpha = 0.5) +
  facet_wrap(~ subtype) +
  theme_classic() +
  labs(title = "Number of Non-NA Values Across Datasets",
       x = "Dataset",
       y = "Number of Non-NA Values",
       color = "Category") +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Save as high-quality PNG
ggsave("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/effect_num_nonNA.png",
       plot = p, width = 10, height = 6, dpi = 1200)


# create barplot for # significant hits
library(ggplot2)
library(reshape2)
df_long <- melt(all_TOP_plink2_compare_numsig, id.vars = c("subtype", "dataset"),
                measure.vars = c("num_sig_TOP", "num_sig_plink2"),
                variable.name = "Method", value.name = "Significant_Counts")

# Plot with increased font sizes
p <- ggplot(df_long, aes(x = dataset, y = Significant_Counts, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +  # Side-by-side bars
  facet_wrap(~ subtype, scales = "free_y") +  # Facet by subtype
  theme_classic(base_size = 20) +  # Increase overall font size
  labs(
    x = "Dataset",
    y = "Number of Significant Variants",
    fill = "Method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Rotate x-axis labels and increase size
    axis.text.y = element_text(size = 18),  # Increase y-axis text size
    axis.title.x = element_text(size = 22, face = "bold"),  # Increase x-axis title size
    axis.title.y = element_text(size = 22, face = "bold"),  # Increase y-axis title size
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),  # Title size and centering
    legend.text = element_text(size = 18),  # Increase legend text size
    legend.title = element_text(size = 20, face = "bold"),  # Increase legend title size
    legend.position="top"
  )

# Save as a high-quality PNG (1200 DPI)
ggsave("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/num_significant.png",
       plot = p, width = 12, height = 10, dpi = 1200)