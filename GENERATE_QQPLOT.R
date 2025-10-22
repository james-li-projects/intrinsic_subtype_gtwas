library(dplyr)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(plotrix)
library(data.table)
library(RColorBrewer)
library(optparse)
library(tidyr)
library(qqman)

# importing list of problematic variants
problematic_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/problematic_variant_list/problematic_variant_list.txt",header=F))$V1

# running for loop
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_sumstats"))

all_lamba_values <- c()  
all_lamba1000_values <- c()  

for (subtype in c("HRNEG_HER2NEG","HRNEG_HER2POS","HRPOS_HER2POS","HRPOS_HER2NEG")) {
  ###################################
  # load afr_sumstats and removing problematic variants
  afr_sumstats<-fread(paste0("./afr/",subtype,".tsv")) %>% filter(!(ID%in%problematic_variant_list))
  
  # Extract p-values from afr_sumstats and remove NAs
  p_values <- as.numeric(afr_sumstats$P)
  p_values[p_values == 0] <- 4.9e-324
  p_values <- p_values[!is.na(p_values)]  # Remove NA p-values
  
  # Compute expected -log10(p) under uniform distribution
  n <- length(p_values)
  expected <- -log10((1:n) / (n + 1))
  
  # Compute observed -log10(p)
  observed <- -log10(sort(p_values))
  
  # Compute genomic inflation factor (lambda)
  chisq <- qchisq(p_values, df = 1, lower.tail = FALSE)  # Convert p-values to chi-square
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)  # Compute inflation factor
  
  # computing afr_N_eff and lambda1000
  afr_sumstats <- afr_sumstats %>% mutate(MAF=ifelse(EAF>0.5,1-EAF,EAF))
  afr_sumstats <- afr_sumstats %>% mutate(updated_SE = SE*sqrt(2*MAF*(1-MAF))) %>% mutate(inv_updated_SE=1/updated_SE)
  afr_N_eff <- mean(1/afr_sumstats$updated_SE^2)    
  lambda_1000 = 1+(lambda-1)*500/afr_N_eff
  
  # Create data frame for plotting
  qqplot_df <- data.frame(Expected = expected, Observed = observed)
  
  # Create QQ plot
  p <- ggplot(qqplot_df, aes(x = Expected, y = Observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_classic() +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    annotate("text", x = max(expected) * 0.1, y = max(observed) * 0.9,
             label = bquote(lambda == .(round(lambda, 3))),
             size = 5, color = "black", hjust = 0) + 
    annotate("text", x = max(expected) * 0.1, y = max(observed) * 0.8,label = bquote(lambda[1000] == .(round(lambda_1000, 3))),size = 5, color = "black", hjust = 0)
  
  # Define output file path
  output_path <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/qqplots/afr_", subtype,".png")
  
  # Save plot
  ggsave(output_path, plot = p, width = 6, height = 6, dpi = 300)
  message("QQ plot saved to: ", output_path)
  
  
  ###################################
  # load eur_sumstats
  eur_sumstats<-fread(paste0("./eur/",subtype,".tsv"))
  
  # Extract p-values from eur_sumstats and remove NAs
  p_values <- as.numeric(eur_sumstats$P)
  p_values[p_values == 0] <- 4.9e-324
  p_values <- p_values[!is.na(p_values)]  # Remove NA p-values
  
  # Compute expected -log10(p) under uniform distribution
  n <- length(p_values)
  expected <- -log10((1:n) / (n + 1))
  
  # Compute observed -log10(p)
  observed <- -log10(sort(p_values))
  
  # Compute genomic inflation factor (lambda)
  chisq <- qchisq(p_values, df = 1, lower.tail = FALSE)  # Convert p-values to chi-square
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)  # Compute inflation factor
  
  # computing eur_N_eff and lambda1000
  eur_sumstats <- eur_sumstats %>% mutate(MAF=ifelse(EAF>0.5,1-EAF,EAF))
  eur_sumstats <- eur_sumstats %>% mutate(updated_SE = SE*sqrt(2*MAF*(1-MAF))) %>% mutate(inv_updated_SE=1/updated_SE)
  eur_N_eff <- mean(1/eur_sumstats$updated_SE^2)    
  lambda_1000 = 1+(lambda-1)*500/eur_N_eff
  
  # Create data frame for plotting
  qqplot_df <- data.frame(Expected = expected, Observed = observed)
  
  # Create QQ plot
  p <- ggplot(qqplot_df, aes(x = Expected, y = Observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_classic() +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    annotate("text", x = max(expected) * 0.1, y = max(observed) * 0.9,
             label = bquote(lambda == .(round(lambda, 3))),
             size = 5, color = "black", hjust = 0) + 
    annotate("text", x = max(expected) * 0.1, y = max(observed) * 0.8,label = bquote(lambda[1000] == .(round(lambda_1000, 3))),size = 5, color = "black", hjust = 0)
  
  # Define output file path
  output_path <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/qqplots/eur_", subtype,".png")
  
  # Save plot
  ggsave(output_path, plot = p, width = 6, height = 6, dpi = 300)
  message("QQ plot saved to: ", output_path)
  
  ###################################
  # load xancestry sumstats
  sumstats<-fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/xancestry/",subtype,".tsv")) %>% filter(!(ID%in%problematic_variant_list))
  
  # Extract p-values from eur_sumstats and remove NAs
  p_values <- as.numeric(sumstats$P)
  p_values[p_values == 0] <- 4.9e-324
  p_values <- p_values[!is.na(p_values)]  # Remove NA p-values
  
  # Compute expected -log10(p) under uniform distribution
  n <- length(p_values)
  expected <- -log10((1:n) / (n + 1))
  
  # Compute observed -log10(p)
  observed <- -log10(sort(p_values))
  
  # Compute genomic inflation factor (lambda)
  chisq <- qchisq(p_values, df = 1, lower.tail = FALSE)  # Convert p-values to chi-square
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)  # Compute inflation factor
  
  # computing lambda 1000
  lambda_1000 = 1+(lambda-1)*500/(eur_N_eff+afr_N_eff)
  
  # Create data frame for plotting
  qqplot_df <- data.frame(Expected = expected, Observed = observed)
  
  # Create QQ plot
  p <- ggplot(qqplot_df, aes(x = Expected, y = Observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_classic() +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    annotate("text", x = max(expected) * 0.1, y = max(observed) * 0.9,label = bquote(lambda == .(round(lambda, 3))),size = 5, color = "black", hjust = 0) + 
    annotate("text", x = max(expected) * 0.1, y = max(observed) * 0.8,label = bquote(lambda[1000] == .(round(lambda_1000, 3))),size = 5, color = "black", hjust = 0)
  
  
  # Define output file path
  output_path <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_sumstats_remove_problematic_variants/qqplots/xancestry_", subtype,".png")
  
  # Save plot
  ggsave(output_path, plot = p, width = 6, height = 6, dpi = 300)
  message("QQ plot saved to: ", output_path)
  
}


