library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)

# import european sumstats
EUR_SUMSTATS <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/EUR_SUMSTATS/BCAC_4subtypes_sum_data.txt")

# Rename beta and se columns
colnames(EUR_SUMSTATS) <- gsub("log_or_meta", "BETA", colnames(EUR_SUMSTATS))
colnames(EUR_SUMSTATS) <- gsub("se_meta", "SE", colnames(EUR_SUMSTATS))

# computing p-values for each subtype
EUR_SUMSTATS <- EUR_SUMSTATS %>% 
  mutate(
    Luminal_A_P=2*pnorm(-abs(Luminal_A_BETA/Luminal_A_SE), lower.tail = T),
    Luminal_B_P=2*pnorm(-abs(Luminal_B_BETA/Luminal_B_SE), lower.tail = T),
    HER2_Enriched_P=2*pnorm(-abs(HER2_Enriched_BETA/HER2_Enriched_SE), lower.tail = T),
    Triple_Neg_P=2*pnorm(-abs(Triple_Neg_BETA/Triple_Neg_SE), lower.tail = T)
  )

# Rename subtype names in columns
colnames(EUR_SUMSTATS) <- gsub("Luminal_A", "HRPOS_HER2NEG", colnames(EUR_SUMSTATS))
colnames(EUR_SUMSTATS) <- gsub("Luminal_B", "HRPOS_HER2POS", colnames(EUR_SUMSTATS))
colnames(EUR_SUMSTATS) <- gsub("HER2_Enriched", "HRNEG_HER2POS", colnames(EUR_SUMSTATS))
colnames(EUR_SUMSTATS) <- gsub("Triple_Neg", "HRNEG_HER2NEG", colnames(EUR_SUMSTATS))

  
# modify variant names format from underscores to colons
EUR_SUMSTATS$var_name <- gsub("_", ":", EUR_SUMSTATS$var_name) 




###################################
# clumping variants identified using each approach
system("
for subtype in HRNEG_HER2NEG HRPOS_HER2POS HRNEG_HER2POS HRPOS_HER2NEG
do
  for approach in plink2 TOP
  do
    # performing LD clumping using plink2
    gwas_sumstats=/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/${approach}/METAANALYSIS_${subtype}_1.tbl
    plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS --clump-p1 0.00000005 --clump-p2 0.00000005 --clump-r2 0.5 --clump-kb 500 --clump ${gwas_sumstats} --clump-id-field MarkerName --clump-p-field P-value --clump-unphased --memory 100000 --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/clumped/${approach}/METAANALYSIS_${subtype}_1.tbl
  done
done 
")


library(data.table)
library(dplyr)
library(ggplot2)
library(broom)

# Importing significant hits and clumped variants for each subtype and approach 
clumps_df <- data.frame()
for (subtype in c("HRNEG_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRPOS_HER2NEG")) {
  
  ###################################
  # Import TOP clump lists
  num_clumps_TOP <- nrow(fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/clumped/TOP/METAANALYSIS_",subtype,"_1.tbl.clumps")))
  
  # Import plink2 clump lists
  num_clumps_plink2 <- nrow(fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/clumped/plink2/METAANALYSIS_",subtype,"_1.tbl.clumps")))
  
  # Storing number of clumps as a DF
  tmp_clumps_df <- data.frame(num_clumps=c(num_clumps_plink2,num_clumps_TOP), approach=c("plink2","TOP"), subtype=c(subtype,subtype))
  clumps_df <- rbind(clumps_df, tmp_clumps_df)
  
  ###################################
  # Import TOP sumstats
  tmp_TOP <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/TOP/METAANALYSIS_",subtype,"_1.tbl"))
  colnames(tmp_TOP)[4:ncol(tmp_TOP)] <- paste0("TOP_", colnames(tmp_TOP)[4:ncol(tmp_TOP)])
  
  # Import plink2 sumstats
  tmp_plink2 <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/plink2/METAANALYSIS_",subtype,"_1.tbl"))
  colnames(tmp_plink2)[4:ncol(tmp_plink2)] <- paste0("plink2_", colnames(tmp_plink2)[4:ncol(tmp_plink2)])
  
  # Identify all significant variant names and filter sumstat DFs only for these variants
  tmp_sig_variant_list <- unique(c(
    (tmp_TOP %>% filter(`TOP_P-value` < 5e-8))$MarkerName,
    (tmp_plink2 %>% filter(`plink2_P-value` < 5e-8))$MarkerName
  ))
  tmp_TOP <- tmp_TOP %>% filter(MarkerName %in% tmp_sig_variant_list)
  tmp_plink2 <- tmp_plink2 %>% filter(MarkerName %in% tmp_sig_variant_list)
  
  # Joining these sumstats for variants significant in either dataset
  merged_approach_sumstats <- full_join(tmp_plink2, tmp_TOP, by = c("MarkerName", "Allele1", "Allele2"))
  write.table(merged_approach_sumstats,paste0(subtype,"_merged_approach_meta.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  
  # Perform linear regression
  lm_model <- lm(TOP_Effect ~ plink2_Effect, data = merged_approach_sumstats)
  lm_summary <- summary(lm_model)
  
  # Extract statistics
  beta_coef <- lm_summary$coefficients["plink2_Effect", "Estimate"]
  p_value <- lm_summary$coefficients["plink2_Effect", "Pr(>|t|)"]
  r_squared <- lm_summary$r.squared
  
  ###################################
  # Create scatter plot
  scatterplot_path <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/scatterplot_effectsize_", subtype, ".png")
  
  plot <- ggplot(merged_approach_sumstats, aes(x = plink2_Effect, y = TOP_Effect)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(
      x = expression("Logistic regression (" * beta * ")"),
      y = expression("Two-stage polytomous logistic regression (" * beta * ")")
    ) +
    annotate("text", 
             x = min(merged_approach_sumstats$plink2_Effect, na.rm = TRUE), 
             y = max(merged_approach_sumstats$TOP_Effect, na.rm = TRUE) * 1.05,  # Position above highest point
             label = bquote(beta == .(round(beta_coef, 2)) * 
                              "\n p: " * .(signif(p_value, 3)) * 
                              "\n " * r^2 * " = " * .(round(r_squared, 2))), 
             hjust = 0, size = 5) +  # Match axis text size
    theme_bw(base_size = 14)
  
  ggsave(scatterplot_path, plot, width = 6, height = 6, dpi = 300)
  
  # Create table of consistency and save as a TSV
  consistency_table <- as.data.frame.matrix(table(merged_approach_sumstats$TOP_Direction, merged_approach_sumstats$plink2_Direction))
  consistency_table_path <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/tabulate_study_consistency_", subtype, ".tsv")
  
  write.table(consistency_table, consistency_table_path, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  
  ###################################
  # extracting associations for these significant variants among all the datasets
  gwas_methods <- c("TOP", "plink2")  
  
  # Initialize a list to store merged data for each GWAS method
  merged_list <- list()
  
  for (gwas_method in gwas_methods) {
    
    # Define file path
    file_prefix <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_", gwas_method, "/", subtype)
    
    # List all files that match the pattern
    file_list <- list.files(path = dirname(file_prefix), pattern = paste0("^", basename(file_prefix)), full.names = TRUE)
    
    # Function to process each file
    process_file <- function(file) {
      # Read the file
      dt <- fread(file)
      # Filter based on ID column
      dt <- dt[ID %in% tmp_sig_variant_list]
      # Extract the dataset label from the filename
      if (gwas_method == "TOP") {
        dataset_label <- gsub(paste0(subtype, "\\."), "", gsub("^.*sumstats_TOP/|\\.sumstats$", "", file))
      } else if (gwas_method == "plink2") {
        dataset_label <- gsub(paste0(subtype, "\\."), "", gsub("^.*sumstats_plink2/|\\.sumstats$", "", file))
      } else {}
      
      # Rename columns
      setnames(dt, old = c("BETA", "SE", "P"),
               new = c(paste0(dataset_label, "_BETA_", gwas_method),
                       paste0(dataset_label, "_SE_", gwas_method),
                       paste0(dataset_label, "_P_", gwas_method)))
      return(dt)
    }
    
    # Process all files and store results in a list
    dt_list <- lapply(file_list, process_file)
    
    # Perform full join on ID, EffectAllele, and BaselineAllele
    merged_dt <- Reduce(function(x, y) merge(x, y, by = c("ID", "EffectAllele", "BaselineAllele"), all = TRUE), dt_list)
    
    # Store merged result for this method
    merged_list[[gwas_method]] <- merged_dt
  }
  
  # Perform final full join across all GWAS methods
  final_merged_dt <- Reduce(function(x, y) merge(x, y, by = c("ID", "EffectAllele", "BaselineAllele"), all = TRUE), merged_list)
  
  # Reorder the columns in final_merged_dt
  fixed_cols <- c("ID", "EffectAllele", "BaselineAllele")
  data_cols <- setdiff(names(final_merged_dt), fixed_cols)
  prefixes <- sub("_.*", "", data_cols)
  sorted_data_cols <- data_cols[order(prefixes)]
  final_merged_dt <- final_merged_dt[, c(fixed_cols, sorted_data_cols), with = FALSE]
  write.table(final_merged_dt,file=paste0("merged_effect_approach_dataset_",subtype,".tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  
  #################################
  # create scatterplot of effect sizes for each method and dataset for each variant
  output_file <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/scatterplot_method_", subtype, ".png")
  beta_cols <- grep("_BETA_", names(final_merged_dt), value = TRUE)
  long_dt <- melt(final_merged_dt, id.vars = "ID", measure.vars = beta_cols, variable.name = "Method", value.name = "BETA")
  long_dt[, Suffix := sub(".*_BETA_", "", Method)]
  p <- ggplot(long_dt, aes(x = ID, y = BETA, color = Suffix)) +
    geom_point(alpha = 0.4) +  # Scatter points with some transparency
    geom_hline(yintercept = 0, linetype = "dotted", color = "red", alpha = 0.6) +  # Horizontal reference line
    labs(title = paste("BETA Values by Method for", subtype),
         x = "Variant Index",
         y = "GWAS Effect Size",
         color = "GWAS Method") +
    theme_classic() +
    theme(axis.text.x = element_blank(), # Hide x-axis text for readability
          axis.ticks.x = element_blank())
  ggsave(output_file, plot = p, width = 10, height = 6, dpi = 1500)
}


#####################################################
# plotting barplot of the number of independent clumps by approach and subtype
library(ggplot2)

# Sample data
clumps_df <- data.frame(
  num_clumps = c(21, 36, 5, 9, 1, 18, 14, 28),
  approach = c("plink2", "TOP", "plink2", "TOP", "plink2", "TOP", "plink2", "TOP"),
  subtype = c("HRNEG_HER2NEG", "HRNEG_HER2NEG", "HRPOS_HER2POS", "HRPOS_HER2POS", 
              "HRNEG_HER2POS", "HRNEG_HER2POS", "HRPOS_HER2NEG", "HRPOS_HER2NEG")
)

# Create the bar plot
plot <- ggplot(clumps_df, aes(x = subtype, y = num_clumps, fill = approach)) +
  geom_bar(stat = "identity", position = "dodge") +  # "dodge" ensures bars are side-by-side
  labs(
    x = "Subtype",
    y = "Number of independent loci"
  ) +
  scale_fill_manual(values = c("plink2" = "blue", "TOP" = "red")) +  # Custom colors
  theme_classic(base_size = 14) +  # Apply classic theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    plot.title = element_blank()  # Remove main title
  )

# Save the plot as a high-quality PNG
ggsave(
  filename = "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/number_clumps_approach.png",
  plot = plot,
  width = 8,   # Set width of the image
  height = 6,  # Set height of the image
  dpi = 300    # High resolution for better quality
)


##################################
# COMPARING EUR and AFR SUMSTATS #
##################################
subtype="HRNEG_HER2POS"
for (subtype in c("HRNEG_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRPOS_HER2NEG")) {

  # importing EUR sumstats 
  tmp_eur_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_sumstats/", subtype, ".tsv")) %>% separate(ID1,into=c("chr","pos","a2","a1"),sep="\\:",remove=F)
  tmp_eur_sumstats <- tmp_eur_sumstats %>% mutate(
    chr=as.numeric(chr),
    pos=as.numeric(pos)
  )
  tmp_eur_sumstats <- tmp_eur_sumstats %>% filter(!is.na(chr)) %>% filter(!is.na(pos))
  
  # importing AFR sumstats and making alleles consistent with IDs
  tmp_afr_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/TOP/METAANALYSIS_",subtype,"_1.tbl")) 
  tmp_afr_sumstats <- tmp_afr_sumstats %>% separate(MarkerName,into=c("chr","pos","a2","a1"),remove=F,sep="\\:")
  tmp_afr_sumstats <- tmp_afr_sumstats %>% mutate(
    chr=as.numeric(chr),
    pos=as.numeric(pos)
  )
 
  ##############################################
  # making uppercase alleles
  parsed_tmp_afr_sumstats <- tmp_afr_sumstats %>% 
    mutate(Allele1=toupper(Allele1),
           Allele2=toupper(Allele2))
  
  # altering direction and effect columns if A1 and Allele1 do not match
  parsed_tmp_afr_sumstats <- parsed_tmp_afr_sumstats %>%
    mutate(
      Direction = if_else(Allele1 != a1, chartr('+-', '-+', Direction), Direction),
      Effect=ifelse(Allele1==a1,Effect,-Effect),
      Allele1=a1,
      Allele2=a2
    )
  parsed_tmp_afr_sumstats <- parsed_tmp_afr_sumstats %>% mutate(ID=MarkerName,EffectAllele=a1,BaselineAllele=toupper(a2),BETA=Effect,SE=StdErr,P=`P-value`) %>% select(ID,EffectAllele,BaselineAllele,BETA,SE,P)
  
  # joining data
  join_ancestry_df <- inner_join(parsed_tmp_afr_sumstats,tmp_eur_sumstats%>%mutate(ID=ID1),by=c("ID","EffectAllele","BaselineAllele")) 
  join_ancestry_df <- join_ancestry_df %>% filter(P.x<5e-8) %>% mutate(BETA.x=as.numeric(BETA.x),BETA.y=as.numeric(BETA.y))
  
  #####################################
  # Perform linear regression
  regression_model <- lm(BETA.y ~ BETA.x, data = join_ancestry_df)
  
  # Extract regression statistics
  beta_coef <- coef(regression_model)["BETA.x"]
  r_squared <- summary(regression_model)$r.squared
  p_value <- summary(regression_model)$coefficients["BETA.x", "Pr(>|t|)"]
  
  # Create annotation label
  label <- bquote(beta == .(round(beta_coef, 2)) * 
                    "\n p: " * .(signif(p_value, 3)) * 
                    "\n " * r^2 * " = " * .(round(r_squared, 2)))
  
  # Generate scatter plot
  p <- ggplot(join_ancestry_df, aes(x = BETA.x, y = BETA.y)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "red", alpha = 0.5) +
    xlab(expression("African ancestry (" * beta * ")")) +
    ylab(expression("European ancestry (" * beta * ")")) +
    theme_classic() +
    theme(axis.text = element_text(size = 14),   # Increase axis text size
          axis.title = element_text(size = 16)) + # Increase axis title size
    annotate("text", x = min(join_ancestry_df$BETA.x, na.rm = TRUE),
             y = max(join_ancestry_df$BETA.y, na.rm = TRUE), 
             label = label, hjust = 0, vjust = 1, size = 5)
  
  # Define the file path
  file_path <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/scatterplot_compare_ancestry_effectsize_", subtype, ".png")
  
  # Save the plot
  ggsave(file_path, plot = p, width = 6, height = 6, dpi = 300)

}










head(tmp_eur_sumstats)
head(parsed_tmp_afr_sumstats)

tmp_eur_sumstats %>% filter(P<5e-8) %>% separate(ID, into=c("chr","pos","a2","a1"),sep="\\:") %>% select(chr) %>% unique() %>% arrange(chr)
tmp_afr_sumstats %>% filter(P<5e-8) %>% separate(ID, into=c("chr","pos","a2","a1"),sep="\\:") %>% select(chr) %>% unique() %>% arrange(chr)

joined_ancestry_sumstats_1 <- inner_join(parsed_tmp_afr_sumstats,tmp_eur_sumstats%>%mutate(ID=ID1),by=c("ID","EffectAllele","BaselineAllele"))
joined_ancestry_sumstats_2 <- inner_join(parsed_tmp_afr_sumstats,tmp_eur_sumstats%>%mutate(ID=ID2),by=c("ID","EffectAllele","BaselineAllele"))
dim(joined_ancestry_sumstats_1)
dim(joined_ancestry_sumstats_2)



dim(inner_join(tmp_afr_sumstats,tmp_eur_sumstats%>%mutate(ID=ID1),by=c("ID","EffectAllele","BaselineAllele")))

dim(inner_join(tmp_afr_sumstats,tmp_eur_sumstats%>%mutate(ID=ID2),by=c("ID","EffectAllele","BaselineAllele")))


