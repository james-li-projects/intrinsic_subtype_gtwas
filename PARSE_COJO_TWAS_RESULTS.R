library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)

####################################
# Importing bonferroni significant TWAS results
bonf_significant_twas <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/bonf_significant_twas.tsv")
bonf_significant_twas <- bonf_significant_twas %>% separate(ucsc_cytoband,into=c("CHR","band"),sep="\\:",remove=F)

####################################
# Importing COJO TWAS results
raw_cojo_twas <- fread("/gpfs/data/huo-lab/Vanderbilt/Julian/040_parse_through_initial_twas/output/spredixcan_combined/cojo/james_vars_cojo_condTWAS_summary.tsv") %>% filter(celltype!="Stromal_and_Immune_cells") %>%
  unique() %>%
  rename(zscore=cond_zscore,pvalue=cond_pvalue) %>%
  mutate(category = paste(subtype, celltype, ancestry, sep = ".")) %>%
  select(-pvalue, -ancestry, -celltype, -subtype) %>%
  arrange(category)

####################################
# Importing mean effective sample sizes
load("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/mean_n_eff/n_eff_df.RData")

# Clean up subtype and race
n_eff_df$subtype <- trimws(as.character(n_eff_df$subtype))
n_eff_df$race <- trimws(as.character(n_eff_df$race))

####################################
# Make wide format
wide_cojo_twas <- data.frame(raw_cojo_twas %>%
                          pivot_wider(
                            names_from = category,
                            values_from = zscore
                          ))

####################################
# Get all column names that end in .BLACK or .WHITE
hr_cols <- grep("^HR.*\\.(BLACK|WHITE)$", names(wide_cojo_twas), value = TRUE)

# Extract unique prefixes
prefixes <- unique(sub("\\.(BLACK|WHITE)$", "", hr_cols))

# Loop over prefixes
for (prefix in prefixes) {
  black_col <- paste0(prefix, ".BLACK")
  white_col <- paste0(prefix, ".WHITE")
  
  if (black_col %in% names(wide_cojo_twas) && white_col %in% names(wide_cojo_twas)) {
    df_pair <- data.frame(
      wide_cojo_twas[[black_col]],
      wide_cojo_twas[[white_col]]
    )
    
    # Extract subtype string from prefix
    subtype_str <- trimws(strsplit(prefix, "\\.")[[1]][1])
    
    # Correct subtype filtering
    n_black <- n_eff_df %>%
      filter(race == "BLACK", subtype == subtype_str) %>%
      pull(n_eff)
    
    n_white <- n_eff_df %>%
      filter(race == "WHITE", subtype == subtype_str) %>%
      pull(n_eff)
    
    if (length(n_black) != 1 | length(n_white) != 1) {
      warning(paste("Check subtype:", subtype_str, "- unexpected number of matches in n_eff_df"))
      wide_cojo_twas[[paste0(prefix, ".XANCESTRY")]] <- rep(NA, nrow(wide_cojo_twas))
    } else {
      n_black <- n_black[1]
      n_white <- n_white[1]
      
      meta_z <- apply(df_pair, 1, function(z) {
        z1 <- z[1]
        z2 <- z[2]
        
        if (is.na(z1) && is.na(z2)) return(NA)
        if (is.na(z1)) return(z2)
        if (is.na(z2)) return(z1)
        
        numerator <- sqrt(n_black) * z1 + sqrt(n_white) * z2
        denominator <- sqrt(n_black + n_white)
        return(numerator / denominator)
      })
      
      new_col_name <- paste0(prefix, ".XANCESTRY")
      wide_cojo_twas[[new_col_name]] <- meta_z
    }
  }
}


####################################
# converting COJO TWAS results back to long format
long_cojo_twas <- data.frame(wide_cojo_twas %>%
                          pivot_longer(
                            cols = starts_with("HR"),  # all your z-score columns
                            names_to = "category",     # name for the category column
                            values_to = "zscore"       # name for the zscore column
                          )) %>% filter(!is.na(zscore)) %>% separate(category,into=c("subtype","celltype","ancestry"),sep="\\.") %>% mutate(pvalue=2*pnorm(-abs(zscore)))
# renaming pvalue columns
long_cojo_twas <- long_cojo_twas %>% rename(cond_pvalue=pvalue,cond_zscore=zscore)


####################################
# obtaining parsed version of target COJO results
joined_COJO_results <- inner_join(bonf_significant_twas,long_cojo_twas%>%filter(status=="successful"),by=c("ensg_id","subtype","celltype","ancestry")) %>% separate(ensg_id,into=c("ENSG","decimal"),sep="\\.") %>% select(-decimal) %>% select(ENSG,subtype,celltype,ancestry,num_gene_test,cond_zscore,cond_pvalue)

# Importing bonferroni significant TWAS results that have been annotated with loci and subtype
annotated_bonf_significant_twas_LOCI_SUBTYPE <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI_SUBTYPE.tsv") 
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO <- left_join(annotated_bonf_significant_twas_LOCI_SUBTYPE,joined_COJO_results,by=c("ENSG","subtype","celltype","ancestry")) %>% mutate(COJO_sig=ifelse(cond_pvalue < 0.0125 / num_gene_test / 4,"Yes","No"))
# Adding a modified genomic position values
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO <- annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% mutate(MODIFY_GENE_START = numeric_chr*1e14 + gene_start,MODIFY_GENE_END = numeric_chr*1e14 + gene_end)

# Importing lead variants for each subtype
base_path <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_lead_variants"
ancestry_map <- c(eur = "WHITE", afr = "BLACK")
all_lead_variant_df <- list.files(base_path, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE) |>
  lapply(function(file) {
    ancestry <- ancestry_map[basename(dirname(file))]
    subtype <- tools::file_path_sans_ext(basename(file))
    read.delim(file)[, c("ID", "CHR", "POS")] |>
      mutate(subtype = subtype, ancestry = ancestry)
  }) |>
  bind_rows()
# Vector of problematic AFR variants
exclude_ids <- c("3:197512284:C:T", "13:110604015:T:C", "6:20534230:G:GA")
# Excluding these problematic variants
all_lead_variant_df <- all_lead_variant_df %>%
  filter(!(ID %in% exclude_ids & ancestry == "BLACK"))
# Adding a modified genomic position value
all_lead_variant_df <- all_lead_variant_df %>% mutate(MODIFY_POS = CHR*1e14 + POS)


# Identifying closest lead variant to each gene as well as the distance
library(data.table)

setDT(annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO)
setDT(all_lead_variant_df)

# Initialize result columns
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO[, `:=`(
  closest_lead_variant = NA_character_,
  distance_closest_lead_variant = NA_real_
)]

# Iterate through each row and compute distance to nearest index variant
for (i in seq_len(nrow(annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO))) {
  row <- annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO[i]
  
  ancestry_filter <- if (row$ancestry == "XANCESTRY") c("BLACK", "WHITE") else row$ancestry
  lv_sub <- all_lead_variant_df[subtype == row$subtype & ancestry %in% ancestry_filter]
  
  if (nrow(lv_sub) > 0) {
    # Check if MODIFY_POS is within the gene region — assign distance = 0 if so
    distances <- ifelse(
      lv_sub$MODIFY_POS >= row$MODIFY_GENE_START & lv_sub$MODIFY_POS <= row$MODIFY_GENE_END,
      0,
      pmin(
        abs(lv_sub$MODIFY_POS - row$MODIFY_GENE_START),
        abs(lv_sub$MODIFY_POS - row$MODIFY_GENE_END)
      )
    )
    min_idx <- which.min(distances)
    
    annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO[i, `:=`(
      closest_lead_variant = lv_sub$ID[min_idx],
      distance_closest_lead_variant = distances[min_idx]
    )]
  }
}

# writing out this fully annotated table
write.table(annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO.tsv",quote=F,row.names=F,col.names=T,sep="\t")

###################################
# identifying genes >2Mb away from GWAS lead variants
not_explained_by_GWAS_genes_1 <- annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>%
  filter(is.na(COJO_sig))

# identifying COJO significant genes within 2Mb of GWAS lead variants
not_explained_by_GWAS_genes_2 <- annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>%
  filter(COJO_sig=="Yes")

# compiling DF of TWAS signals not explained by GWAS variants
not_explained_by_GWAS_genes <- rbind(
  not_explained_by_GWAS_genes_2,
  not_explained_by_GWAS_genes_1
) %>% select(-gene_end,-subtype_specific,-num_gene_test,-MODIFY_GENE_START,-MODIFY_GENE_END)

############################
# creating an organized table of variants not explained by GWAS signals
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# Mappings
celltype_map <- c(
  "Adipocytes" = "Adip", "Endothelial_cells" = "Endo", 
  "Epithelial_cells" = "Epi", "Stromal_and_Immune_cells" = "SIC", 
  "Breast_tissue" = "Bulk", "Fibroblasts" = "Fib"
)

ancestry_map <- c(
  "WHITE" = "E", "BLACK" = "A", "XANCESTRY" = "C"
)

# Assuming your DF is named not_explained_by_GWAS_genes
summary_df <- not_explained_by_GWAS_genes %>%
  group_by(subtype, loci, gene_name, ENSG, gene_type, numeric_chr,gene_start,COJO_sig) %>%
  summarize(
    min_pvalue = min(pvalue, na.rm = TRUE),
    min_zscore = zscore[which.min(pvalue)],
    summary_stat = sprintf("%.2f (%.2E)", min_zscore, min_pvalue),
    
    # For conditional summary
    cond_min_pvalue = if (all(is.na(COJO_sig))) NA_real_ else min(cond_pvalue, na.rm = TRUE),
    cond_min_zscore = if (all(is.na(COJO_sig))) NA_real_ else cond_zscore[which.min(cond_pvalue)],
    cond_summary_stat = if (all(is.na(COJO_sig))) NA_character_ else sprintf("%.2f (%.2E)", cond_min_zscore, cond_min_pvalue),
    
    CT = celltype_map[unique(celltype)] %>% unique() %>% sort() %>% paste(collapse = ", "),
    ancestries = ancestry_map[unique(ancestry)] %>% unique() %>% sort() %>% paste(collapse = ", "),
    
    closest_lead_variant = first(closest_lead_variant),
    distance_closest_lead_variant = first(distance_closest_lead_variant)
  ) %>%
  ungroup()
summary_df <- data.frame(summary_df) %>%
  mutate(distance_closest_lead_variant_kb = round(distance_closest_lead_variant / 1000, 2)) %>% arrange(COJO_sig,desc(subtype),numeric_chr,gene_start) %>%
  select(-distance_closest_lead_variant,-numeric_chr,-gene_start,-COJO_sig,-min_pvalue,-min_zscore,-cond_min_pvalue,-cond_min_zscore)

# writing output table
write.table(summary_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/TABLE2.tsv",quote=F,row.names=F,col.names=T,sep="\t")

