# Load required libraries
library(data.table)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)

########################################
# importing list of lead variants
HRPOS_HER2POS_afr_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRPOS_HER2POS.txt") %>% filter(Ancestry=="afr"))$ID
HRPOS_HER2NEG_afr_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRPOS_HER2NEG.txt") %>% filter(Ancestry=="afr"))$ID
HRNEG_HER2NEG_afr_lead_variant_list <- (fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/HRNEG_HER2NEG.txt") %>% filter(Ancestry=="afr"))$ID
afr_lead_variant_list <- c(
  HRPOS_HER2POS_afr_lead_variant_list,
  HRPOS_HER2NEG_afr_lead_variant_list,
  HRNEG_HER2NEG_afr_lead_variant_list
)

########################################
# Set directory path
dir_path <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats"
# List all TSV files
tsv_files <- list.files(dir_path, pattern = "\\.tsv$", full.names = TRUE)
# Initialize an empty list to store data frames
filtered_dfs <- list()
# Loop through files
for (file in tsv_files) {
  prefix <- str_remove(basename(file), "\\.tsv$")
  df <- read_tsv(file, show_col_types = FALSE) %>%
    filter(ID %in% afr_lead_variant_list) %>%
    select(ID, CHR, POS, EffectAllele, BaselineAllele, BETA, SE) %>%
    rename(!!paste0(prefix, "_BETA") := BETA,
           !!paste0(prefix, "_SE") := SE)
  filtered_dfs[[prefix]] <- df
}
# Full join all data frames
joined_df <- data.frame(reduce(filtered_dfs, full_join, by = c("ID", "CHR", "POS", "EffectAllele", "BaselineAllele")))
# Convert BETA and SE columns to formatted OR (95% CI: L-U)
formatted_df <- data.frame(joined_df %>%
  pivot_longer(
    cols = matches("_BETA$|_SE$"),
    names_to = c("Prefix", ".value"),
    names_pattern = "(.*)_(BETA|SE)"
  ) %>%
  mutate(
    OR = round(exp(BETA), 2),
    CI_lower = round(exp(BETA - 1.96 * SE), 2),
    CI_upper = round(exp(BETA + 1.96 * SE), 2),
    Value = paste0(OR, " (", CI_lower, "-", CI_upper, ")")
  ) %>%
  select(ID, CHR, POS, EffectAllele, BaselineAllele, Prefix, Value) %>%
  pivot_wider(names_from = Prefix, values_from = Value)
)

########################################
# importing heterogeneity test p-values 
het_p <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/subtype_heterogeneity/het_p_afr.tsv") %>% rename(ID=variant,het_p=p_value)

########################################
# assembling final df with ORs and heterogeneity p-values
final_df <- full_join(formatted_df,het_p,by=c("ID"))

# write output 
write.table(final_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/subtype_heterogeneity/heterogeneity_table_afr.tsv",quote=F,row.names=F,col.names=T,sep="\t")
