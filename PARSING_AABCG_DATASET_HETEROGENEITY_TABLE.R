library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

# Set the directory path
dir_path <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/parsed_output/"

# List all matching files
file_list <- list.files(path = dir_path, pattern = "_FINAL_LEAD_VARIANT_DF\\.tsv$", full.names = TRUE)

# Read and combine all files
combined_dt <- rbindlist(lapply(file_list, fread), use.names = TRUE, fill = TRUE) %>% mutate(Alleles=paste(EffectAllele,BaselineAllele,sep="/"))

# selecting columns
selected_cols <- c(
  "SUBTYPE", "LOCI", "rsid", "ID", "Alleles","AFR_Direction",
  "OR_AABC", "P_AABC", "AF_AABC",
  "OR_AMBER", "P_AMBER", "AF_AMBER",
  "OR_BCAC_OncoArray", "P_BCAC_OncoArray", "AF_BCAC_OncoArray",
  "OR_GBHS", "P_GBHS", "AF_GBHS",
  "OR_MEGA", "P_MEGA", "AF_MEGA",
  "OR_ROOT", "P_ROOT", "AF_ROOT",
  "OR_WGS", "P_WGS", "AF_WGS"
)
filtered_dt <- combined_dt[, ..selected_cols]
filtered_dt[, names(filtered_dt) := lapply(.SD, as.character)]

# specifying columns to keep as identifiers for wide to long conversion
id_cols <- c("SUBTYPE", "LOCI", "rsid", "ID", "Alleles","AFR_Direction")

# Melt OR, P, and AF columns separately to long format
OR_filtered_dt <- filtered_dt[, !grepl("^(P_|AF_)", names(filtered_dt)), with = FALSE]
OR_long <- melt(OR_filtered_dt, measure.vars = patterns("^OR_"), variable.name = "Platform", value.name = "OR") %>% mutate(Platform=gsub("AF_|P_|OR_", "",Platform))
P_filtered_dt <- filtered_dt[, !grepl("^(OR_|AF_)", names(filtered_dt)), with = FALSE]
P_long <- melt(P_filtered_dt, measure.vars = patterns("^P_"), variable.name = "Platform", value.name = "P") %>% mutate(Platform=gsub("AF_|P_|OR_", "",Platform))
AF_filtered_dt <- filtered_dt[, !grepl("^(P_|OR_)", names(filtered_dt)), with = FALSE]
AF_long <- melt(AF_filtered_dt, measure.vars = patterns("^AF_"), variable.name = "Platform", value.name = "AF") %>% mutate(Platform=gsub("AF_|P_|OR_", "",Platform))

# joining these tables
final_long_df<-inner_join(OR_long,P_long) %>% inner_join(AF_long) 

# processing final platform heterogeneity table
final_long_df <- final_long_df %>%
  mutate(SUBTYPE = trimws(SUBTYPE)) %>%
  mutate(SUBTYPE = factor(SUBTYPE, levels = c("HRPOS_HER2NEG", "HRPOS_HER2POS", "HRNEG_HER2NEG")))
final_long_df <- final_long_df %>% separate(ID,into=c("chr","pos","ref","alt"),sep="\\:",remove=F) %>% select(-ref,-alt) %>% mutate(chr=as.numeric(chr),pos=as.numeric(pos)) %>% arrange(SUBTYPE,chr,pos,Platform) %>% select(-ID,-chr,-pos)

# writing out table
write.table(final_long_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/platform_heterogeneity/platform_het.tsv",quote=F,row.names=F,col.names=T,sep="\t")
















































"AFR_Direction"