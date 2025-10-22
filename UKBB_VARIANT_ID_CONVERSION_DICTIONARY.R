library(data.table)
library(dplyr)
library(tidyr)

chr_bim <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/xancestry_ld_ref/eur/CrossMap_chr",i))

# extracting IDs from current meta-analysis results
dirs <- c(
  "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats",
  "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_meta_sumstats"
)
# Initialize an empty list to store ID columns
id_list <- list()
# Loop through directories and read ID columns from all files
for (dir in dirs) {
  files <- list.files(dir, full.names = TRUE)
  for (file in files) {
    dt <- fread(file, select = "ID")
    id_list[[length(id_list) + 1]] <- dt
  }
}
# Combine all IDs into one data.table
all_ids <- rbindlist(id_list, use.names = FALSE)
# Get unique IDs
unique_ids <- unique(all_ids)
# obtaining different orientations

unique_ids <- 

