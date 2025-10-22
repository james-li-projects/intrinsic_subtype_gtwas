library(dplyr)
library(data.table)
library(stringr)

# Set the directory where your files are located
input_dir <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry"

# Define ancestries
ancestries <- c("afr", "eur", "xancestry")

# Initialize list to store data.tables by ancestry
all_by_ancestry <- list()

for (anc in ancestries) {
  # List files for the current ancestry
  files <- list.files(path = input_dir, pattern = paste0("_", anc, "\\.tsv$"), full.names = TRUE)
  
  # Read and annotate each file
  dt_list <- lapply(files, function(file) {
    # Extract the subtype from the file name
    subtype <- str_extract(basename(file), "^[^_]+_[^_]+")
    
    # Read the file and add the subtype column
    dt <- fread(file)
    dt[, subtype := subtype]
    return(dt)
  })
  
  # Combine all data.tables for this ancestry
  all_by_ancestry[[anc]] <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
}

# You now have:
# all_by_ancestry$afr
# all_by_ancestry$eur
# all_by_ancestry$xancestry

# defining list of problematic variants to remove (identified by performing analysis and subsequently removing)
problem_variant_list <- c("3:197512284:C:T","13:110604015:T:C","6:20534230:G:GA")

# import ID to RSID conversion dictionary
ID_RSID_dict <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/rsid_dictionary/ID_RSID_dict.tsv")










library(data.table)
library(stringr)

# Define directories
input_dir <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry"
output_dir <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry_parsed"

# Ancestries
ancestries <- c("afr", "eur", "xancestry")

# Desired subtype sorting order
subtype_levels <- c("HRPOS_HER2NEG", "HRPOS_HER2POS", "HRNEG_HER2POS", "HRNEG_HER2NEG")

# Ensure rsid mapping is available
ID_RSID_dict <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/rsid_dictionary/ID_RSID_dict.tsv")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Process each ancestry
for (anc in ancestries) {
  files <- list.files(input_dir, pattern = paste0("_", anc, "\\.tsv$"), full.names = TRUE)
  
  dt_list <- lapply(files, function(file) {
    subtype <- str_extract(basename(file), "^[^_]+_[^_]+")
    dt <- fread(file)
    dt[, subtype := subtype]
    return(dt)
  })
  
  dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  
  # Merge rsid
  dt <- merge(dt, ID_RSID_dict, by = "ID", all.x = TRUE)
  
  # ALT/REF 
  dt[, Alleles := paste(ALT, REF, sep = "/")]
  
  # Calculate OR and CI
  dt[, `:=`(
    OR = exp(BETA),
    OR_lower = exp(BETA - 1.96 * SE),
    OR_upper = exp(BETA + 1.96 * SE),
    P = 2 * pnorm(-abs(BETA / SE))
  )]
  
  # Format OR + 95% CI
  dt[, OR_95CI := sprintf("%.2f (%.2f-%.2f)", OR, OR_lower, OR_upper)]
  
  # Numeric CHROM for sorting
  dt[, CHROM_num := as.numeric(gsub("chr", "", CHROM))]
  
  # Factorize subtype with desired order
  dt[, subtype := factor(subtype, levels = subtype_levels)]
  
  # Sort and select columns
  dt <- dt[order(subtype, CHROM_num, POS)]
  dt <- dt[, .(subtype, LOCI, rsid, ID, Alleles, OR_95CI, P)]
  
  # remove chr string from loci
  dt <- dt %>% mutate(LOCI=gsub("chr","",LOCI))
  
  # Write to file
  fwrite(
    dt,
    file = file.path(output_dir, paste0("parsed_", anc, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}
