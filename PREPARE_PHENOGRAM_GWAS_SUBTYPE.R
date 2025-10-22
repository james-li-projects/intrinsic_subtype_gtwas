library(data.table)
library(dplyr)

# Set directory
input_dir <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry"
output_file <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/phenogram/input_subtype.tsv"

# List all .tsv files
files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

# Function to read and process each file
process_file <- function(file_path) {
  dt <- fread(file_path)
  
  # Add subtype based on filename
  if (grepl("HRPOS_HER2NEG", file_path)) {
    dt[, phenotype := "Luminal-A-like"]
  } else if (grepl("HRPOS_HER2POS", file_path)) {
    dt[, phenotype := "Luminal-B-like"]
  } else if (grepl("HRNEG_HER2POS", file_path)) {
    dt[, phenotype := "HER2-enriched-like"]
  } else if (grepl("HRNEG_HER2NEG", file_path)) {
    dt[, phenotype := "Triple-negative"]
  } else {
    dt[, phenotype := NA_character_]
  }
  
  # Remove files that don't match any subtype
  dt <- dt[!is.na(phenotype)]
  
  # Rename POS to pos
  setnames(dt, "POS", "pos")
  
  # Extract chr number from CHROM
  dt[, chr := sub("chr", "", CHROM)]
  
  # Retain only required columns
  dt[, .(LOCI, chr, pos, phenotype)]
}

# Apply processing to all files and combine
combined_dt <- rbindlist(lapply(files, process_file), use.names = TRUE, fill = TRUE)

# Retain only unique records
combined_dt <- combined_dt %>% 
  mutate(annotation = gsub("chr", "", LOCI)) %>% 
  select(annotation, chr, pos, phenotype) %>% 
  unique()

# Make sure your data is a data.table
setDT(combined_dt)

# Optional: sort for reproducibility
setorder(combined_dt, annotation, phenotype, pos)

# Step 1: Randomly thin records within (annotation, phenotype)
thin_within_group <- function(dt, max_distance = 2e6) {
  keep <- integer()
  remaining <- dt
  
  while (nrow(remaining) > 0) {
    selected_index <- sample(nrow(remaining), 1)
    selected_pos <- remaining[selected_index, pos]
    keep <- c(keep, remaining[selected_index, .I])
    remaining <- remaining[abs(pos - selected_pos) > max_distance]
  }
  
  dt[keep]
}

filtered_dt <- combined_dt[, thin_within_group(.SD), by = .(annotation, phenotype)]

# Step 2: Unify pos for all records with the same annotation
adjusted_dt <- copy(filtered_dt)

adjusted_dt[, {
  pos_new <- min(pos)
  adjusted_dt[.I, pos := pos_new]
  NULL
}, by = annotation]

# Final result
result <- adjusted_dt %>% unique()

# Write to output file
fwrite(result, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# writing out color conversion table
color_conversion <- data.frame(
  phenotype = c(
    "Luminal-A-like",
    "Luminal-B-like",
    "HER2-enriched-like",
    "Triple-negative"
  ),
  color = c("royalblue", "lightgreen", "tan1", "violetred"),
  stringsAsFactors = FALSE
)
write.table(color_conversion%>%select(color,phenotype),file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/phenogram/input_subtype_color_conversion.tsv",quote=F,row.names=F,col.names=T,sep="\t")
