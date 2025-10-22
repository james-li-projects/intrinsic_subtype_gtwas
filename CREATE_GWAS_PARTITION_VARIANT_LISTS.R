###########################
# importing libraries
library(dplyr)
library(data.table)

###########################
# importing pvar file
pvar <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS.pvar") 
pvar <- pvar %>% select(ID)

###########################
# creating partitions
set.seed(123) # Set seed for reproducibility

# Directory to save the partition files
output_dir <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/partition/variant_lists/"

# Ensure the output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Number of rows per partition
rows_per_partition <- 5000

# Total number of rows in pvar
n_rows <- nrow(pvar)

# Create a sequence of partition indices
partition_indices <- split(seq_len(n_rows), ceiling(seq_len(n_rows) / rows_per_partition))

# Loop through each partition and save as a file
for (i in seq_along(partition_indices)) {
  # Subset the dataframe
  partition <- pvar[partition_indices[[i]], , drop = FALSE]
  
  # Create a file name with zero-padded partition index
  file_name <- sprintf("partition_%05d.txt", i)
  
  # Write the partition to a file
  write.table(
    partition,
    file = file.path(output_dir, file_name),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}
