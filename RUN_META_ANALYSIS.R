library(data.table)
library(meta)

# Define the directory path
dir_path <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_plink2/"

# Define the input string for file filtering
input_string <- "HRNEG_HER2NEG"

# List all files in the directory matching the pattern
file_list <- list.files(path = dir_path, pattern = input_string, full.names = TRUE)

# Read all files, rename BETA and SE columns, and store them in a list
data_list <- lapply(file_list, function(file) {
  cohort_name <- sub(".*HRNEG_HER2NEG\\.", "", basename(file))
  cohort_name <- sub("\\.sumstats$", "", cohort_name)
  
  dt <- fread(file)
  dt[, P := NULL]  # Remove P column
  
  setnames(dt, c("BETA", "SE"), c(paste0("BETA.", cohort_name), paste0("SE.", cohort_name)))
  
  return(dt)
})

# Merge all data tables on ID, EffectAllele, and BaselineAllele
joined_sumstats <- Reduce(function(x, y) merge(x, y, by = c("ID", "EffectAllele", "BaselineAllele"), all = TRUE), data_list)

# Initialize empty columns
joined_sumstats[, `:=`(meta_beta = NA_real_, meta_se = NA_real_, meta_z = NA_real_, meta_p_value = NA_real_)]

# Get column names for BETA and SE
beta_cols <- grep("^BETA\\.", names(joined_sumstats), value = TRUE)
se_cols <- grep("^SE\\.", names(joined_sumstats), value = TRUE)

# Ensure SE column order matches BETA
se_suffixes <- sub("BETA\\.", "", beta_cols)
se_cols <- se_cols[order(match(sub("SE\\.", "", se_cols), se_suffixes))]

# Function to perform meta-analysis on a row
meta_analysis_function <- function(beta_values, se_values) {
  valid_idx <- !is.na(beta_values) & !is.na(se_values) & se_values > 0
  
  if (sum(valid_idx) > 1) {  # Need at least two valid studies
    meta_result <- tryCatch({
      metagen(TE = beta_values[valid_idx], 
              seTE = se_values[valid_idx], 
              sm = "OR", 
              common = TRUE, 
              random = FALSE)
    }, error = function(e) return(NULL))
    
    if (!is.null(meta_result)) {
      meta_beta <- meta_result$TE.common
      meta_se <- meta_result$seTE.common
      meta_z <- meta_beta / meta_se
      meta_p_value <- 2 * pnorm(-abs(meta_z))
      return(c(meta_beta, meta_se, meta_z, meta_p_value))
    }
  }
  return(rep(NA_real_, 4))
}

# Apply meta-analysis function across all rows (Vectorized approach)
result_matrix <- t(mapply(meta_analysis_function, 
                          as.data.frame(t(joined_sumstats[, ..beta_cols])),
                          as.data.frame(t(joined_sumstats[, ..se_cols]))))

# Assign results back to the data.table
joined_sumstats[, c("meta_beta", "meta_se", "meta_z", "meta_p_value") := as.data.table(result_matrix)]

# Save the result
output_file <- file.path(dir_path, "joined_sumstats_meta.tsv")
fwrite(joined_sumstats, output_file, sep = "\t")

cat("Meta-analysis completed. Output saved to:", output_file, "\n")
