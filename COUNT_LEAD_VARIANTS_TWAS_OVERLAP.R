library(topr)
library(dplyr)
library(data.table)
library(tidyr)

# set working directory to that containing all output lead variants
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry/")

# List all .tsv files in the current directory
all_files <- list.files(pattern = "\\.tsv$")

# Extract unique subtypes from filenames (e.g., HRPOS_HER2NEG)
subtypes <- unique(gsub("_(eur|afr|xancestry)\\.tsv$", "", all_files))

# Create a named list to store combined data.tables per subtype
combined_list <- list()

# Loop through each subtype and combine corresponding files
for (subtype in subtypes) {
  # Get matching files for the current subtype
  subtype_files <- list.files(pattern = paste0("^", subtype, "_(eur|afr|xancestry)\\.tsv$"))
  
  # Read and bind files with ancestry column
  combined_dt <- rbindlist(lapply(subtype_files, function(file) {
    dt <- fread(file)
    dt[, ancestry := sub(paste0("^", subtype, "_|\\.tsv$"), "", file)]
    return(dt)
  }), use.names = TRUE, fill = TRUE)
  
  # Store the result in the list
  combined_list[[subtype]] <- combined_dt
}

ALL_LEAD <- rbind(
  combined_list$HRPOS_HER2NEG,
  combined_list$HRPOS_HER2POS,
  combined_list$HRNEG_HER2POS,
  combined_list$HRNEG_HER2NEG
)


# across all meta-analyses
nrow(get_lead_snps(combined_list$HRPOS_HER2NEG%>%mutate(P = 2 * pnorm(-abs(BETA / SE)))%>%mutate(P=ifelse(P<1e-300,1e-300,P)),thresh = 1.25e-08,region_size = 4e+06))
     
nrow(get_lead_snps(combined_list$HRPOS_HER2POS%>%mutate(P = 2 * pnorm(-abs(BETA / SE)))%>%mutate(P=ifelse(P<1e-300,1e-300,P)),thresh = 1.25e-08,region_size = 4e+06))

nrow(get_lead_snps(combined_list$HRNEG_HER2POS%>%mutate(P = 2 * pnorm(-abs(BETA / SE)))%>%mutate(P=ifelse(P<1e-300,1e-300,P)),thresh = 1.25e-08,region_size = 4e+06))

nrow(get_lead_snps(combined_list$HRNEG_HER2NEG%>%mutate(P = 2 * pnorm(-abs(BETA / SE)))%>%mutate(P=ifelse(P<1e-300,1e-300,P)),thresh = 1.25e-08,region_size = 4e+06))

nrow(get_lead_snps(ALL_LEAD%>%mutate(P = 2 * pnorm(-abs(BETA / SE)))%>%mutate(P=ifelse(P<1e-300,1e-300,P)),thresh = 1.25e-08,region_size = 4e+06))

     
# each meta-analysis
system("for f in *; do [ -f \"$f\" ] && echo $(( $(wc -l < \"$f\") - 1 )) $f; done")

#################################################
# Assessing how many GWAS loci are in proximity with TWAS genes #
#################################################
library(data.table)
library(dplyr)
library(tidyr)

##########################################
# Step 1: Prepare pruned lead SNPs
##########################################
pruned_lead <- get_lead_snps(
  ALL_LEAD %>%
    mutate(P = 2 * pnorm(-abs(BETA / SE))) %>%
    mutate(P = ifelse(P < 1e-300, 1e-300, P)),
  thresh = 1.25e-08,
  region_size = 4e6
) %>%
  mutate(
    CHR = as.numeric(gsub("chr", "", CHROM)),
    POS = as.numeric(POS),
    MODIFY_POS = CHR * 1e14 + POS
  )

##########################################
# Step 2: Load gene lists
##########################################
# Previous TWAS genes
previous_intrinsic_twas <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/intrinsic_subtype_gene_list_long.tsv")
previous_intrinsic_twas_genes <- unique(previous_intrinsic_twas$ensg_id)

# Current TWAS results
annotated_twas <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO.tsv")
current_intrinsic_twas_genes <- unique(annotated_twas$ENSG)

##########################################
# Step 3: Load GENCODE and filter for relevant genes
##########################################
gencode26_all <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/gencode_v26_all.txt", sep = "\t") %>%
  rename(gene_symbol = gene_name) %>%
  filter(chromosome %in% paste0("chr", 1:22)) %>%
  mutate(
    CHR = as.numeric(gsub("chr", "", chromosome)),
    MODIFY_GENE_START = CHR * 1e14 + start_location,
    MODIFY_GENE_END = CHR * 1e14 + end_location
  ) %>%
  separate(gene_id, into = c("ENSG", "decimal"), sep = "\\.")

# Filter for previous and current TWAS genes separately
gencode_prev_genes <- gencode26_all %>%
  filter(ENSG %in% previous_intrinsic_twas_genes)

gencode_curr_genes <- gencode26_all %>%
  filter(ENSG %in% current_intrinsic_twas_genes)

##########################################
# Step 4: Identify nearby current TWAS genes (±2Mb)
##########################################
setDT(pruned_lead)
pruned_lead[, nearby_curr_twas := FALSE]

for (i in 1:nrow(gencode_curr_genes)) {
  start_pos <- gencode_curr_genes$MODIFY_GENE_START[i] - 2e6
  end_pos   <- gencode_curr_genes$MODIFY_GENE_END[i] + 2e6
  pruned_lead[MODIFY_POS >= start_pos & MODIFY_POS <= end_pos, nearby_curr_twas := TRUE]
}

##########################################
# Step 5: Identify nearby previous TWAS genes (±2Mb)
##########################################
pruned_lead[, nearby_prev_twas := FALSE]

for (i in 1:nrow(gencode_prev_genes)) {
  start_pos <- gencode_prev_genes$MODIFY_GENE_START[i] - 2e6
  end_pos   <- gencode_prev_genes$MODIFY_GENE_END[i] + 2e6
  pruned_lead[MODIFY_POS >= start_pos & MODIFY_POS <= end_pos, nearby_prev_twas := TRUE]
}

##########################################
# Step 6: Combined summary for current + previous
##########################################
pruned_lead[, nearby_any_twas := nearby_curr_twas | nearby_prev_twas]

##########################################
# Step 7: Summary
##########################################
cat("Nearby current TWAS genes:\n")
print(table(pruned_lead$nearby_curr_twas))

cat("Nearby previous TWAS genes:\n")
print(table(pruned_lead$nearby_prev_twas))

cat("Nearby any TWAS gene (current or previous):\n")
print(table(pruned_lead$nearby_any_twas))

##########################################
# Step 8: Preparing for phenotype plotting
##########################################
# parsing the input table for phenogram visualization
phenogram_input <- pruned_lead %>% mutate(chrom=gsub("chr","",CHROM),pos=POS) %>% select(chrom,pos,nearby_prev_twas,nearby_curr_twas)

# apply conditions
phenogram_input[, phenotype := fifelse(
  nearby_prev_twas == FALSE & nearby_curr_twas == FALSE, 
  "GWAS loci without nearby TWAS signals",
  fifelse(nearby_prev_twas == FALSE & nearby_curr_twas == TRUE,
          "GWAS loci with nearby signals from current TWAS",
          fifelse(nearby_prev_twas == TRUE & nearby_curr_twas == FALSE,
                  "GWAS loci with nearby signals from previous TWAS",
                  "GWAS loci with nearby signals in current and previous TWAS")))]

# writing out table
write.table(phenogram_input%>%select(-nearby_prev_twas,-nearby_curr_twas),file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/relation_gwas/input.tsv",quote=F,row.names=F,col.names=T,sep="\t")

# writing out color conversion table
color_conversion <- data.frame(
  phenotype = c(
    "GWAS loci with nearby signals in current and previous TWAS",
    "GWAS loci with nearby signals from previous TWAS",
    "GWAS loci without nearby TWAS signals",
    "GWAS loci with nearby signals from current TWAS"
  ),
  color = c("purple", "red", "white", "blue"),
  stringsAsFactors = FALSE
)
write.table(color_conversion%>%select(color,phenotype),file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/relation_gwas/color_conversion.tsv",quote=F,row.names=F,col.names=T,sep="\t")
