library(data.table)
library(dplyr)
library(readxl)
library(tidyr)

# processing genes in a recently published TWAS of BC intrinsic subtypes [PMID: 38723630, 38697998] 
gencode26_raw<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/gencode_v26_all.txt",sep="\t") %>% rename(gene_symbol=gene_name) %>% dplyr::select(gene_symbol,gene_id)
pmid_38723630_38697998 <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/unconverted_pmid_38723630_38697998_genes.txt")
joined_pmid_38723630_38697998 <- left_join(pmid_38723630_38697998,gencode26_raw) %>% filter(!is.na(gene_id)) %>% separate(gene_id,into=c("ensg_id","decimal"),sep="\\.") %>% dplyr::select(-decimal)

# subsetting recently published TWAS by whether the genes were reported for overall BC or subtypes
overall_pmid_38723630_38697998 <- joined_pmid_38723630_38697998 %>% filter(subtype=="overall") %>% select(gene_symbol,ensg_id,pmid)
subtype_pmid_38723630_38697998 <- joined_pmid_38723630_38697998 %>% filter(subtype=="subtype") %>% select(gene_symbol,ensg_id,pmid)

# importing other previous TWAS gene lists
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/")
older_twas_genes <- fread("older_twas_genes.txt")
recent_twas_genes <- fread("recent_twas_genes.txt")
# adding in the above imported genes from PMID: 38723630_38697998
recent_twas_genes <- rbind(
  recent_twas_genes,
  overall_pmid_38723630_38697998,
  subtype_pmid_38723630_38697998
)

# compiling an intrinsic subtype TWAS gene list 
subtype_pmid_38400758 <- recent_twas_genes %>% filter(pmid=="38400758")
intrinsic_subtype_gene_list <- rbind(
  subtype_pmid_38400758,
  subtype_pmid_38723630_38697998
) %>% unique()

# assembling combined twas gene lists in both long and wide formats
setDT(older_twas_genes)
setDT(recent_twas_genes)

recent_collapsed <- recent_twas_genes[, .(pmid_list = paste(unique(pmid), collapse = ", ")), 
                                      by = .(gene_symbol, ensg_id)]

combined_gene_list_wide <- rbind(older_twas_genes, recent_collapsed)[
  , .(pmid_list = paste(sort(unique(unlist(strsplit(pmid_list, ",\\s*")))), collapse = ", ")),
  by = .(gene_symbol, ensg_id)
]

older_expanded <- older_twas_genes[
  , .(pmid = as.integer(unlist(strsplit(pmid_list, ",\\s*")))), 
  by = .(gene_symbol, ensg_id)
]

combined_gene_list_long <- unique(rbind(older_expanded, recent_twas_genes))

# writing out gene list tsv files
fwrite(combined_gene_list_wide, "combined_gene_list_wide.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
fwrite(combined_gene_list_long, "combined_gene_list_long.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
fwrite(intrinsic_subtype_gene_list, "intrinsic_subtype_gene_list_long.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


# making a finalized, aggregated wide table
library(data.table)
# Define a function to clean gene_symbol vector
remove_ensg_ids <- function(symbols) {
  clean <- unique(symbols[!grepl("^ENSG\\d{11}$", symbols)])
  if (length(clean) == 0) clean <- unique(symbols)  # fallback if only ENSG-style names
  paste(clean, collapse = ", ")
}
# Aggregate
combined_gene_list_wide_final <- combined_gene_list_wide[
  , .(
    gene_symbol = remove_ensg_ids(gene_symbol),
    pmid_list = paste(unique(unlist(strsplit(pmid_list, ",\\s*"))), collapse = ", ")
  ),
  by = ensg_id
]
# annotate whether identified previously in intrinsic subtype TWAS
combined_gene_list_wide_final <- combined_gene_list_wide_final %>%
  mutate(previous_intrinsic_twas = if_else(ensg_id %in% intrinsic_subtype_gene_list$ensg_id, "Yes", "No"))
# writing out this finalized wide list
fwrite(combined_gene_list_wide_final, "combined_gene_list_wide_final.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

