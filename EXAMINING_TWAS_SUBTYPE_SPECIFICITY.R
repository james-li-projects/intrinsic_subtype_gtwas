library(data.table)
library(dplyr)
library(tidyr)
library(venn)

# importing previously identified TWAS genes for BC intrinsic subtypes
previous_TWAS <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/intrinsic_subtype_gene_list_long.tsv")
previous_TWAS_genes <- unique(previous_TWAS$ensg_id)

# importing genes identified in this study
genes_to_report_df <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/bonf_significant_twas.tsv") %>% separate(ensg_id,sep="\\.",into=c("ENSG_nodecimal","decimal")) %>%
  mutate(subtype = recode(subtype,"HRPOS_HER2NEG" = 1,"HRPOS_HER2POS" = 2,"HRNEG_HER2POS" = 3,"HRNEG_HER2NEG" = 4)) %>% select(ENSG_nodecimal,subtype) %>% unique()

########################################
# initializing list to store genes identified for each subtype in our TWAS
subtype_gene_list <- vector(mode = "list", length = 5)

# add the subtype genes to the subtype gene list
for (i in 1:4) {
  subtype_gene_list[[i]] <- (genes_to_report_df %>% filter(subtype == i))$ENSG_nodecimal
}

subtype_gene_list[[5]]<- previous_TWAS_genes

# labeling these lists
names(subtype_gene_list) <- c("Luminal-A-like", "Luminal-B-like", "HER2-enriched-like", "TNBC", "Previously reported")
subtype_gene_list_4 <- subtype_gene_list[1:4]
names(subtype_gene_list_4) <- c("Luminal-A-like", "Luminal-B-like", "HER2-enriched-like", "TNBC")


########################################
# plotting venn diagram of subtypes and previous genes
library(VennDiagram)
library(grid)

png("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/subtype_specificity/venn_genes_subtypes_previousTWAS.png", width = 1200, height = 1200, units = "px", pointsize = 25)
grid.draw(venn.diagram(
  x = subtype_gene_list,
  category.names = names(subtype_gene_list),
  filename = NULL,
  fill = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6"),
  col = "black", lwd = 2, alpha = 0.7,
  cat.col = "black", cat.cex = 1, cat.fontfamily = "sans",
  cex = 1, fontfamily = "sans",
  margin = 0.3
))
dev.off()

########################################
# plotting venn diagram of subtypes only from this study
library(VennDiagram)
library(grid)

# Use only the first 4 subtypes
png("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/subtype_specificity/venn_genes_subtypes_this_study.png", width = 1200, height = 1200, units = "px", pointsize = 25)

grid.draw(venn.diagram(
  x = subtype_gene_list_4,
  category.names = names(subtype_gene_list_4),
  filename = NULL,
  fill = c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4"),
  col = "black", lwd = 2, alpha = 0.7,
  cat.col = "black", cat.cex = 1, cat.fontfamily = "sans",
  cex = 1, fontfamily = "sans",
  margin = 0.3
))
dev.off()


########################################
# plotting venn diagram of subtypes and previous genes
library(eulerr)
library(tidyverse)
# saving/plotting output
png("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/subtype_specificity/venn_genes_subtypes_previousTWAS.png",width = 12, height = 12, units = "in", res = 1500)
venn_plot<-euler(subtype_gene_list,counts=T)
plot(
  venn_plot,
  fills = c("lightpink", "lightgreen", "skyblue", "lavender", "lightyellow"),
  labels = list(labels = names(subtype_gene_list), font = 2, cex = 0),
  quantities = list(cex = 2, font = 1))
dev.off()

########################################
# plotting venn diagram of subtypes only from this study
library(eulerr)
library(tidyverse)
# saving/plotting output
png("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/subtype_specificity/venn_genes_subtypes_this_study.png",width = 12, height = 12, units = "in", res = 1500)
venn_plot<-euler(subtype_gene_list_4,counts=T)
plot(
  venn_plot,
  fills = c("lightpink", "lightgreen", "skyblue", "lavender"),
  labels = list(labels = names(subtype_gene_list_4), font = 2, cex = 0),
  quantities = list(cex = 2, font = 1))
dev.off()


########################################
# examining subtype specificity, excluding genes with a stringent P value threshold (P<0.05) in other subtypes
long_twas <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/long_twas.tsv")
bonf_significant_twas <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/bonf_significant_twas.tsv")

# Ensure data.table format
setDT(long_twas)
setDT(bonf_significant_twas)

# Lowest p-value per ensg_id + subtype, keeping only relevant columns
lowest_p_long <- long_twas[
  , .SD[which.min(pvalue), .(pvalue)], 
  by = .(ensg_id, subtype)
]

lowest_p_bonf <- bonf_significant_twas[
  , .SD[which.min(pvalue), .(pvalue)], 
  by = .(ensg_id, subtype)
]

# Ensure data.table format
setDT(lowest_p_bonf)
setDT(lowest_p_long)

# 1. Identify ensg_ids that are subtype-specific in lowest_p_bonf
bonf_counts <- lowest_p_bonf[, .N, by = ensg_id]
subtype_specific_ids <- bonf_counts[N == 1, ensg_id]  # Keep only those appearing once

# 2. Get their associated subtype
subtype_specific_genes <- lowest_p_bonf[ensg_id %in% subtype_specific_ids, .(ensg_id, subtype)]

# 3. Initialize result list
results <- list()

# 4. For each gene, verify that in lowest_p_long it has p > 0.05 for all *other* subtypes
for (i in seq_len(nrow(subtype_specific_genes))) {
  gene <- subtype_specific_genes[i, ensg_id]
  this_subtype <- subtype_specific_genes[i, subtype]
  
  gene_pvalues <- lowest_p_long[ensg_id == gene]
  
  # Subtypes other than the one it's "specific" to
  other_subtypes <- gene_pvalues[subtype != this_subtype]
  
  # Check condition: all p-values in other subtypes > 0.05
  if (nrow(other_subtypes) == 0 || all(other_subtypes$pvalue > 0.05, na.rm = TRUE)) {
    results[[length(results) + 1]] <- data.table(
      ensg_id = gene,
      subtype = this_subtype
    )
  }
}

# 5. Combine into one data.table
subtype_specific_df <- rbindlist(results) %>% separate(ensg_id,into=c("ENSG_nodecimal","decimal"),sep="\\.")

# identifying all and novel genes that are subtype specific
subtype_specific_gene_list <- subtype_specific_df$ENSG_nodecimal
novel_gene_list <- setdiff(genes_to_report_df$ENSG_nodecimal,previous_TWAS_genes)
novel_subtype_specific_gene_list <- intersect(novel_gene_list,subtype_specific_gene_list)


########################################
# annotating all bonferroni significant hits by their subtype specific status
annotated_bonf_significant_twas_LOCI_SUBTYPE <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI.tsv") %>% mutate(subtype_specific=ifelse(ENSG%in%subtype_specific_gene_list,"Yes","No"))

# writing out these dually annotated significant TWAS hits
write.table(annotated_bonf_significant_twas_LOCI_SUBTYPE,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI_SUBTYPE.tsv",quote=F,row.names=F,col.names=T,sep="\t")


########################################
# compiling organized table
library(dplyr)
library(tidyr)
library(stringr)

# Mappings
celltype_map <- c(
  "Adipocytes" = "Adip", "Endothelial_cells" = "Endo", 
  "Epithelial_cells" = "Epi", "Stromal_and_Immune_cells" = "SIC", 
  "Breast_tissue" = "Bulk", "Fibroblasts" = "Fib"
)

ancestry_map <- c(
  "WHITE" = "E", "BLACK" = "A", "XANCESTRY" = "C"
)

# Processing the annotated bonferroni significant TWAS hits
annotated_bonf_significant_twas <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI.tsv")
summary_df <- annotated_bonf_significant_twas %>%
  group_by(subtype,loci,numeric_chr,gene_start,gene_end,gene_name,ENSG,gene_type) %>%
  summarize(
    min_pvalue = min(pvalue, na.rm = TRUE),
    min_zscore = zscore[which.min(pvalue)],
    
    summary_stat = sprintf("%.2f (%.2E)", min_zscore, min_pvalue),
    
    CT = celltype_map[unique(celltype)] %>% unique() %>% sort() %>% paste(collapse = ", "),
    ancestries = ancestry_map[unique(ancestry)] %>% unique() %>% sort() %>% paste(collapse = ", ")
  ) %>%
  ungroup() 
summary_df <- data.frame(summary_df) %>% select(-min_pvalue,-min_zscore) %>% arrange(desc(subtype),numeric_chr,gene_start) %>% select(-numeric_chr,-gene_start,-gene_end)

# filtering this summary DF to only novel genes
TABLE1 <- summary_df %>% filter(ENSG%in%novel_gene_list) %>% mutate(subtype_specific=ifelse(ENSG%in%novel_subtype_specific_gene_list,"Yes","No"))
write.table(TABLE1,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/TABLE1.tsv",quote=F,row.names=F,col.names=T,sep="\t")

