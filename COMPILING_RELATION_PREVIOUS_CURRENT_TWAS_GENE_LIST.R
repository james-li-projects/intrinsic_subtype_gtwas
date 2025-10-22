library(data.table)
library(dplyr)

# importing bonferroni-significant TWAS hits
## all
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO.tsv")
## parsed
TABLE1 <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/TABLE1.tsv")

# importing previously identified TWAS genes
combined_gene_list_wide_final<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/combined_gene_list_wide_final.tsv")

# creating a column indicating whether a previously reported gene was identified in the current TWAS
combined_gene_list_wide_final <- combined_gene_list_wide_final %>%
  mutate(current_intrinsic_twas = if_else(ensg_id %in% annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO$ENSG, "Yes", "No"))

# Previous TWAS hits newly identified as intrinsic subtype genes
new_intrinsic_genes <- (combined_gene_list_wide_final %>% filter(previous_intrinsic_twas=="No",current_intrinsic_twas=="Yes"))$ensg_id
TABLE1 %>% filter(ENSG %in% new_intrinsic_genes) %>% filter(subtype_specific=="Yes")

# writing out results
write.table(combined_gene_list_wide_final,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/relation_previous_twas/compiled_previous_current_gene_list.tsv",quote=F,row.names=F,col.names=T,sep="\t")
