
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO.tsv")
TABLE1<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/TABLE1.tsv")
TABLE2<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/TABLE2.tsv")

library(dplyr)
library(tidyr)

# all genes
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>%
  group_by(subtype) %>%
  summarise(n_genes = n_distinct(ENSG))
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>%
  group_by(subtype) %>%
  summarise(n_loci = n_distinct(loci))

ensg_counts_table_subtype_celltype <- annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>%
  group_by(subtype, celltype) %>%
  summarise(n_unique_ENSG = n_distinct(ENSG), .groups = "drop") %>%
  pivot_wider(names_from = celltype, values_from = n_unique_ENSG, values_fill = 0) %>%
  as.data.frame()
ensg_counts_table_subtype_celltype

ensg_counts_by_celltype <- annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>%
  group_by(celltype) %>%
  summarise(n_unique_ENSG = n_distinct(ENSG), .groups = "drop") %>%
  as.data.frame()
print(ensg_counts_by_celltype)

ensg_counts_by_ancestry <- annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>%
  group_by(ancestry) %>%
  summarise(n_unique_ENSG = n_distinct(ENSG), .groups = "drop") %>%
  as.data.frame()
print(ensg_counts_by_ancestry)


# novel genes
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% filter(ENSG%in%TABLE1$ENSG) %>%
  group_by(subtype) %>%
  summarise(n_genes = n_distinct(ENSG))
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% filter(ENSG%in%TABLE1$ENSG) %>%
  group_by(subtype) %>%
  summarise(n_loci = n_distinct(loci))
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% filter(ENSG%in%TABLE1$ENSG) %>% select(ENSG) %>% unique() %>% nrow()
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% filter(ENSG%in%TABLE1$ENSG) %>% select(loci) %>% unique() %>% nrow()


# previously reported genes
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% filter(!(ENSG%in%TABLE1$ENSG)) %>%
  group_by(subtype) %>%
  summarise(n_genes = n_distinct(ENSG))
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% filter(!(ENSG%in%TABLE1$ENSG)) %>%
  group_by(subtype) %>%
  summarise(n_loci = n_distinct(loci))
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% filter(!(ENSG%in%TABLE1$ENSG)) %>% select(ENSG) %>% unique() %>% nrow()
annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO %>% filter(!(ENSG%in%TABLE1$ENSG)) %>% select(loci) %>% unique() %>% nrow()


# subtype specific genes







