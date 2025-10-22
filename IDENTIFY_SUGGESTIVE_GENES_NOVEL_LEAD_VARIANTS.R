library(data.table)
library(tidyr)

novel_variants <- c(
  "1:220613960:G:A",
  "10:5707486:G:A",
  "10:9058859:T:G",
  "18:32351771:A:G",
  "20:59573644:T:G",
  "10:123814971:T:C",
  "12:96591322:G:A",
  "2:173354390:C:G",
  "1:9156852:C:T",
  "1:31087245:A:G",
  "11:133440533:A:C",
  "2:19121207:C:T",
  "5:6120548:T:C"
)

long_twas <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/long_twas.tsv")

############################################
# Identify genes nearby novel GWAS regions #
############################################
integrate_gwas_twas_hits <- data.frame()
# initializing subtype
for (tmp_subtype in c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2NEG")) {
  # importing lead variants for this subtype and creating a modified position for evaluating distance to TWAS hits
  tmp_lead_variant_df <- data.frame(ID=setdiff(
    fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/",tmp_subtype,".txt"))$ID,
    c("3:197512284:C:T","13:110604015:T:C","6:20534230:G:GA")
  )) %>% 
    separate(ID,into=c("CHR","POS","REF","ALT"),remove=F,sep="\\:") %>% 
    mutate(CHR=as.numeric(CHR),POS=as.numeric(POS)) %>%
    mutate(MODIFY_POS=CHR*1e14 + POS) %>% 
    select(ID,MODIFY_POS)
  
  # modifying df of bonferroni-significant TWAS hits to compare with lead variant positions 
  
  modified_bonf_significant_twas <- long_twas %>% filter(pvalue<1e-3) %>% separate(ucsc_cytoband,into=c("CHR","band"),remove=F,sep="\\:") %>% unique() %>% mutate(CHR=as.numeric(CHR)) %>% mutate(MODIFY_GENE_START = CHR*1e14 + gene_start,MODIFY_GENE_END = CHR*1e14 + gene_end) %>% filter(subtype==tmp_subtype)
  
  # modified_bonf_significant_twas <- bonf_significant_twas %>% separate(ucsc_cytoband,into=c("CHR","band"),remove=F,sep="\\:") %>% unique() %>% mutate(CHR=as.numeric(CHR)) %>% mutate(MODIFY_GENE_START = CHR*1e14 + gene_start,MODIFY_GENE_END = CHR*1e14 + gene_end) %>% filter(subtype==tmp_subtype)
  
  # Identify genes near novel lead variants
  tmp_integrate_gwas_twas_hits <- data.frame(
    modified_bonf_significant_twas %>%
      rowwise() %>%
      mutate(variant_id = paste(
        tmp_lead_variant_df %>%
          filter(
            between(MODIFY_POS, MODIFY_GENE_START - 2e6, MODIFY_GENE_START + 2e6) |
              between(MODIFY_POS, MODIFY_GENE_END - 2e6, MODIFY_GENE_END + 2e6)
          ) %>%
          pull(ID),
        collapse = ","
      )) %>%
      ungroup() %>%
      filter(variant_id != "")
  )
  
  integrate_gwas_twas_hits <- rbind(
    integrate_gwas_twas_hits,
    tmp_integrate_gwas_twas_hits
  )
}

###########################################
# creating a summary DF for these gene associations
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# Mappings
celltype_map <- c(
  "Adipocytes" = "Adip", "Endothelial_cells" = "Endo", 
  "Epithelial_cells" = "Epi", "Stromal_and_Immune_cells" = "SIC", 
  "Breast_tissue" = "Bulk", "Fibroblasts" = "Fib"
)

ancestry_map <- c(
  "WHITE" = "E", "BLACK" = "A", "XANCESTRY" = "C"
)

# Assume `integrate_gwas_twas_hits` is already loaded
summary_df <- integrate_gwas_twas_hits %>%
  group_by(gene_name, ensg_id, gene_type, ucsc_cytoband, CHR, band) %>%
  summarize(
    min_pvalue = min(pvalue, na.rm = TRUE),
    min_zscore = zscore[which.min(pvalue)],
    
    summary_stat = sprintf("%.2f (%.2E)", min_zscore, min_pvalue),
    
    subtypes = unique(subtype) %>% sort() %>% paste(collapse = ", "),
    CT = celltype_map[unique(celltype)] %>% unique() %>% sort() %>% paste(collapse = ", "),
    ancestries = ancestry_map[unique(ancestry)] %>% unique() %>% sort() %>% paste(collapse = ", ")
  ) %>%
  ungroup()

# Final selection and arrangement
summary_df <- data.frame(summary_df %>%
  select(-min_pvalue, -min_zscore) %>%
  arrange(CHR, ucsc_cytoband, gene_name)) %>% rename(subtype=subtypes)

# View summary
head(summary_df)




###########################################
# specifying target gene/variant hits
target_hits <- integrate_gwas_twas_hits %>% filter(variant_id %in% novel_variants)


########################################
# extracting prediction model variants #
########################################
library(RSQLite)
library(dplyr)
library(stringr)

# Set working directory
setwd("/gpfs/data/huo-lab/Vanderbilt/Julian/03_make_models/output/predixcan_dbs/")

# List all gene model database files
db_paths <- list.files(pattern = "gene_models.*\\.db$")

# Remove version suffix for safe variable naming
target_hits <- target_hits %>%
  mutate(ensg_id_clean = str_replace(ensg_id, "\\.\\d+$", ""))

# Unique gene x ancestry x celltype combos
unique_genes <- target_hits %>%
  select(ensg_id, ensg_id_clean, ancestry, celltype) %>%
  distinct()

# Storage for merged gene weights
gene_weight_list <- list()

# Loop through all combinations
for (i in seq_len(nrow(unique_genes))) {
  ensg_id     <- unique_genes$ensg_id[i]
  gene_var    <- unique_genes$ensg_id_clean[i]
  ancestry    <- unique_genes$ancestry[i]
  celltype    <- unique_genes$celltype[i]
  
  # Determine which ancestries to query
  ancestry_pool <- if (ancestry == "XANCESTRY") c("WHITE", "BLACK") else ancestry
  
  # Try each ancestry in the pool
  for (anc in ancestry_pool) {
    db_file_pattern <- paste0("^", anc, "\\.gene_models\\.", celltype, "\\.db$")
    db_file <- db_paths[grepl(db_file_pattern, db_paths)]
    
    if (length(db_file) == 1) {
      con <- dbConnect(SQLite(), db_file)
      weights_df <- dbReadTable(con, "weights")
      dbDisconnect(con)
      
      gene_weights <- weights_df %>%
        filter(gene == ensg_id) %>%
        select(varID, weight) %>%
        mutate(ancestry = anc, celltype = celltype)
      
      # Combine entries
      gene_weight_list[[gene_var]] <- bind_rows(
        gene_weight_list[[gene_var]], gene_weights
      )
      
      message(sprintf("Appended %s (%s) from %s", gene_var, ensg_id, db_file))
    } else if (length(db_file) == 0) {
      warning(sprintf("No DB found for ancestry=%s, celltype=%s", anc, celltype))
    } else {
      warning(sprintf("Multiple DBs matched for ancestry=%s, celltype=%s", anc, celltype))
    }
  }
}

# Optional: unpack into individual variables
list2env(gene_weight_list, envir = .GlobalEnv)

# specifying genes to analyze
analyzed_targets <- target_hits %>% select(subtype,gene_name,ensg_id,variant_id) %>% unique()
analyzed_targets$maxR2 <- NA 

# computing correlations using plink between lead variant and model variants
for (i in 1:nrow(analyzed_targets)) {
  
  model_variants <- (gene_weight_list[[i]])%>%select(varID)%>%unique()
  
  target_lead_variant <- data.frame(varID=c(analyzed_targets$variant_id[i]))
  
  containslead=intersect(model_variants$varID,target_lead_variant$varID)
  
  if (length(containslead)==1) {
    analyzed_targets$maxR2[i] <- 1
  } else {
    total_variant_list <- rbind(model_variants,target_lead_variant)
    
    write.table(total_variant_list,file="/scratch/jll1/tmp/variant.list",quote=F,row.names=F,col.names=F)
    
    max_snp_window=nrow((gene_weight_list[[i]])%>%select(varID)%>%unique())+1
    
    system(paste0("module load plink/1.9; plink -bfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/afr --extract /scratch/jll1/tmp/variant.list --ld-snp ",analyzed_targets$variant_id[i]," --ld-window-kb 4000 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --memory 100000 --out /scratch/jll1/tmp/tmp_LD"))
    pairwise_cor_1 <- fread("/scratch/jll1/tmp/tmp_LD.ld")
    
    system(paste0("module load plink/1.9; plink -bfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/eur --extract /scratch/jll1/tmp/variant.list --ld-snp ",analyzed_targets$variant_id[i]," --ld-window-kb 4000 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --memory 100000 --out /scratch/jll1/tmp/tmp_LD"))
    pairwise_cor_2 <- fread("/scratch/jll1/tmp/tmp_LD.ld")
    pairwise_cor <- rbind(pairwise_cor_1,pairwise_cor_2)
    
    pairwise_cor <- pairwise_cor %>%
      filter(!(SNP_A == analyzed_targets$variant_id[i] & SNP_B == analyzed_targets$variant_id[i]))
    analyzed_targets$maxR2[i]<-max(pairwise_cor$R2)
    print(max(pairwise_cor$R2))
  }
}

# examining max correlations for nearby genes identified at 1e-3
print(analyzed_targets)

# joining with summary DF
joined_df <- inner_join(summary_df,analyzed_targets,by=c("subtype","gene_name","ensg_id")) %>% select(subtype,variant_id,gene_name,ensg_id,summary_stat,CT,ancestries,maxR2) %>% filter(maxR2>0.1)
write.table(joined_df,"/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/suggestive_genes/suggestive_genes.tsv",quote=F,row.names=F,col.names=T,sep="\t")



#####################################
#####################################
#####################################
# creating a summary DF for NEAREST gene associations that are not necessarily relevant for these GWAS lead variants
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# Mappings
celltype_map <- c(
  "Adipocytes" = "Adip", "Endothelial_cells" = "Endo", 
  "Epithelial_cells" = "Epi", "Stromal_and_Immune_cells" = "SIC", 
  "Breast_tissue" = "Bulk", "Fibroblasts" = "Fib"
)

ancestry_map <- c(
  "WHITE" = "E", "BLACK" = "A", "XANCESTRY" = "C"
)

# importing nearest genes
nearest_gene_df <- data.frame(
  subtype = c(
    "HRPOS_HER2NEG", "HRPOS_HER2NEG", "HRPOS_HER2NEG",
    "HRPOS_HER2POS", "HRPOS_HER2POS", "HRPOS_HER2POS",
    "HRNEG_HER2NEG", "HRNEG_HER2NEG", "HRNEG_HER2NEG"
  ),
  gene_name = c(
    "MARK1", "FAM208B", "GAREM1",
    "CPXM2", "CFAP54", "AC092573.2",
    "MIR34AHG", "OPCML", "LINC01376"
  ),
  stringsAsFactors = FALSE
)

nearest_gene_summary_df <- long_twas %>% inner_join(nearest_gene_df,by=c("subtype","gene_name")) %>% separate(ucsc_cytoband,into = c("CHR","band"),sep="\\:",remove=F) %>%
  group_by(gene_name, ensg_id, gene_type, ucsc_cytoband, CHR, band) %>%
  summarize(
    min_pvalue = min(pvalue, na.rm = TRUE),
    min_zscore = zscore[which.min(pvalue)],
    
    summary_stat = sprintf("%.2f (%.2E)", min_zscore, min_pvalue),
    
    subtypes = unique(subtype) %>% sort() %>% paste(collapse = ", "),
    CT = celltype_map[unique(celltype)] %>% unique() %>% sort() %>% paste(collapse = ", "),
    ancestries = ancestry_map[unique(ancestry)] %>% unique() %>% sort() %>% paste(collapse = ", ")
  ) %>%
  ungroup()

# Final selection and arrangement
nearest_gene_summary_df <- data.frame(nearest_gene_summary_df %>%
                                        select(-min_pvalue, -min_zscore) %>%
                                        arrange(CHR, ucsc_cytoband, gene_name)) %>% rename(subtype=subtypes)
