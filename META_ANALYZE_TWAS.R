library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)

####################################
# Importing TWAS results
raw_twas <- fread("/gpfs/data/huo-lab/Vanderbilt/Julian/040_parse_through_initial_twas/output/spredixcan_combined/ancestry.subtype.celltype.TWAS.results.txt") %>% filter(celltype!="Stromal_and_Immune_cells") %>%
  unique() %>%
  mutate(category = paste(subtype, celltype, ancestry, sep = ".")) %>%
  select(-pvalue, -ancestry, -celltype, -subtype) %>%
  arrange(category)

####################################
# Importing mean effective sample sizes
load("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/mean_n_eff/n_eff_df.RData")

# Clean up subtype and race
n_eff_df$subtype <- trimws(as.character(n_eff_df$subtype))
n_eff_df$race <- trimws(as.character(n_eff_df$race))

####################################
# Make wide format
wide_twas <- data.frame(raw_twas %>%
                          pivot_wider(
                            names_from = category,
                            values_from = zscore
                          ))

####################################
# Get all column names that end in .BLACK or .WHITE
hr_cols <- grep("^HR.*\\.(BLACK|WHITE)$", names(wide_twas), value = TRUE)

# Extract unique prefixes
prefixes <- unique(sub("\\.(BLACK|WHITE)$", "", hr_cols))

# Loop over prefixes
for (prefix in prefixes) {
  black_col <- paste0(prefix, ".BLACK")
  white_col <- paste0(prefix, ".WHITE")
  
  if (black_col %in% names(wide_twas) && white_col %in% names(wide_twas)) {
    df_pair <- data.frame(
      wide_twas[[black_col]],
      wide_twas[[white_col]]
    )
    
    # Extract subtype string from prefix
    subtype_str <- trimws(strsplit(prefix, "\\.")[[1]][1])
    
    # Correct subtype filtering
    n_black <- n_eff_df %>%
      filter(race == "BLACK", subtype == subtype_str) %>%
      pull(n_eff)
    
    n_white <- n_eff_df %>%
      filter(race == "WHITE", subtype == subtype_str) %>%
      pull(n_eff)
    
    if (length(n_black) != 1 | length(n_white) != 1) {
      warning(paste("Check subtype:", subtype_str, "- unexpected number of matches in n_eff_df"))
      wide_twas[[paste0(prefix, ".XANCESTRY")]] <- rep(NA, nrow(wide_twas))
    } else {
      n_black <- n_black[1]
      n_white <- n_white[1]
      
      meta_z <- apply(df_pair, 1, function(z) {
        z1 <- z[1]
        z2 <- z[2]
        
        # Only perform meta-analysis if both are non-NA
        if (is.na(z1) || is.na(z2)) return(NA)
        
        numerator <- sqrt(n_black) * z1 + sqrt(n_white) * z2
        denominator <- sqrt(n_black + n_white)
        return(numerator / denominator)
      })
      
      new_col_name <- paste0(prefix, ".XANCESTRY")
      wide_twas[[new_col_name]] <- meta_z
    }
  }
}


####################################
# converting TWAS results back to long format
long_twas <- data.frame(wide_twas %>%
  pivot_longer(
    cols = starts_with("HR"),  # all your z-score columns
    names_to = "category",     # name for the category column
    values_to = "zscore"       # name for the zscore column
  )) %>% filter(!is.na(zscore)) %>% separate(category,into=c("subtype","celltype","ancestry"),sep="\\.") %>% mutate(pvalue=2*pnorm(-abs(zscore)))

# filtering out genes that haven't been experimentally confirmed
long_twas <- long_twas %>% filter(gene_type!="TEC")

# writing out all TWAS results
write.table(long_twas,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/long_twas.tsv",quote=F,row.names=F,col.names=T,sep="\t")


####################################
# identify Bonferroni-significant TWAS hits, accounting for number of genes, subtypes, and cell types
num_gene_test_XANCESTRY <- length(unique((long_twas)$ensg_id))
num_gene_test_BLACK <- length(unique((long_twas%>%filter(ancestry=="BLACK"))$ensg_id))
num_gene_test_WHITE <- length(unique((long_twas%>%filter(ancestry=="WHITE"))$ensg_id))
bonf_significant_twas <- long_twas %>%
  mutate(num_gene_test = case_when(
    ancestry == "BLACK" ~ num_gene_test_BLACK,
    ancestry == "WHITE" ~ num_gene_test_WHITE,
    TRUE ~ num_gene_test_XANCESTRY)) %>% filter(pvalue < 0.0125 / num_gene_test / 4)

# removing genes within 2 Mb of problematic lead variants 
bonf_significant_twas <- bonf_significant_twas %>% separate(ucsc_cytoband,into=c("CHR","band"),remove=F,sep="\\:") %>% unique() %>% mutate(CHR=as.numeric(CHR)) %>% mutate(MODIFY_GENE_START = CHR*1e14 + gene_start,MODIFY_GENE_END = CHR*1e14 + gene_end)

# identifying problematic lead variants and parsing them
problematic_lead_variant_df <- data.frame(ID=c(
  "3:197512284:C:T",
  "13:110604015:T:C",
  "6:20534230:G:GA"
)) %>% 
  separate(ID,into=c("CHR","POS","REF","ALT"),remove=F,sep="\\:") %>% 
  mutate(CHR=as.numeric(CHR),POS=as.numeric(POS)) %>%
  mutate(MODIFY_POS=CHR*1e14 + POS) %>% 
  select(ID,MODIFY_POS)

# performing the removal based on a 2Mb window - no genes are removed
threshold <- 2000000
remove_idx <- apply(bonf_significant_twas[, c("MODIFY_GENE_START", "MODIFY_GENE_END")], 1, function(pos_pair) {
  any(abs(pos_pair[1] - problematic_lead_variant_df$MODIFY_POS) <= threshold) |
    any(abs(pos_pair[2] - problematic_lead_variant_df$MODIFY_POS) <= threshold)
})
bonf_significant_twas_filtered <- bonf_significant_twas[!remove_idx, ]
bonf_significant_twas <- bonf_significant_twas_filtered %>% select(-CHR,-band,-MODIFY_GENE_START,-MODIFY_GENE_END)

# writing out these significant TWAS hits
write.table(bonf_significant_twas,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/bonf_significant_twas.tsv",quote=F,row.names=F,col.names=T,sep="\t")

####################################
# identify genes uniquely identified in each analysis
XANCESTRY_twas<-bonf_significant_twas %>% filter(ancestry=="XANCESTRY")
WHITE_twas<-bonf_significant_twas %>% filter(ancestry=="WHITE")
BLACK_twas<-bonf_significant_twas %>% filter(ancestry=="BLACK")
# unique xancestry TWAS genes
unique_xancestry_ensg <- setdiff(
  unique(XANCESTRY_twas$ensg_id),
  unique(WHITE_twas$ensg_id)
)
XANCESTRY_twas %>% filter(ensg_id %in% unique_xancestry_ensg)
# unique black TWAS genes
unique_BLACK_ensg <- setdiff(
  unique(BLACK_twas$ensg_id),
  unique(WHITE_twas$ensg_id)
)
BLACK_twas %>% filter(ensg_id %in% unique_BLACK_ensg)

####################################
# filter genes where the lowest p-value for the gene is from XANCESTRY
lowest_pval_xancestry <- data.frame(bonf_significant_twas %>%
  group_by(ensg_id) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  slice(1) %>%
  filter(ancestry == "XANCESTRY")
)


################################################
#### Identifying genes to perform COJO TWAS ####
################################################
subtype_str_list = c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRNEG_HER2NEG")
race_str_list = c("WHITE","BLACK")
ancestry_str_list = c("eur","afr")

for (i in 1:length(subtype_str_list)) {
  for (j in 1:length(race_str_list)) {
    # initializing current subtype, race, and ancestry to examine
    tmp_subtype_str = subtype_str_list[i]
    tmp_race_str = race_str_list[j]
    tmp_ancestry_str = ancestry_str_list[j]
    
    print(paste("Processing COJO TWAS input for",tmp_subtype_str,tmp_race_str,tmp_ancestry_str))
    
    # generate modified positions to harmonize the TWAS and lead variant results to scan for potential genes to run COJO analysis on
    tmp_twas_coordinates <- bonf_significant_twas %>% filter(ancestry%in%c(tmp_race_str,"XANCESTRY")) %>% filter(subtype==tmp_subtype_str) %>% separate(ucsc_cytoband,into=c("CHR","band"),remove=F,sep="\\:") %>% select(ensg_id,CHR,gene_start,gene_end) %>% unique() %>% mutate(CHR=as.numeric(CHR)) %>% mutate(MODIFY_GENE_START = CHR*1e14 + gene_start,MODIFY_GENE_END = CHR*1e14 + gene_end) %>% select(-gene_start,-gene_end,-CHR)
    # importing lead variants for the specified ancestry
    current_lead_variants <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_lead_variants/",tmp_ancestry_str,"/",tmp_subtype_str,".tsv")) %>% mutate(MODIFY_POS=CHR*1e14 + POS) %>% select(ID,MODIFY_POS)
    
    # Ensure both are data.tables
    setDT(tmp_twas_coordinates)
    setDT(current_lead_variants)
    
    # Precompute the +/- 2Mb window
    tmp_twas_coordinates[, gene_start_window := MODIFY_GENE_START - 2e6]
    tmp_twas_coordinates[, gene_end_window := MODIFY_GENE_END + 2e6]
    
    # Now do the non-equi join using those new columns
    tmp_target_twas_coordinates <- current_lead_variants[
      tmp_twas_coordinates,
      on = .(MODIFY_POS >= gene_start_window,
             MODIFY_POS <= gene_end_window),
      nomatch = 0,
      allow.cartesian = TRUE
    ]
    
    # Keep unique matched gene coordinates (i.e. within +/- 2Mb window)
    tmp_target_cojo_twas_genes <- (unique(tmp_target_twas_coordinates[, .(ensg_id, MODIFY_GENE_START, MODIFY_GENE_END)]))$ensg_id
    write.table(data.frame(tmp_target_cojo_twas_genes),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_twas_genes/gene_lists/",tmp_subtype_str,".",tmp_race_str,".target_cojo_twas_genes.txt"),quote=F,row.names=F,col.names=F)
    
    # Get genes outside the 2Mb window
    tmp_notarget_cojo_twas_genes <- (tmp_twas_coordinates[!ensg_id %in% tmp_target_twas_coordinates$ensg_id])$ensg_id
    write.table(data.frame(tmp_notarget_cojo_twas_genes),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_twas_genes/gene_lists/",tmp_subtype_str,".",tmp_race_str,".notarget_cojo_twas_genes.txt"),quote=F,row.names=F,col.names=F)
    
    # extract all TWAS gene hits for a given subtype
    tmp_cojo_input_table <- long_twas %>% filter(subtype==tmp_subtype_str) %>% filter(ensg_id %in% tmp_target_cojo_twas_genes) %>% filter(ancestry==tmp_race_str)
    write.table(tmp_cojo_input_table,paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_twas_genes/parsed_twas_output/",tmp_subtype_str,".",tmp_race_str,".cojo_input_table.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
  }
}

############################################
# Identify genes nearby novel GWAS regions #
############################################
integrate_gwas_twas_hits <- data.frame()

for (tmp_subtype in c("HRPOS_HER2NEG", "HRPOS_HER2POS", "HRNEG_HER2NEG")) {
  
  # Import and process lead variants
  tmp_lead_variant_df <- data.frame(ID = setdiff(
    fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/", tmp_subtype, ".txt"))$ID,
    c("3:197512284:C:T", "13:110604015:T:C", "6:20534230:G:GA")
  )) %>%
    separate(ID, into = c("CHR", "POS", "REF", "ALT"), sep = ":", remove = FALSE) %>%
    mutate(CHR = as.numeric(CHR), POS = as.numeric(POS)) %>%
    mutate(MODIFY_POS = CHR * 1e14 + POS) %>%
    select(ID, MODIFY_POS)
  
  # Process Bonferroni-significant TWAS hits
  modified_bonf_significant_twas <- bonf_significant_twas %>%
    separate(ucsc_cytoband, into = c("CHR", "band"), sep = ":", remove = FALSE) %>%
    mutate(CHR = as.numeric(CHR)) %>%
    mutate(
      MODIFY_GENE_START = CHR * 1e14 + gene_start,
      MODIFY_GENE_END = CHR * 1e14 + gene_end
    ) %>%
    filter(subtype == tmp_subtype)
  
  # Match TWAS genes to nearby lead variants
  tmp_integrate_gwas_twas_hits <- modified_bonf_significant_twas %>%
    rowwise() %>%
    mutate(
      variant_id = paste(
        tmp_lead_variant_df %>%
          filter(between(MODIFY_POS, MODIFY_GENE_START - 2e6, MODIFY_GENE_END + 2e6)) %>%
          pull(ID),
        collapse = ","
      )
    ) %>%
    ungroup() %>%
    filter(variant_id != "")
  
  integrate_gwas_twas_hits <- bind_rows(integrate_gwas_twas_hits, tmp_integrate_gwas_twas_hits)
}

######################################
# Annotating cytobands for TWAS hits #
######################################
# import cytoband dictionary 
cytoband_dictionary <- fread("/gpfs/data/huo-lab/splice_TWAS/gencode/ucsc.hg38.cytoband.txt",header=T)

# annotate cytoband based on whether the genes are contained within a single cytoband or span multiple cytobands
annotated_bonf_significant_twas <- bonf_significant_twas %>% separate(ucsc_cytoband,into=c("chr","band"),sep="\\:") %>% select(-band) %>% mutate(numeric_chr=as.numeric(chr))
annotated_bonf_significant_twas$loci <- NA
for (i in 1:nrow(annotated_bonf_significant_twas)) {
  chr_str <- paste0("chr",annotated_bonf_significant_twas$numeric_chr[i])
  start_index <- annotated_bonf_significant_twas$gene_start[i] 
  end_index <- annotated_bonf_significant_twas$gene_end[i] 
  
  band1 <- tail((cytoband_dictionary %>% filter(Chr == chr_str,Start < start_index))$Band,1)
  chrband1 <- paste0(annotated_bonf_significant_twas$numeric_chr[i],band1)
  
  band2 <- head((cytoband_dictionary %>% filter(Chr == chr_str, End > end_index))$Band,1)
  chrband2 <- paste0(annotated_bonf_significant_twas$numeric_chr[i],band2)
  
  band_str<-NA
  
  if (chrband1==chrband2) {
    band_str <- chrband1
  } else {
    band_str <- paste(chrband1,band2,sep="-")
  }
  annotated_bonf_significant_twas$loci[i] <- band_str
}

# rearranging columns
annotated_bonf_significant_twas <- annotated_bonf_significant_twas %>% arrange(desc(subtype),numeric_chr,gene_start) %>% select(subtype,loci,numeric_chr,gene_start,gene_end,gene_name,ensg_id,gene_type,celltype,ancestry,zscore,pvalue) %>%
  separate(ensg_id, into = c("ENSG", "decimal"), remove = TRUE) %>%
  select(-decimal)

# writing out these annotated significant TWAS hits
write.table(annotated_bonf_significant_twas,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI.tsv",quote=F,row.names=F,col.names=T,sep="\t")


