library(topr)
library(data.table)
library(dplyr)
library(tidyr)
# library(GenomicRanges)
# library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

##################################
# import cytoband dictionary
cytoband_dictionary <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/ucsc.hg38.cytoband.txt",header=T)
colnames(cytoband_dictionary) <- c("Chr","Start","End","Band","Stain")

# define cytoband annotation function
annotate_cytoband <- function(input_df, cytoband_dictionary) {
  # Rename and convert position to integer
  df_results <- input_df %>%
    dplyr::rename(chr = CHROM, pos = POS) %>%
    dplyr::mutate(pos = as.integer(pos))
  
  # Preallocate LOCI column for performance
  df_results$LOCI <- NA
  
  for (i in seq_len(nrow(df_results))) {
    message(paste("Extracting cytoband:", i, "of", nrow(df_results)))
    chr_str <- df_results$chr[i]
    current_variant_pos <- df_results$pos[i]
    band <- dplyr::filter(cytoband_dictionary, Chr == chr_str, Start < current_variant_pos, End > current_variant_pos)$Band
    band <- tail(band, 1)
    cytoband <- paste0(chr_str, band)
    df_results$LOCI[i] <- gsub("chr","",cytoband)
  }
  
  # revert chr/pos column names in df_results
  df_results <- df_results %>% rename(CHROM=chr,POS=pos)
  return(df_results)
}

# importing rsid conversion dictionary
ID_RSID_dict <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/rsid_dictionary/ID_RSID_dict.tsv")

# import previous EUR variants identified in NG 2020 to identify variants for exclusion -- setting their pvalues very low so they'll get selected for overlapping loci with current variants below
EUR_lead_variants <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LIFTOVER_LEAD_VARIANTS_HG38/liftover_alleles.bed") 
EUR_lead_variants <- EUR_lead_variants %>% select(-V3) %>% dplyr::rename(CHROM=V1,POS=V2,REF=V4,ALT=V5) %>% mutate(BETA=0.1,SE=0.02,P=1e-300,CHR=gsub("chr","",CHROM)) %>% mutate(ID=paste(CHR,POS,REF,ALT,sep=":")) %>% select(CHROM,POS,ID,REF,ALT,BETA,SE,P) %>% mutate(ancestry="EUR")


##################################
# importing subtype lead variants for all subtypes except "HRNEG_HER2POS", which had no lead variants
for (subtype in c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2NEG")) {

# joining variants from our study and the 32 variants identified in NG 2020 -- setting pvalues to be more marginal
tmp_subtype_lead_variants <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/",subtype,".txt")) %>% separate(ID,sep="\\:",into=c("CHROM","POS","REF","ALT"),remove=F) %>% mutate(BETA=0.05,SE=0.02,P=1e-8) %>% select(CHROM,POS,ID,REF,ALT,BETA,SE,P) %>% mutate(ancestry="Novel")

# rbinding variants to input to identify lead variants
tmp_combined_variants = rbind(
  EUR_lead_variants,
  tmp_subtype_lead_variants
)

# identifying lead variants
tmp_novel_variant_list <- setdiff((get_lead_snps(tmp_combined_variants,thresh = 5e-08,region_size = 4e+06) %>% filter(ancestry=="Novel"))$ID,c("3:197512284:C:T","13:110604015:T:C","6:20534230:G:GA"))

# extracting lead variants from all study effects
tmp_filtered_subtype_lead_variants <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/",subtype,".txt")) %>% filter(ID%in%tmp_novel_variant_list)

# Function to convert BETA and SE to OR with 95% CI
logit_to_OR <- function(beta, se) {
  or <- exp(beta)
  lower <- exp(beta - 1.96 * se)
  upper <- exp(beta + 1.96 * se)
  sprintf("%.2f (%.2f, %.2f)", or, lower, upper)
}

# Use your actual data.table
dt <- tmp_filtered_subtype_lead_variants

# Extract all suffixes from BETA_* columns
suffixes <- unique(gsub("^BETA_", "", grep("^BETA_", names(dt), value = TRUE)))

# Loop through each suffix to create OR_* and reorder columns
for (suf in suffixes) {
  beta_col <- paste0("BETA_", suf)
  se_col <- paste0("SE_", suf)
  or_col <- paste0("OR_", suf)
  p_col <- paste0("P_", suf)
  
  # Only compute if both BETA and SE columns are present
  if (all(c(beta_col, se_col) %in% names(dt))) {
    dt[[or_col]] <- logit_to_OR(dt[[beta_col]], dt[[se_col]])
    
    # Reorder: place OR_* right before the matching P_*, if P_* exists
    if (p_col %in% names(dt)) {
      setcolorder(dt, c(setdiff(names(dt), or_col), or_col))  # Add at end first
      setcolorder(dt, append(setdiff(names(dt), or_col), or_col, after = which(names(dt) == p_col) - 1))
    }
  }
}

# removing columns with BETA or SE
dt <- dt %>%
  select(-matches("BETA|SE"))
rsid_annotated_dt <- inner_join(dt,ID_RSID_dict,by=c("ID"))
dt <- rsid_annotated_dt

# importing allele frequencies for ancestries and AABCG datasets
AFR_AF_DF <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/afr_eaf_n.txt") %>% select(-N) %>% rename(AF_AFR_META=EAF)
EUR_AF_DF <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/eur_eaf_n.txt") %>% select(-N) %>% rename(AF_EUR_META=EAF)
dt_by_ancestry_AF <- left_join(dt,AFR_AF_DF,by=c("ID")) %>% left_join(EUR_AF_DF,by=c("ID"))
by_study_AF <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/by_study_AF/merged_AF.tsv")
start_col <- which(names(by_study_AF) == "ALT") + 1
names(by_study_AF)[start_col:ncol(by_study_AF)] <- paste0("AF_", names(by_study_AF)[start_col:ncol(by_study_AF)])
by_study_AF <- by_study_AF%>% select(-`#CHROM`,-ALT) %>% rename(BaselineAllele=REF)
dt_by_ancestry_by_study_AF <- left_join(dt_by_ancestry_AF,by_study_AF,by=c("ID"))


# importing ancestry heterogeneity test p-values
heterogeneity_ancestry_df <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/ancestry_heterogeneity/heterogeneity_ancestry.tsv") %>% rename(ancestry_het_p=HetPVal)
dt_by_ancestry_by_study_het_AF <- left_join(dt_by_ancestry_by_study_AF,heterogeneity_ancestry_df,by=c("ID"))

# assigning cytobands
dt_by_ancestry_by_study_het_AF_CHROM_POS <- dt_by_ancestry_by_study_het_AF %>% separate(ID,into=c("CHROM","POS","REF","ALT"),sep="\\:",remove=F) %>% select(-REF,-ALT) %>% mutate(CHROM=paste0("chr",CHROM))
dt_by_ancestry_by_study_het_AF_CHROM_POS_CYTOBAND <- annotate_cytoband(dt_by_ancestry_by_study_het_AF_CHROM_POS,cytoband_dictionary) %>% select(-CHROM,-POS) 

# parse final lead variant data.frame to output 
all_studies_lead_variant_df <- copy(dt_by_ancestry_by_study_het_AF_CHROM_POS_CYTOBAND)
all_studies_lead_variant_df <- all_studies_lead_variant_df %>% mutate(SUBTYPE=subtype)
lead_cols <- c("SUBTYPE","LOCI", "rsid", "ID", "EffectAllele", "BaselineAllele")
suffixes <- c("AABC", "AMBER", "BCAC_OncoArray", "GBHS", "MEGA", "ROOT", "WGS",
              "AFR_META", "EUR_META", "XANCESTRY_META")
grouped_cols <- unlist(lapply(suffixes, function(s) {
  cols <- grep(paste0("_(?i)", s, "$"), names(all_studies_lead_variant_df), value = TRUE)
  c(grep("^OR_", cols, value = TRUE),
    grep("^P_" , cols, value = TRUE),
    grep("^AF_", cols, value = TRUE))
}))
setcolorder(all_studies_lead_variant_df, c(intersect(lead_cols, names(all_studies_lead_variant_df)), grouped_cols, setdiff(names(all_studies_lead_variant_df), c(lead_cols, grouped_cols))))

# importing subtype heterogeneity test p-values
AFR_SUBTYPE_HET_DF<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/subtype_heterogeneity/heterogeneity_table_afr.tsv") %>% select(-CHR,-POS) %>% rename(subtype_het_p=het_p)
XANCESTRY_SUBTYPE_HET_DF<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/subtype_heterogeneity/heterogeneity_table_xancestry.tsv") %>% select(-CHR,-POS) %>% rename(subtype_het_p=het_p)
combined_subtype_het_df <- rbind(
  AFR_SUBTYPE_HET_DF,
  XANCESTRY_SUBTYPE_HET_DF
)

# compiling final lead variant df
FINAL_LEAD_VARIANT_DF <- inner_join(all_studies_lead_variant_df,combined_subtype_het_df,by=c("ID","EffectAllele","BaselineAllele"))

# removing the NA (NA, NA) string 
FINAL_LEAD_VARIANT_DF[FINAL_LEAD_VARIANT_DF=="NA (NA, NA)"] <- "-"
FINAL_LEAD_VARIANT_DF[is.na(FINAL_LEAD_VARIANT_DF)] <- "-"

# writing out this final lead variant df 
write.table(FINAL_LEAD_VARIANT_DF,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/parsed_output/",subtype,"_FINAL_LEAD_VARIANT_DF.tsv"),quote=F,row.names=F,col.names=T,sep="\t")

}
