# initializing packages
library(data.table)
library(dplyr)
library(topr)
library(locusplotr)
library(ggplot2)
library(patchwork)
library(ggrepel)

# defining list of problematic variants to remove (identified by performing analysis and subsequently removing)
problem_variant_list <- c("3:197512284:C:T","13:110604015:T:C","6:20534230:G:GA")

# import ID to RSID conversion dictionary
ID_RSID_dict <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/rsid_dictionary/ID_RSID_dict.tsv")

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
    df_results$LOCI[i] <- cytoband
  }
  
  # revert chr/pos column names in df_results
  df_results <- df_results %>% rename(CHROM=chr,POS=pos)
  return(df_results)
}

# importing bim files for LD references
afr_bim <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/afr.bim")
eur_bim <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/eur.bim")

# importing recombination rates
recombination_rates <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/recombination_rates/recomb-hg38/genetic_map_GRCh38_merged.tab") %>% arrange(chrom,pos)

# subtype="HRPOS_HER2NEG"
# subtype="HRNEG_HER2NEG"
# subtype="HRPOS_HER2POS"
# subtype="HRNEG_HER2POS" # No variants

for (subtype in c(
  "HRPOS_HER2NEG",
  "HRNEG_HER2NEG",
  "HRPOS_HER2POS",
  "HRNEG_HER2POS"
)) {

# import sumstats
raw_eur_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_meta_sumstats/",subtype,".tsv"))
raw_afr_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats/",subtype,".tsv"))
raw_xancestry_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_xancestry_meta_sumstats/",subtype,".tsv"))

# filter xancestry sumstats to only have variants present in both ancestries
available_variants <- intersect(raw_eur_sumstats$ID,raw_afr_sumstats$ID)
raw_xancestry_sumstats <- raw_xancestry_sumstats %>% filter(ID%in%available_variants)

# modifying very small p-values to a defined p-value 
eur_sumstats <- raw_eur_sumstats %>% mutate(P=as.numeric(P)) %>%
  mutate(P = if_else(P < 3e-308, 3e-308, P))
afr_sumstats <- raw_afr_sumstats %>% mutate(P=as.numeric(P)) %>%
  mutate(P = if_else(P < 3e-308, 3e-308, P))
xancestry_sumstats <- raw_xancestry_sumstats %>% mutate(P=as.numeric(P)) %>%
  mutate(P = if_else(P < 3e-308, 3e-308, P))

# parsing summary statistics as input into the "topr" package which will be utilized to obtain lead variants 
if (nrow(afr_sumstats %>% filter(P<1.25e-08)) > 0) {
  input_afr <- afr_sumstats %>% mutate(CHROM=paste0("chr",CHR),REF=BaselineAllele,ALT=EffectAllele) %>% select(CHROM,POS,ID,REF,ALT,BETA,SE,P)
  lead_variants_afr <- get_lead_snps(input_afr,thresh = 1.25e-08,region_size = 4e+06)
  dim(lead_variants_afr)
  lead_variants_afr <- annotate_cytoband(lead_variants_afr, cytoband_dictionary) %>% filter(!(ID%in%problem_variant_list))
  dim(lead_variants_afr)
  write.table(lead_variants_afr %>% select(-P),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry/",subtype,"_afr.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
}

if (nrow(xancestry_sumstats %>% filter(P<1.25e-08)) > 0) {
  input_xancestry <- xancestry_sumstats %>% mutate(CHROM=paste0("chr",CHR),REF=BaselineAllele,ALT=EffectAllele) %>% select(CHROM,POS,ID,REF,ALT,BETA,SE,P)
  lead_variants_xancestry <- get_lead_snps(input_xancestry,thresh = 1.25e-08,region_size = 4e+06)
  dim(lead_variants_xancestry)
  lead_variants_xancestry <- annotate_cytoband(lead_variants_xancestry, cytoband_dictionary) %>% filter(!(ID%in%problem_variant_list))
  dim(lead_variants_xancestry)
  write.table(lead_variants_xancestry %>% select(-P),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry/",subtype,"_xancestry.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
}

if (nrow(eur_sumstats %>% filter(P<1.25e-08)) > 0) {
  input_eur <- eur_sumstats %>% mutate(CHROM=paste0("chr",CHR),REF=BaselineAllele,ALT=EffectAllele) %>% select(CHROM,POS,ID,REF,ALT,BETA,SE,P)
  lead_variants_eur <- get_lead_snps(input_eur,thresh = 1.25e-08,region_size = 4e+06)
  dim(lead_variants_eur)
  lead_variants_eur <- annotate_cytoband(lead_variants_eur, cytoband_dictionary) %>% filter(!(ID%in%problem_variant_list))
  dim(lead_variants_eur)
  write.table(lead_variants_eur %>% select(-P),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry/",subtype,"_eur.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
}

if (subtype!="HRNEG_HER2POS") {

# collating lead variant results to identify loci discovered in the meta-analysis that were not identified in EUR
input_lead_variants_parse_xancestry <- rbind(
  lead_variants_eur %>% select(-LOCI) %>% mutate(ancestry="eur"),
  lead_variants_xancestry %>% select(-LOCI) %>% mutate(ancestry="xancestry")
)
input_lead_variants_parse_xancestry <- input_lead_variants_parse_xancestry %>% mutate(P=ifelse(ancestry=="xancestry",1E-8,3e-308))
lead_variant_parse_xancestry <- get_lead_snps(input_lead_variants_parse_xancestry,thresh = 1.25e-08,region_size = 4e+06)

dim(lead_variant_parse_xancestry)
lead_variant_parse_xancestry <- annotate_cytoband(lead_variant_parse_xancestry, cytoband_dictionary)
dim(lead_variant_parse_xancestry)

# further collating lead variant results to identify loci discovered in only AFR
input_lead_variants_parse_afr <- rbind(
  lead_variants_afr %>% select(-LOCI) %>% mutate(ancestry="afr"),
  lead_variant_parse_xancestry %>% select(-LOCI)
)
input_lead_variants_parse_afr <- input_lead_variants_parse_afr %>% mutate(P=ifelse(ancestry=="afr",1.2499999E-8,P))
lead_variant_parse_afr_xancestry <-get_lead_snps(input_lead_variants_parse_afr,thresh = 1.25e-08,region_size = 4e+06)

dim(lead_variant_parse_afr_xancestry)
lead_variant_parse_afr_xancestry <- annotate_cytoband(lead_variant_parse_afr_xancestry, cytoband_dictionary)
dim(lead_variant_parse_afr_xancestry)

# only examining these variants identified uniquely in AFR or XAncestry meta-analyses
focus_variant_afr_xancestry <- lead_variant_parse_afr_xancestry %>% filter(ancestry %in% c("xancestry","afr"))

# output list of variants identified uniquely in AFR or XAncestry meta-analyses
write.table(focus_variant_afr_xancestry%>%filter(ancestry=="xancestry")%>%select(ID),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/lead_variant_lists/",subtype,".xancestry.txt"),quote=F,row.names=F,col.names=F)
write.table(focus_variant_afr_xancestry%>%filter(ancestry=="afr")%>%select(ID),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/lead_variant_lists/",subtype,".afr.txt"),quote=F,row.names=F,col.names=F)

###################################
# importing all study DFs to parse out study-specific effect sizes
AFR_META <- raw_afr_sumstats %>% rename(BETA_AFR_META=BETA,SE_AFR_META=SE,P_AFR_META=P) %>% select(-CHR,-POS) %>% select(-Direction) #%>% rename(AFR_Direction=Direction) 
EUR_META <- raw_eur_sumstats %>% rename(BETA_EUR_META=BETA,SE_EUR_META=SE,P_EUR_META=P) %>% select(-CHR,-POS)
XANCESTRY_META <- raw_xancestry_sumstats %>% rename(BETA_XANCESTRY_META=BETA,SE_XANCESTRY_META=SE,P_XANCESTRY_META=P) %>% select(-CHR,-POS) %>% select(-Direction) #%>% rename(XANCESTRY_Direction=Direction)
AABC <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_COLLATED/",subtype,".AABC.sumstats")) %>% rename(BETA_AABC=BETA,SE_AABC=SE,P_AABC=P)
AMBER <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_COLLATED/",subtype,".AMBER.sumstats")) %>% rename(BETA_AMBER=BETA,SE_AMBER=SE,P_AMBER=P)
BCAC_OncoArray <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_COLLATED/",subtype,".BCAC_OncoArray.sumstats")) %>% rename(BETA_BCAC_OncoArray=BETA,SE_BCAC_OncoArray=SE,P_BCAC_OncoArray=P)
GBHS <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_COLLATED/",subtype,".GBHS.sumstats")) %>% rename(BETA_GBHS=BETA,SE_GBHS=SE,P_GBHS=P)
MEGA <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_COLLATED/",subtype,".MEGA.sumstats")) %>% rename(BETA_MEGA=BETA,SE_MEGA=SE,P_MEGA=P)
ROOT <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_COLLATED/",subtype,".ROOT.sumstats")) %>% rename(BETA_ROOT=BETA,SE_ROOT=SE,P_ROOT=P)
WGS <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/consistent_allele/sumstats_COLLATED/",subtype,".WGS.sumstats")) %>% rename(BETA_WGS=BETA,SE_WGS=SE,P_WGS=P)

# retrieving association summary statistics for these lead variants from all possible studies 
lead_variant_study_effect_df <- data.frame()
for (j in 1:nrow(focus_variant_afr_xancestry)) {
  tmp_lead_variant_ID=focus_variant_afr_xancestry$ID[j]
  tmp_lead_variant_Ancestry=focus_variant_afr_xancestry$ancestry[j]
  # Helper function to extract alleles from ID
  extract_alleles <- function(variant_id) {
    parts <- unlist(strsplit(variant_id, ":"))
    list(baseline = parts[3], effect = parts[4])
  }
  alleles <- extract_alleles(tmp_lead_variant_ID)
  # Function to filter or fill missing record
  process_dataset <- function(dt, dataset_name, id_col = "ID") {
    result <- dt[get(id_col) == tmp_lead_variant_ID]
    if (nrow(result) == 0) {
      cols <- names(dt)
      new_row <- as.list(rep(NA, length(cols)))
      names(new_row) <- cols
      new_row[[id_col]] <- tmp_lead_variant_ID
      new_row[["EffectAllele"]] <- alleles$effect
      new_row[["BaselineAllele"]] <- alleles$baseline
      
      result <- as.data.table(new_row)
    }
    return(result)
  }
  # Apply to all datasets
  AFR_META_row         <- process_dataset(AFR_META, "AFR_META")
  EUR_META_row         <- process_dataset(EUR_META, "EUR_META")
  XANCESTRY_META_row   <- process_dataset(XANCESTRY_META, "XANCESTRY_META")
  AABC_row             <- process_dataset(AABC, "AABC")
  AMBER_row            <- process_dataset(AMBER, "AMBER")
  BCAC_OncoArray_row   <- process_dataset(BCAC_OncoArray, "BCAC_OncoArray")
  GBHS_row             <- process_dataset(GBHS, "GBHS")
  MEGA_row             <- process_dataset(MEGA, "MEGA")
  ROOT_row             <- process_dataset(ROOT, "ROOT")
  WGS_row              <- process_dataset(WGS, "WGS")
  # Now merge all rows together by ID, EffectAllele, BaselineAllele
  tmp_lead_variant_study_effect <- Reduce(function(x, y) merge(x, y, by = c("ID", "EffectAllele", "BaselineAllele"), all = TRUE), list(
    AFR_META_row,
    EUR_META_row,
    XANCESTRY_META_row,
    AABC_row,
    AMBER_row,
    BCAC_OncoArray_row,
    GBHS_row,
    MEGA_row,
    ROOT_row,
    WGS_row
  ))
  # Helper function to get sign string (+/-/?)
  get_sign_string <- function(vec) {
    sapply(vec, function(val) {
      if (is.na(val)) return("?")
      else if (val > 0) return("+")
      else if (val < 0) return("-")
      else return("0")
    }) |> paste0(collapse = "")
  }
  # Get AFR_Direction from selected BETA columns
  afr_beta_cols <- c("BETA_AABC", "BETA_AMBER", "BETA_BCAC_OncoArray", "BETA_GBHS", "BETA_MEGA", "BETA_ROOT", "BETA_WGS")
  tmp_lead_variant_study_effect[, AFR_Direction := get_sign_string(.SD), .SDcols = afr_beta_cols]
  # Get XANCESTRY_Direction from selected BETA columns
  xancestry_beta_cols <- c("BETA_AFR_META", "BETA_EUR_META")
  tmp_lead_variant_study_effect[, XANCESTRY_Direction := get_sign_string(.SD), .SDcols = xancestry_beta_cols]
  # adding ancestry column
  tmp_lead_variant_study_effect <- cbind(tmp_lead_variant_study_effect,data.frame(Ancestry=tmp_lead_variant_Ancestry))
  # joining results with combined DF
  lead_variant_study_effect_df <- rbind(lead_variant_study_effect_df,tmp_lead_variant_study_effect)
}

# writing out DF
write.table(lead_variant_study_effect_df,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/",subtype,".txt"),quote=F,row.names=F,col.names=T,sep="\t")


###################################
## PLOTTING CODE: +/- 2000000 BP ##
###################################
set_window_size=4000000
half_window_size=set_window_size/2

for (j in 1:nrow(focus_variant_afr_xancestry)) {
  print(paste("Processing variant",j,"out of",nrow(focus_variant_afr_xancestry),"variants"))
  # extracting features for a given variant
  tmp_CHROM=focus_variant_afr_xancestry$CHROM[j]
  tmp_POS=as.integer(focus_variant_afr_xancestry$POS[j])
  tmp_Ancestry=toupper(focus_variant_afr_xancestry$ancestry[j])
  tmp_lead_variant_ID=focus_variant_afr_xancestry$ID[j]
  
  ###########################
  # Assembling plotting dfs #
  ###########################
  plot_df_afr <- input_afr %>% filter(CHROM==tmp_CHROM) %>% filter(POS>tmp_POS-half_window_size) %>% filter(POS<tmp_POS+half_window_size) %>% mutate(Ancestry="African-ancestry")
  plot_df_eur <- input_eur %>% filter(CHROM==tmp_CHROM) %>% filter(POS>tmp_POS-half_window_size) %>% filter(POS<tmp_POS+half_window_size) %>% mutate(Ancestry="European-ancestry")
  plot_df_xancestry <- input_xancestry %>% filter(CHROM==tmp_CHROM) %>% filter(POS>tmp_POS-half_window_size) %>% filter(POS<tmp_POS+half_window_size) %>% mutate(Ancestry="Cross-ancestry")
  
  #######################################
  # GENERATING afr LD WITH LEAD VARIANT #
  #######################################
  if (tmp_lead_variant_ID %in% afr_bim$V2) {
    # writing out variant list 
    write.table(rbind((plot_df_afr%>%select(ID)),data.frame(ID=c(tmp_lead_variant_ID)))%>%unique(),file="/scratch/jll1/tmp/variant.list",quote=F,row.names=F,col.names=F)
    
    # computing LD
    max_snp_window<-nrow(plot_df_afr)
    system(paste0("module load plink/1.9; plink -bfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/afr --extract /scratch/jll1/tmp/variant.list --ld-snp ",tmp_lead_variant_ID," --ld-window-kb 4000 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --memory 100000 --out /scratch/jll1/tmp/tmp_LD"))
    # import pairwise correlations
    pairwise_cor <- fread("/scratch/jll1/tmp/tmp_LD.ld")
    # Filter the rows where either SNP_A or SNP_B equals tmp_lead_variant_ID
    filtered_df <- pairwise_cor[pairwise_cor$SNP_A == tmp_lead_variant_ID | pairwise_cor$SNP_B == tmp_lead_variant_ID, ]
    # Create a new column 'OTHER_VARIANT'
    filtered_df$OTHER_VARIANT <- ifelse(filtered_df$SNP_A == tmp_lead_variant_ID, filtered_df$SNP_B, filtered_df$SNP_A)
    # finalizing r2 df
    r2_df <- filtered_df %>% select(OTHER_VARIANT,R2) %>% rename(ID=OTHER_VARIANT,r2=R2)
    r2_df <- rbind(r2_df,data.frame(ID=tmp_lead_variant_ID,r2=1))
    # joining correlations with main df
    tmp_plot_df_afr <- left_join(plot_df_afr,r2_df,by=c("ID"))
    plot_df_afr <- tmp_plot_df_afr
  } else {
    plot_df_afr <- plot_df_afr %>% mutate(r2=NA)
  }
  
  #######################################
  # GENERATING eur LD WITH LEAD VARIANT #
  #######################################
  if (tmp_lead_variant_ID %in% eur_bim$V2) {
    # writing out variant list 
    write.table(rbind((plot_df_eur%>%select(ID)),data.frame(ID=c(tmp_lead_variant_ID)))%>%unique(),file="/scratch/jll1/tmp/variant.list",quote=F,row.names=F,col.names=F)
    
    # computing LD
    max_snp_window<-nrow(plot_df_eur)
    system(paste0("module load plink/1.9; plink -bfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/eur --extract /scratch/jll1/tmp/variant.list --ld-snp ",tmp_lead_variant_ID," --ld-window-kb 4000 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --memory 100000 --out /scratch/jll1/tmp/tmp_LD"))
    
    # import pairwise correlations
    pairwise_cor <- fread("/scratch/jll1/tmp/tmp_LD.ld")
    # Filter the rows where either SNP_A or SNP_B equals tmp_lead_variant_ID
    filtered_df <- pairwise_cor[pairwise_cor$SNP_A == tmp_lead_variant_ID | pairwise_cor$SNP_B == tmp_lead_variant_ID, ]
    # Create a new column 'OTHER_VARIANT'
    filtered_df$OTHER_VARIANT <- ifelse(filtered_df$SNP_A == tmp_lead_variant_ID, filtered_df$SNP_B, filtered_df$SNP_A)
    # finalizing r2 df
    r2_df <- filtered_df %>% select(OTHER_VARIANT,R2) %>% rename(ID=OTHER_VARIANT,r2=R2)
    r2_df <- rbind(r2_df,data.frame(ID=tmp_lead_variant_ID,r2=1))
    # joining correlations with main df
    tmp_plot_df_eur <- left_join(plot_df_eur,r2_df,by=c("ID"))
    plot_df_eur <- tmp_plot_df_eur
  } else {
    plot_df_eur <- plot_df_eur %>% mutate(r2=NA)
  }
  
  ############################################
  # SETTING META-ANALYZED CORRELATIONS TO NA #
  ############################################
  plot_df_xancestry <- plot_df_xancestry %>% mutate(r2=NA)
  
  ######################
  # RBIND ALL PLOT DFS #   
  ######################
  plot_df <- rbind(plot_df_afr,plot_df_eur,plot_df_xancestry)
  # initialize the tmp_recombination_rates data.frame
  tmp_recombination_rates <- recombination_rates %>% filter(chrom==tmp_CHROM)
  
  ##########################################
  # GENERATING COMPREHENSIVE REGIONAL PLOT # 
  ##########################################
  
  # Ensure Ancestry facet order is fixed
  plot_df <- plot_df %>%
    mutate(
      Ancestry = factor(Ancestry, levels = c("African-ancestry", "European-ancestry", "Cross-ancestry")),
      is_lead = ID == tmp_lead_variant_ID
    )
  
  # Create is_lead indicator
  plot_df <- plot_df %>%
    mutate(is_lead = ID == tmp_lead_variant_ID)
  
  # Define x-axis domain (in base pairs)
  pos_min <- min(plot_df$POS, na.rm = TRUE)
  pos_max <- max(plot_df$POS, na.rm = TRUE)
  
  # Filter recombination rates: include extra buffer points just outside POS range
  tmp_recomb_filtered <- tmp_recombination_rates %>%
    filter(pos >= pos_min | pos <= pos_max) %>%
    arrange(pos)
  
  extra_low <- tmp_recombination_rates %>%
    filter(pos < pos_min) %>%
    slice_max(pos, n = 1)
  
  extra_high <- tmp_recombination_rates %>%
    filter(pos > pos_max) %>%
    slice_min(pos, n = 1)
  
  tmp_recomb_final <- bind_rows(extra_low, tmp_recomb_filtered, extra_high) %>%
    arrange(pos)
  
  # Scale recombination rate to match -log10(P) axis
  logp_min <- 0
  logp_max <- ceiling(max(-log10(plot_df$P), na.rm = TRUE))
  
  # FIXED recombination axis range
  recomb_min <- 0
  recomb_max <- 100
  scale_factor <- (logp_max - logp_min) / (recomb_max - recomb_min)
  
  # Add scaled recombination rate for plotting
  tmp_recomb_final <- tmp_recomb_final %>%
    mutate(recomb_rate_scaled = (recomb_rate - recomb_min) * scale_factor + logp_min)

  # defining offsets for lead variant label text
  x_offset <- 0.02 * ((pos_max - pos_min) / 1e6)  # 2% of x-axis width (in Mb)
  y_offset <- 0.00 * (logp_max - logp_min)        # 0% of y-axis height
  
  # Create a separate dataframe with one row per ancestry for the lead SNP that includes rsID
  lead_labels_df <- plot_df %>%
    filter(is_lead) %>%
    group_by(Ancestry) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      x_label = POS / 1e6 + x_offset,
      y_label = -log10(P) + y_offset
    )
  current_rsid = (ID_RSID_dict %>% filter(ID==head(lead_labels_df$ID,1)))$rsid
  # lead_labels_df <- lead_labels_df %>% mutate(rsid=current_rsid)
  

  # Create base plot
  p <- ggplot(plot_df, aes(x = POS / 1e6, y = -log10(P))) +
    # Main points colored by continuous r2 values
    geom_point(data = filter(plot_df, !is_lead & !is.na(r2)),
               aes(color = r2),
               shape = 16, size = 1.8, alpha = 0.8) +
    # Points with NA r2 (black)
    geom_point(data = filter(plot_df, !is_lead & is.na(r2)),
               color = "black",
               shape = 16, size = 1.8, alpha = 0.8) +
    # Lead SNP (diamond, blue-violet)
    geom_point(data = filter(plot_df, is_lead),
               shape = 18, size = 3, color = "blueviolet") +
    # Recombination rate line (scaled to y-axis), neon bright blue, thinner
    geom_line(data = tmp_recomb_final,
              aes(x = pos / 1e6, y = recomb_rate_scaled),
              color = "#00FFFF", linewidth = 0.3, alpha = 0.9) +
    # Significance threshold line
    geom_hline(yintercept = -log10(1.25e-8), color = "red", alpha = 0.5) +
    # Facet by ancestry
    facet_wrap(~Ancestry, ncol = 1) +
    # Continuous color scale for r²
    scale_color_gradientn(
      name = expression("Linkage Disequilibrium (r"^2*")"),
      colours = c("darkblue", "lightblue", "chartreuse3", "orange", "red"),
      values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1)),
      limits = c(0, 1),
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 8,
        barheight = 0.6
      )
    ) +
    # Axis labels and secondary axis with fixed limits
    scale_y_continuous(
      name = expression(-log[10](italic(P))),
      limits = c(logp_min, logp_max),  # Force full height of y-axis
      sec.axis = sec_axis(~ (. - logp_min) / scale_factor + recomb_min,
                          name = "Recombination rate (cM/Mb)",
                          breaks = seq(0, 100, by = 20))
    ) +
    # Limit x-axis to plot_df range (but recomb line has buffer)
    coord_cartesian(xlim = c(pos_min, pos_max) / 1e6) +
    # Labels
    labs(
      x = paste("Position on Chromosome", gsub("chr", "", tmp_CHROM), "(in Mb)")
    ) +
    # Theme
    theme_classic() +
    theme(
      strip.background = element_rect(color = "black", fill = NA),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "top"
    )  +
    # Lead SNP label with a white box and line
    # Add segment (line) connecting the label to the lead SNP
    geom_segment(
      data = lead_labels_df,
      aes(x = x_label, xend = POS / 1e6,
          y = y_label, yend = -log10(P)),
      color = "black",
      linewidth = 0.3
    ) +
    # Add label in a white box
    geom_label(
      data = lead_labels_df,
      aes(x = x_label, y = y_label, label = current_rsid),
      size = 3.25,
      fontface = "bold",
      fill = "white",
      color = "black",
      label.size = 0.3,  # Border thickness around label
      hjust = 0,
      vjust = 0
    )
  
  # Save plot
  output_path <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/regional_plots/", subtype, "_", tmp_Ancestry, "_", tmp_CHROM, "_", tmp_POS, ".png")
  ggsave(output_path, plot = p, width = 8, height = 8, dpi = 1200)
  
  ####################################################
  # Creating a combined manhattan + local genes plot #
  ####################################################
  # modifying manhattan plot 
  p_variant <- p + theme(strip.text = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text.x  = element_text(size = 12), axis.text.y  = element_text(size = 12)) +
    # Add label in a white box
    geom_label(
      data = lead_labels_df,
      aes(x = x_label, y = y_label, label = current_rsid),
      size = 4,
      fontface = "bold",
      fill = "white",
      color = "black",
      label.size = 0.3,  # Border thickness around label
      hjust = 0,
      vjust = 0
    ) + theme(
      legend.text = element_text(size = 11),       # Legend labels
      legend.title = element_text(size = 11)      # Legend title
    )
  
  # defining gene plot
  p_gene <- gg_geneplot(chr=as.numeric(gsub("chr","",unique(plot_df$CHROM))), start=min(plot_df$POS), end=max(plot_df$POS), genome_build = "GRCh38", max_levels = 3) + theme_classic() + theme(
    strip.background = element_rect(color = "black", fill = NA),
    panel.border = element_rect(color = "black", fill = NA)
  ) + theme(
    axis.title.y = element_blank(),   # remove y-axis title
    axis.text.y = element_blank(),    # remove y-axis tick labels
    axis.ticks.y = element_blank(),    # remove y-axis tick marks
    axis.title.x = element_blank(),   # remove x-axis title
    axis.text.x = element_blank(),    # remove x-axis tick labels
    axis.ticks.x = element_blank()    # remove x-axis tick marks
  ) +
    facet_wrap(~"Genes in this region") +
    theme(
      strip.text = element_text(size = 10),
      strip.background = element_rect(
        fill = "white",         # background color
        color = "black",        # border color
        linewidth = 0.5         # border thickness
      ),
      strip.placement = "outside",
      panel.spacing = unit(0, "lines")
    )
  
  # outputting combined plot, stacking the manhattan and gene plots vertically
  output_path_combined <- paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/regional_plots/", subtype, "_", tmp_Ancestry, "_", tmp_CHROM, "_", tmp_POS, "_combined_plot.png")
  
  # Save them stacked vertically
  png(file=output_path_combined, width = 10, height = 12, res = 1200, unit="in")
  print((p_variant / p_gene)  + plot_layout(heights = c(2, 1/3)))
  dev.off()
  
}
} else {
  print("Subtype is HRNEG_HER2POS and has no novel lead variants")
}
}
