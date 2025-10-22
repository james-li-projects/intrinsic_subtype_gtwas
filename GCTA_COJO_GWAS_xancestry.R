# initializing packages
library(data.table)
library(dplyr)
library(ggplot2)
library(topr)
library(meta)
library(tidyr)

# set working directory
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis")

# importing EAF and N annotations for each ancestry's GWAS meta-analysis studies
afr_eaf_n<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/afr_eaf_n.txt") %>% select(-N)
eur_eaf_n<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/eur_eaf_n.txt") %>% select(-N)
xancestry_eaf_n <- inner_join(afr_eaf_n%>%rename(afr_EAF=EAF),eur_eaf_n%>%rename(eur_EAF=EAF),by=c("ID"))
xancestry_eaf_n <- xancestry_eaf_n %>% mutate(
  EAF=(eur_EAF*91477+afr_EAF*18800)/
    (91477+18800)) %>% select(-afr_EAF,-eur_EAF)

# initializing subtype
# subtype="HRPOS_HER2NEG"
# subtype="HRNEG_HER2NEG"
# subtype="HRPOS_HER2POS"
# subtype="HRNEG_HER2POS" # No variants

for (subtype in c(
  "HRPOS_HER2NEG",
  "HRNEG_HER2NEG",
  "HRPOS_HER2POS"
)) {
#  
#}

# importing summary statistics
eur_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_meta_sumstats/",subtype,".tsv")) %>% inner_join(eur_eaf_n,by=c("ID"))
afr_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats/",subtype,".tsv")) %>% inner_join(afr_eaf_n,by=c("ID"))

# computing effective sample sizes
eur_sumstats <- eur_sumstats %>% mutate(N= (-((1/91477)-2*SE^2*EAF*(1-EAF))^-1) + 91477) 
afr_sumstats <- afr_sumstats %>% mutate(N= (-((1/18800)-2*SE^2*EAF*(1-EAF))^-1) + 18800) 

# filtering summary statistics for available variants across all datasets
available_variants <- intersect(eur_sumstats$ID,afr_sumstats$ID)
input_eur_sumstats <- eur_sumstats %>% filter(ID%in%available_variants) 
input_afr_sumstats <- afr_sumstats %>% filter(ID%in%available_variants)
# input_xancestry_sumstats <- xancestry_sumstats %>% filter(ID%in%available_variants)

# adjusting the minimum p-values to be processable
input_eur_sumstats <- input_eur_sumstats %>% mutate(P=as.numeric(P)) %>% mutate(P = if_else(P < 3e-308, 3e-308, P))
input_afr_sumstats <- input_afr_sumstats %>% mutate(P=as.numeric(P)) %>% mutate(P = if_else(P < 3e-308, 3e-308, P))
# input_xancestry_sumstats <- input_xancestry_sumstats %>% mutate(P=as.numeric(P)) %>% mutate(P = if_else(P < 3e-308, 3e-308, P))


######################################
# IMPORTING LEAD VARIANTS TO ANALYZE #
######################################
target_xancestry_variant_df <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/",subtype,".txt")) %>% filter(Ancestry=="xancestry") %>% separate(ID,into=c("CHR","POS","A2","A1"),sep="\\:",remove=F)


##########################################
# INITIALIZING ANALYSIS FOR EACH VARIANT #
##########################################
# initializing DFs to store final COJO results
FINAL_cojo_lead_variant_output <- data.frame()
FINAL_cojo_joint_model_output <- data.frame()


# performing COJO for each variant 
for (i in 1:nrow(target_xancestry_variant_df)) {
# obtaining arguments for current lead variant being analyzed
current_variant_ID<-target_xancestry_variant_df$ID[i]
current_variant_POS<-as.numeric(target_xancestry_variant_df$POS[i])
current_variant_CHR<-as.numeric(target_xancestry_variant_df$CHR[i])

#########################
# preparing COJO inputs #
#########################
#######
# afr #
#######
# assembling COJO input file
tmp_cojo_input_afr <- input_afr_sumstats %>% filter(CHR==current_variant_CHR) %>% filter(POS>current_variant_POS-2e6) %>% filter(POS<current_variant_POS+2e6) %>% rename(A1=EffectAllele,A2=BaselineAllele,freq=EAF,b=BETA,se=SE,p=P) %>% select(ID,A1,A2,freq,b,se,p,N) 
write.table(tmp_cojo_input_afr,file="tmp_cojo_input_afr.input",quote=F,row.names=F,col.names=T,sep="\t")

# writing out targetted variant list for COJO 
write.table(tmp_cojo_input_afr%>%select(ID),"cojo_variant_list_afr.list",quote=F,row.names=F,col.names=F)

# subset bfile with only relevant variants for COJO
system(paste0("plink2 -bfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/afr --extract cojo_variant_list_afr.list --make-bed --out tmp_cojo_afr"))

#######
# eur #
#######
# assembling COJO input file
tmp_cojo_input_eur <- input_eur_sumstats %>% filter(CHR==current_variant_CHR) %>% filter(POS>current_variant_POS-2e6) %>% filter(POS<current_variant_POS+2e6) %>% rename(A1=EffectAllele,A2=BaselineAllele,freq=EAF,b=BETA,se=SE,p=P) %>% select(ID,A1,A2,freq,b,se,p,N) 
write.table(tmp_cojo_input_eur,file="tmp_cojo_input_eur.input",quote=F,row.names=F,col.names=T,sep="\t")

# writing out targetted variant list for COJO 
write.table(tmp_cojo_input_eur%>%select(ID),"cojo_variant_list_eur.list",quote=F,row.names=F,col.names=F)

# subset bfile with only relevant variants for COJO
system(paste0("plink2 -bfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/eur --extract cojo_variant_list_eur.list --make-bed --out tmp_cojo_eur"))


##################################################
# outputting initial cojo condition variant list #
##################################################
tmp_cojo_cond_file <- data.frame(ID=c(current_variant_ID))
write.table(tmp_cojo_cond_file,file="tmp_cojo_cond_variant_list.txt",quote=F,row.names=F,col.names=F)

#############################################################
# Running first-step COJO to condition only on lead variant #
#############################################################
# afr
system(paste0("gcta --bfile tmp_cojo_afr --cojo-file tmp_cojo_input_afr.input --cojo-cond tmp_cojo_cond_variant_list.txt --out ","tmp",".output"))
cojo_output_afr <- fread(paste0("tmp",".output.cma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bC,bC_se,pC) %>% rename(afr_freq=freq,afr_bC=bC,afr_bC_se=bC_se,afr_pC=pC)
# eur
system(paste0("gcta --bfile tmp_cojo_eur --cojo-file tmp_cojo_input_eur.input --cojo-cond tmp_cojo_cond_variant_list.txt --out ","tmp",".output"))
cojo_output_eur <- fread(paste0("tmp",".output.cma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bC,bC_se,pC) %>% rename(eur_freq=freq,eur_bC=bC,eur_bC_se=bC_se,eur_pC=pC)
# combining conditional results from both ancestries
cojo_output_both_ancestry <- inner_join(cojo_output_afr,cojo_output_eur,by=c("Chr","SNP","bp","refA"))

# performing meta-analysis of conditional results and initializing columns to store results 
cojo_output_both_ancestry$meta_cond_beta <- NA
cojo_output_both_ancestry$meta_cond_se <- NA
cojo_output_both_ancestry$meta_cond_p <- NA
for (m in 1:nrow(cojo_output_both_ancestry)) {
  # parsing inputs for meta-analysis of conditional estimates
  log_odds_ratios <- c(
    cojo_output_both_ancestry$afr_bC[m],
    cojo_output_both_ancestry$eur_bC[m])
  standard_errors <- c(
    cojo_output_both_ancestry$afr_bC_se[m],
    cojo_output_both_ancestry$eur_bC_se[m])
  # create a meta-analysis object
  meta_result <- metagen(
    TE = log_odds_ratios, 
    seTE = standard_errors, 
    sm = "OR", # Summary measure: Odds Ratio
    comb.fixed = TRUE, 
    comb.random = FALSE) # Fixed-effect only
  # extract results of interest
  meta_beta <- meta_result$TE.fixed   # Meta-analysis beta
  meta_se <- meta_result$seTE.fixed  # Standard error
  meta_z <- meta_beta / meta_se       # Z-score
  meta_p_value <- 2 * pnorm(-abs(meta_z)) # Two-sided p-value
  # add these values to the data.frame containing conditional estimates for both ancestries
  cojo_output_both_ancestry$meta_cond_beta[m] <- meta_beta
  cojo_output_both_ancestry$meta_cond_se[m] <- meta_se
  cojo_output_both_ancestry$meta_cond_p[m] <- meta_p_value
}

# saving a copy of this first step of COJO
first_step_cojo <- cojo_output_both_ancestry

# filtering for significant COJO variants following meta-analysis of these conditional estimates
cojo_output_both_ancestry <- cojo_output_both_ancestry %>% filter(meta_cond_p<1e-5) %>% arrange(meta_cond_p)

###############################################
# keep repeating cojo until there are no more independent signals
while (nrow(cojo_output_both_ancestry) > 0) {
  
  # updating the independent variant list and writing it out
  tmp_cojo_cond_file <- rbind(tmp_cojo_cond_file,cojo_output_both_ancestry %>% slice(1) %>% select(SNP) %>% rename(ID = SNP))
  write.table(tmp_cojo_cond_file,file = "tmp_cojo_cond_variant_list.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  ################
  # Running COJO #
  ################
  # afr
  system(paste0("gcta --bfile tmp_cojo_afr --cojo-file tmp_cojo_input_afr.input --cojo-cond tmp_cojo_cond_variant_list.txt --out ","tmp",".output"))
  cojo_output_afr <- fread(paste0("tmp",".output.cma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bC,bC_se,pC) %>% rename(afr_freq=freq,afr_bC=bC,afr_bC_se=bC_se,afr_pC=pC)
  # eur
  system(paste0("gcta --bfile tmp_cojo_eur --cojo-file tmp_cojo_input_eur.input --cojo-cond tmp_cojo_cond_variant_list.txt --out ","tmp",".output"))
  cojo_output_eur <- fread(paste0("tmp",".output.cma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bC,bC_se,pC) %>% rename(eur_freq=freq,eur_bC=bC,eur_bC_se=bC_se,eur_pC=pC)
  # combining conditional results from both ancestries
  cojo_output_both_ancestry <- inner_join(cojo_output_afr,cojo_output_eur,by=c("Chr","SNP","bp","refA"))
  
  # performing meta-analysis of conditional results and initializing columns to store results 
  cojo_output_both_ancestry$meta_cond_beta <- NA
  cojo_output_both_ancestry$meta_cond_se <- NA
  cojo_output_both_ancestry$meta_cond_p <- NA
  for (m in 1:nrow(cojo_output_both_ancestry)) {
    # parsing inputs for meta-analysis of conditional estimates
    log_odds_ratios <- c(
      cojo_output_both_ancestry$afr_bC[m],
      cojo_output_both_ancestry$eur_bC[m])
    standard_errors <- c(
      cojo_output_both_ancestry$afr_bC_se[m],
      cojo_output_both_ancestry$eur_bC_se[m])
    # create a meta-analysis object
    meta_result <- metagen(
      TE = log_odds_ratios, 
      seTE = standard_errors, 
      sm = "OR", # Summary measure: Odds Ratio
      comb.fixed = TRUE, 
      comb.random = FALSE) # Fixed-effect only
    # extract results of interest
    meta_beta <- meta_result$TE.fixed   # Meta-analysis beta
    meta_se <- meta_result$seTE.fixed  # Standard error
    meta_z <- meta_beta / meta_se       # Z-score
    meta_p_value <- 2 * pnorm(-abs(meta_z)) # Two-sided p-value
    # add these values to the data.frame containing conditional estimates for both ancestries
    cojo_output_both_ancestry$meta_cond_beta[m] <- meta_beta
    cojo_output_both_ancestry$meta_cond_se[m] <- meta_se
    cojo_output_both_ancestry$meta_cond_p[m] <- meta_p_value
  }
  
  # filtering for significant COJO variants following meta-analysis of these conditional estimates
  cojo_output_both_ancestry <- cojo_output_both_ancestry %>% filter(meta_cond_p<1e-5) %>% arrange(meta_cond_p)
}

print(tmp_cojo_cond_file)


##############################################
# obtaining joint model estimates for all conditionally independent variants
#######
# afr #
#######
system(paste0("gcta --bfile tmp_cojo_afr --cojo-file tmp_cojo_input_afr.input --extract tmp_cojo_cond_variant_list.txt --cojo-joint --out ","tmp",".output"))
cojo_joint_model_afr <- fread(paste0("tmp",".output.jma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bJ,bJ_se) %>% rename(afr_freq=freq,afr_bJ=bJ,afr_bJ_se=bJ_se)
#######
# eur #
#######
system(paste0("gcta --bfile tmp_cojo_eur --cojo-file tmp_cojo_input_eur.input --extract tmp_cojo_cond_variant_list.txt --cojo-joint --out ","tmp",".output"))
cojo_joint_model_eur <- fread(paste0("tmp",".output.jma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bJ,bJ_se) %>% rename(eur_freq=freq,eur_bJ=bJ,eur_bJ_se=bJ_se)

# combining joint model results from both ancestries
cojo_joint_model_both_ancestry <- inner_join(cojo_joint_model_afr,cojo_joint_model_eur,by=c("Chr","SNP","bp","refA"))

# performing meta-analysis of conditional results and initializing columns to store results 
cojo_joint_model_both_ancestry$meta_cond_beta <- NA
cojo_joint_model_both_ancestry$meta_cond_se <- NA
cojo_joint_model_both_ancestry$meta_cond_p <- NA
for (m in 1:nrow(cojo_joint_model_both_ancestry)) {
  # parsing inputs for meta-analysis of conditional estimates
  log_odds_ratios <- c(
    cojo_joint_model_both_ancestry$afr_bJ[m],
    cojo_joint_model_both_ancestry$eur_bJ[m])
  standard_errors <- c(
    cojo_joint_model_both_ancestry$afr_bJ_se[m],
    cojo_joint_model_both_ancestry$eur_bJ_se[m])
  # create a meta-analysis object
  meta_result <- metagen(
    TE = log_odds_ratios, 
    seTE = standard_errors, 
    sm = "OR", # Summary measure: Odds Ratio
    comb.fixed = TRUE, 
    comb.random = FALSE) # Fixed-effect only
  # extract results of interest
  meta_beta <- meta_result$TE.fixed   # Meta-analysis beta
  meta_se <- meta_result$seTE.fixed  # Standard error
  meta_z <- meta_beta / meta_se       # Z-score
  meta_p_value <- 2 * pnorm(-abs(meta_z)) # Two-sided p-value
  # add these values to the data.frame containing conditional estimates for both ancestries
  cojo_joint_model_both_ancestry$meta_cond_beta[m] <- meta_beta
  cojo_joint_model_both_ancestry$meta_cond_se[m] <- meta_se
  cojo_joint_model_both_ancestry$meta_cond_p[m] <- meta_p_value
}

######################################
# storing final desired cojo results #
######################################
tmp_FINAL_cojo_lead_variant_output <- data.frame(
  Chr = integer(),
  SNP = character(),
  bp = integer(),
  refA = character(),
  afr_freq = numeric(),
  afr_bC = numeric(),
  afr_bC_se = numeric(),
  afr_pC = numeric(),
  eur_freq = numeric(),
  eur_bC = numeric(),
  eur_bC_se = numeric(),
  eur_pC = numeric(),
  meta_cond_beta = numeric(),
  meta_cond_se = numeric(),
  meta_cond_p = numeric(),
  lead_variant = character(),
  stringsAsFactors = FALSE
)
tmp_FINAL_cojo_joint_model_output <- data.frame(
  Chr = integer(),
  SNP = character(),
  bp = integer(),
  refA = character(),
  afr_freq = numeric(),
  afr_bJ = numeric(),
  afr_bJ_se = numeric(),
  eur_freq = numeric(),
  eur_bJ = numeric(),
  eur_bJ_se = numeric(),
  meta_cond_beta = numeric(),
  meta_cond_se = numeric(),
  meta_cond_p = numeric(),
  lead_variant = character(),
  stringsAsFactors = FALSE
)

# obtain coefficients conditioned on lead variant
tmp_FINAL_cojo_lead_variant_output <- first_step_cojo %>% mutate(lead_variant=current_variant_ID) %>% filter(SNP%in%cojo_joint_model_both_ancestry$SNP)
# if there is at least one independent signal after COJO analysis then proceed
if (nrow(tmp_FINAL_cojo_lead_variant_output)>0) {
  # obtain joint model coefficients and add a column for the lead variant ID to the joint model coefficients DF
  tmp_FINAL_cojo_joint_model_output <- cojo_joint_model_both_ancestry %>% mutate(lead_variant=current_variant_ID)
}
# storing final results
FINAL_cojo_lead_variant_output <- rbind(
  FINAL_cojo_lead_variant_output,
  tmp_FINAL_cojo_lead_variant_output
)
FINAL_cojo_joint_model_output <- rbind(
  FINAL_cojo_joint_model_output,
  tmp_FINAL_cojo_joint_model_output
)
}

# writing out tables
write.table(FINAL_cojo_lead_variant_output%>%select(lead_variant, everything()),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/cojo_output/xancestry/",subtype,".cma.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
write.table(FINAL_cojo_joint_model_output%>%select(lead_variant, everything()),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/cojo_output/xancestry/",subtype,".jma.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
}
