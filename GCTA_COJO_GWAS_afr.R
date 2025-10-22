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
afr_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats/",subtype,".tsv")) %>% inner_join(afr_eaf_n,by=c("ID"))

# computing effective sample sizes
afr_sumstats <- afr_sumstats %>% mutate(N= (-((1/18800)-2*SE^2*EAF*(1-EAF))^-1) + 18800) 

# adjusting the minimum p-values to be processable
input_afr_sumstats <- afr_sumstats %>% mutate(P=as.numeric(P)) %>% mutate(P = if_else(P < 3e-308, 3e-308, P))


######################################
# IMPORTING LEAD VARIANTS TO ANALYZE #
######################################
target_afr_variant_df <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/study_effects/",subtype,".txt")) %>% filter(Ancestry=="afr") %>% separate(ID,into=c("CHR","POS","A2","A1"),sep="\\:",remove=F)


##########################################
# INITIALIZING ANALYSIS FOR EACH VARIANT #
##########################################
# initializing DFs to store final COJO results
FINAL_cojo_lead_variant_output <- data.frame()
FINAL_cojo_joint_model_output <- data.frame()


# performing COJO for each variant 
for (i in 1:nrow(target_afr_variant_df)) {
# obtaining arguments for current lead variant being analyzed
current_variant_ID<-target_afr_variant_df$ID[i]
current_variant_POS<-as.numeric(target_afr_variant_df$POS[i])
current_variant_CHR<-as.numeric(target_afr_variant_df$CHR[i])

#########################
# preparing COJO inputs #
#########################
#######
# afr #
#######
# assembling COJO input file
tmp_cojo_input_afr <- input_afr_sumstats %>% filter(CHR==current_variant_CHR) %>% filter(POS>current_variant_POS-2e6) %>% filter(POS<current_variant_POS+2e6) %>% rename(ALLELE1=EffectAllele,ALLELE0=BaselineAllele,A1FREQ=EAF) %>% select(ID,ALLELE1,ALLELE0,A1FREQ,BETA,SE,P,N) 
write.table(tmp_cojo_input_afr,file="tmp_cojo_input_afr.input",quote=F,row.names=F,col.names=T,sep="\t")

# writing out targetted variant list for COJO 
write.table(tmp_cojo_input_afr%>%select(ID),"cojo_variant_list_afr.list",quote=F,row.names=F,col.names=F)

# subset bfile with only relevant variants for COJO
system(paste0("plink2 -bfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/afr --extract cojo_variant_list_afr.list --make-bed --out tmp_cojo_afr"))


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
cojo_output_afr <- fread(paste0("tmp",".output.cma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bC,bC_se,pC) %>% rename(afr_freq=freq,cond_beta=bC,cond_se=bC_se,cond_p=pC)

# saving a copy of this first step of COJO
first_step_cojo <- cojo_output_afr

# filtering for significant COJO variants following meta-analysis of these conditional estimates
cojo_output_afr <- cojo_output_afr %>% filter(cond_p<1e-5) %>% arrange(cond_p)

###############################################
# keep repeating cojo until there are no more independent signals
while (nrow(cojo_output_afr) > 0) {
  
  # updating the independent variant list and writing it out
  tmp_cojo_cond_file <- rbind(tmp_cojo_cond_file,cojo_output_afr %>% slice(1) %>% select(SNP) %>% rename(ID = SNP))
  write.table(tmp_cojo_cond_file,file = "tmp_cojo_cond_variant_list.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  ################
  # Running COJO #
  ################
  # afr
  system(paste0("gcta --bfile tmp_cojo_afr --cojo-file tmp_cojo_input_afr.input --cojo-cond tmp_cojo_cond_variant_list.txt --out ","tmp",".output"))
  cojo_output_afr <- fread(paste0("tmp",".output.cma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bC,bC_se,pC) %>% rename(afr_freq=freq,cond_beta=bC,cond_se=bC_se,cond_p=pC)
 
  # filtering for significant COJO variants following meta-analysis of these conditional estimates
  cojo_output_afr <- cojo_output_afr %>% filter(cond_p<1e-5) %>% arrange(cond_p)
}

print(tmp_cojo_cond_file)


##############################################
# obtaining joint model estimates for all conditionally independent variants
#######
# afr #
#######
system(paste0("gcta --bfile tmp_cojo_afr --cojo-file tmp_cojo_input_afr.input --extract tmp_cojo_cond_variant_list.txt --cojo-joint --out ","tmp",".output"))
cojo_joint_model_afr <- fread(paste0("tmp",".output.jma.cojo")) %>% select(Chr,SNP,bp,refA,freq,bJ,bJ_se,pJ) %>% rename(afr_freq=freq,cond_beta=bJ,cond_se=bJ_se,cond_p=pJ)

######################################
# storing final desired cojo results #
######################################
tmp_FINAL_cojo_lead_variant_output <- data.frame(
  Chr = integer(),
  SNP = character(),
  bp = integer(),
  refA = character(),
  afr_freq = numeric(),
  cond_beta = numeric(),
  cond_se = numeric(),
  cond_p = numeric(),
  lead_variant = character(),
  stringsAsFactors = FALSE
)
tmp_FINAL_cojo_joint_model_output <- data.frame(
  Chr = integer(),
  SNP = character(),
  bp = integer(),
  refA = character(),
  afr_freq = numeric(),
  cond_beta = numeric(),
  cond_se = numeric(),
  cond_p = numeric(),
  lead_variant = character(),
  stringsAsFactors = FALSE
)

# obtain coefficients conditioned on lead variant
tmp_FINAL_cojo_lead_variant_output <- first_step_cojo %>% mutate(lead_variant=current_variant_ID) %>% filter(SNP%in%cojo_joint_model_afr$SNP)
# if there is at least one independent signal after COJO analysis then proceed
if (nrow(tmp_FINAL_cojo_lead_variant_output)>0) {
  # obtain joint model coefficients and add a column for the lead variant ID to the joint model coefficients DF
  tmp_FINAL_cojo_joint_model_output <- cojo_joint_model_afr %>% mutate(lead_variant=current_variant_ID)
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
write.table(FINAL_cojo_lead_variant_output%>%select(lead_variant, everything()),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/cojo_output/afr/",subtype,".cma.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
write.table(FINAL_cojo_joint_model_output%>%select(lead_variant, everything()),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/cojo_output/afr/",subtype,".jma.tsv"),quote=F,row.names=F,col.names=T,sep="\t")
}
