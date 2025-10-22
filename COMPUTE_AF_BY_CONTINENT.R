###########################
# importing libraries
library(readxl)
library(dplyr)
library(TOP)
library(data.table)
library(tidyr)
###########################

###########################
# subtype name conversion DF
subtype_conv <- data.frame(
  raw = c("HRPOS_HER2NEG", "HRPOS_HER2POS", "HRNEG_HER2POS", "HRNEG_HER2NEG"),
  final = c("Luminal-A-like", "Luminal-B-like", "HER2-enriched-like", "Triple-negative"),
  stringsAsFactors = FALSE
)
###########################

###########################
# reading in lead variants specific to AFR
lead_variant_path="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/all_by_ancestry/"
HRNEG_HER2NEG_afr <- fread(paste0(lead_variant_path,"/HRNEG_HER2NEG_afr.tsv"))
HRNEG_HER2NEG_eur <- fread(paste0(lead_variant_path,"/HRNEG_HER2NEG_eur.tsv"))
HRNEG_HER2NEG_xancestry <- fread(paste0(lead_variant_path,"/HRNEG_HER2NEG_xancestry.tsv"))
lead_variant_list <- rbind(
  HRNEG_HER2NEG_afr,
  HRNEG_HER2NEG_eur,
  HRNEG_HER2NEG_xancestry
  ) %>% unique()
write.table(lead_variant_list %>% select(ID),file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_variant.list",quote=F,row.names=F,col.names=F)
###########################

###########################
# importing covariate data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
# converting breast cancer status to a binary variable
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,`Dataset...3`,Age_GWAS,Status,ER,PR,HER2,AFR_pro,Study)
pheno_data$Status <- 2-pheno_data$Status
colnames(pheno_data)[1] <- "Sample_Name"
colnames(pheno_data)[2] <- "Dataset"

# importing data of principal components
eigenvec <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/PCA/AABCG_PCA_50PCs.eigenvec",header=T)
colnames(eigenvec)[1] <- "Sample_Name"

# joining covariate and PC data.frames by sample IID
combined_cov <- inner_join(pheno_data,eigenvec, by = c("Sample_Name")) 

# recoding continent variable
combined_cov$Continent <- "North America"
combined_cov$Continent[combined_cov$Study=="NBCS" | combined_cov$Study=="WAABCS" | combined_cov$Study=="GBHS"] <- "Africa"

sample_list_NA <- combined_cov %>% filter(Continent=="North America") %>% filter(Status==0) %>% mutate(F1=0) %>% select(F1,Sample_Name)
sample_list_Africa <- combined_cov %>% filter(Continent=="Africa") %>% filter(Status==0) %>% mutate(F1=0) %>% select(F1,Sample_Name)


###########################
# writing out these sample lists 
write.table(sample_list_NA,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_sample_list_NA.list",quote=F,row.names=F,col.names = F,sep="\t")
write.table(sample_list_Africa,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_sample_list_Africa.list",quote=F,row.names=F,col.names = F,sep="\t")


###########################
# computing allele frequencies
system("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS --extract /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_variant.list --keep /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_sample_list_NA.list --freq --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_AF_NA")
system("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/pfile/GWAS --extract /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_variant.list --keep /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_sample_list_Africa.list --freq --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_AF_Africa")
system("plink2 -bfile /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/eur --extract /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_variant.list --freq --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_AF_EUR")

###########################
# importing and joining allele frequencies
African<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_AF_Africa.afreq") %>% select(ID,ALT,REF,ALT_FREQS) %>% rename(African_ALT_FREQ=ALT_FREQS)
NorthAmerican<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/HRNEG_HER2NEG_AF_NA.afreq") %>% select(ID,ALT,REF,ALT_FREQS) %>% rename(NorthAmerican_ALT_FREQ=ALT_FREQS)
European <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/conditional_analysis/prepare_cojo_eaf_n/eur_eaf_n.txt") 
European <- European %>% dplyr::filter(ID%in%lead_variant_list$ID) %>% separate(ID,into=c("chr","pos","REF","ALT"),sep="\\:",remove=F) %>% rename(ALT_FREQS=EAF) %>% select(ID,ALT,REF,ALT_FREQS) %>% rename(European_ALT_FREQ=ALT_FREQS) 

# combining AFs across all continents
annotated_AF_df <- lead_variant_list %>% 
  select(ID,ALT,REF,BETA) %>% 
  full_join(African,by=c("ID","ALT","REF")) %>%
  full_join(NorthAmerican,by=c("ID","ALT","REF")) %>%
  full_join(European,by=c("ID","ALT","REF"))

###########################
library(data.table)

# Get names of all columns that contain "ALT_FREQ"
alt_freq_cols <- grep("ALT_FREQ", names(annotated_AF_df), value = TRUE)

# Identify rows with negative BETA
neg_beta_rows <- annotated_AF_df$BETA < 0

# For these rows, adjust ALT_FREQ columns: 1 - current value
annotated_AF_df[neg_beta_rows, (alt_freq_cols) := lapply(.SD, function(x) 1 - x), .SDcols = alt_freq_cols]

# Swap ALT and REF for those rows
annotated_AF_df[neg_beta_rows, `:=`(ALT = REF, REF = ALT)]

# Make BETA positive
annotated_AF_df[neg_beta_rows, BETA := abs(BETA)]

########################################################
# Replacing missing variant MAF based on manual search #
########################################################
# TNBC variants
annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "6:151634779:A:G"] <- 0.3330
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "6:151634779:A:G"] <- 0.3319

annotated_AF_df$European_ALT_FREQ[annotated_AF_df$ID == "1:31087245:A:G"]  <- 0.0001
annotated_AF_df$European_ALT_FREQ[annotated_AF_df$ID == "12:27719121:G:C"] <- 0.0002
annotated_AF_df$European_ALT_FREQ[annotated_AF_df$ID == "5:6120548:T:C"]   <- 0.002
# Luminal B variants
annotated_AF_df$European_ALT_FREQ[annotated_AF_df$ID == "10:123814971:T:C"] <- 0.0003
annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "22:28681643:C:CA"] <- 0
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "22:28681643:C:CA"] <- 0
# Luminal A variants
annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "1:241870961:A:G"] <- 0.0044
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "1:241870961:A:G"] <- 0.0046

annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "1:201468704:C:T"] <- 0.0124
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "1:201468704:C:T"] <- 0.013

annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "16:80616908:A:G"] <- 0.6627
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "16:80616908:A:G"] <- 0.6599

annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "17:42684371:C:T"] <- 0.004
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "17:42684371:C:T"] <- 0.0042

annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "19:13047463:T:C"] <- 0.0101
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "19:13047463:T:C"] <- 0.0105

annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "3:100037468:G:A"] <- 0.0721
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "3:100037468:G:A"] <- 0.0731

annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "5:50899259:A:AT"] <- 0
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "5:50899259:A:AT"] <- 0

annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "6:149265192:CT:C"] <- 1
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "6:149265192:CT:C"] <- 1

annotated_AF_df$African_ALT_FREQ[annotated_AF_df$ID == "1:113880354:T:A"] <- 1
annotated_AF_df$NorthAmerican_ALT_FREQ[annotated_AF_df$ID == "1:113880354:T:A"] <- 1

# write out table 
write.table(
  annotated_AF_df,
  file = "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/FINAL_HRNEG_HER2NEG_AF.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

library(data.table)
library(ggplot2)
library(dplyr)
library(ggsignif)  # for bracket annotations

# Load and clean data
annotated_AF_df <- fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/FINAL_HRNEG_HER2NEG_AF.tsv") %>% 
  filter(if_all(everything(), ~ !is.na(.)))

# Binomial tests
african_gt_eur <- sum(annotated_AF_df$African_ALT_FREQ > annotated_AF_df$European_ALT_FREQ)
total_rows <- nrow(annotated_AF_df)
african_vs_european_test <- binom.test(african_gt_eur, total_rows, p = 0.5, alternative = "greater")
na_gt_eur <- sum(annotated_AF_df$NorthAmerican_ALT_FREQ > annotated_AF_df$European_ALT_FREQ)
na_vs_european_test <- binom.test(na_gt_eur, total_rows, p = 0.5, alternative = "greater")

# wilcoxon signed rank tests
african_vs_european_test <- wilcox.test(annotated_AF_df$African_ALT_FREQ,annotated_AF_df$European_ALT_FREQ,paired=T)
na_vs_european_test <- wilcox.test(annotated_AF_df$NorthAmerican_ALT_FREQ,annotated_AF_df$European_ALT_FREQ,paired=T)

# Reshape and relabel
freq_long <- melt(annotated_AF_df,
                  measure.vars = c("African_ALT_FREQ", "NorthAmerican_ALT_FREQ", "European_ALT_FREQ"),
                  variable.name = "Population", value.name = "ALT_FREQ"
)[, Population := factor(Population,
                         levels = c("African_ALT_FREQ", "NorthAmerican_ALT_FREQ", "European_ALT_FREQ"),
                         labels = c("African (Africa)", "African (North America)", "European"))
]

# Format p-values in decimal
pval_africa <- paste0("p=", round(african_vs_european_test$p.value, 3))
pval_na <- paste0("p=", round(na_vs_european_test$p.value, 3))

# obtaining pretty subtype name
raw_subtype="HRNEG_HER2NEG"
final_subtype=(subtype_conv %>% filter(raw==raw_subtype))$final

# Generate the plot with binomial test brackets
p <- ggplot(freq_long, aes(x = Population, y = ALT_FREQ)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.2) +
  
  # Higher label for African (Africa) vs European
  geom_signif(comparisons = list(c("African (Africa)", "European")),
              annotations = pval_africa,
              y_position = max(freq_long$ALT_FREQ) + 0.125,  # higher
              tip_length = 0.04, vjust = -0.25, textsize = 5) +
  
  # Lower label for African (North America) vs European
  geom_signif(comparisons = list(c("African (North America)", "European")),
              annotations = pval_na,
              y_position = max(freq_long$ALT_FREQ) + 0.01,  # lower
              tip_length = 0.04, vjust = -0.25, textsize = 5) +
  
  labs(y = "Frequency of Risk Allele", x = "", title=final_subtype) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size=24)  # 0 = left, 0.5 = center, 1 = right
  )

# Save
ggsave("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/lead_variant_output/tables/risk_allele_AF/FINAL_HRNEG_HER2NEG_AF_boxplot.png",
       plot = p, width = 8, height = 6.25, dpi = 300)
