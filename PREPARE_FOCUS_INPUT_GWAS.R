library(data.table)
library(dplyr)

# defining lists to iterate through
subtype_list = c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRNEG_HER2NEG")
race_list = c("BLACK","WHITE")
ancestry_list = c("afr","eur")
celltype_list = c("Adipocytes","Breast_tissue","Endothelial_cells","Epithelial_cells","Stromal_and_Immune_cells")

# process GWAS sumstats for each iteration
for (subtype in subtype_list) {
  for (i in 1:length(race_list)) {
    race=race_list[i]
    ancestry=ancestry_list[i]
    # importing summary statistics
    sumstats_filename=paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_sumstats/",ancestry,"/",subtype,".tsv")
    print(sumstats_filename)
    dt <- fread(sumstats_filename)
    dt <- dt %>% rename(BP=POS,SNP=ID,A1=EffectAllele,A2=BaselineAllele) %>% mutate(Z=BETA/SE) %>% select(CHR,BP,SNP,A1,A2,BETA,SE,Z,N)
    dt <- dt %>% mutate(N=as.integer(round(N)))
    # write the resulting DT as a gzipped tsv
    fwrite(
      dt,
      file = paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/FOCUS/input_gwas/",subtype,".",race,".tsv.gz"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      compress = "gzip"
    )
  }
}

# 11:35 
# 11:50


