library(data.table)
library(dplyr)

# defining ancestry and subtype lists
ancestry_list = c("afr","eur")
race_list = c("BLACK","WHITE")
subtype_list = c("HRPOS_HER2NEG","HRPOS_HER2POS","HRNEG_HER2POS","HRNEG_HER2NEG")

n_eff_df <- data.frame()

# for each ancestry and subype, compute mean effective sample size based on variant standard errors
for (i in 1:length(ancestry_list)) {
  for (tmp_subtype in subtype_list) {
    tmp_ancestry=ancestry_list[i]
    tmp_race=race_list[i]

    tmp_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/cojo_input_sumstats/",tmp_ancestry,"/",tmp_subtype,".tsv"))
    tmp_mean_N_eff <- mean(tmp_sumstats$N) 
    tmp_n_eff_df <- data.frame(ancestry=tmp_ancestry,race=tmp_race,subtype=tmp_subtype,n_eff=tmp_mean_N_eff)
    
    print(tmp_n_eff_df)
    
    n_eff_df <- rbind(
      n_eff_df,
      tmp_n_eff_df
    )
  }
}

# saving output sample sizes
save(n_eff_df,file="/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/mean_n_eff/n_eff_df.RData")
