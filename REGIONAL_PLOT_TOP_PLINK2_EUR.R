
################################################
# PLOTTING REGIONAL PLOTS FOR CLUMPED VARIANTS #
################################################
subtype="HRNEG_HER2POS"

# importing EUR sumstats 
tmp_eur_sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_sumstats/", subtype, ".tsv")) %>% separate(ID1,into=c("chr","pos","a2","a1"),sep="\\:",remove=F)
tmp_eur_sumstats <- tmp_eur_sumstats %>% mutate(
  chr=as.numeric(chr),
  pos=as.numeric(pos)
)
tmp_eur_sumstats <- tmp_eur_sumstats %>% filter(!is.na(chr)) %>% filter(!is.na(pos))

# importing AFR sumstats and making alleles consistent with IDs
tmp_afr_sumstats_TOP <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/TOP/METAANALYSIS_",subtype,"_1.tbl")) 
tmp_afr_sumstats_TOP <- tmp_afr_sumstats_TOP %>% separate(MarkerName,into=c("chr","pos","a2","a1"),remove=F,sep="\\:")
tmp_afr_sumstats_TOP <- tmp_afr_sumstats_TOP %>% mutate(
  chr=as.numeric(chr),
  pos=as.numeric(pos)
)

# importing AFR sumstats and making alleles consistent with IDs
tmp_afr_sumstats_plink2 <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/plink2/METAANALYSIS_",subtype,"_1.tbl")) 
tmp_afr_sumstats_plink2 <- tmp_afr_sumstats_plink2 %>% separate(MarkerName,into=c("chr","pos","a2","a1"),remove=F,sep="\\:")
tmp_afr_sumstats_plink2 <- tmp_afr_sumstats_plink2 %>% mutate(
  chr=as.numeric(chr),
  pos=as.numeric(pos)
)



# importing clumped variants
clumped_df <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/meta_results/clumped/TOP/METAANALYSIS_",subtype,"_1.tbl.clumps"))

# importing MAF
afreq <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/input/afreq/",subtype,".afreq")) 
afreq <- afreq %>% mutate(MAF=ifelse(ALT_FREQS>0.5,1-ALT_FREQS,ALT_FREQS)) %>% rename(MarkerName=ID) %>% select(MarkerName,MAF) %>% filter(MarkerName%in%clumped_df$ID)
print(afreq)




# iterating for lead variants
for (current_variant in clumped_df$ID) {
  #current_variant="5:132367462:T:C"
  current_variant_pos=as.numeric(unlist(current_variant %>% strsplit(x=.,split="\\:"))[2])
  current_variant_chr=as.numeric(unlist(current_variant %>% strsplit(x=.,split="\\:"))[1])
  current_variant_lowerbound <- current_variant_pos - 250000
  current_variant_upperbound <- current_variant_pos + 250000
  
  # regional AFR TOP
  regional_afr_TOP <- tmp_afr_sumstats_TOP %>% filter(chr==current_variant_chr) %>% filter(pos>current_variant_lowerbound) %>% filter(pos<current_variant_upperbound) %>% rename(P=`P-value`) %>% select(pos,P) %>% mutate(Ancestry="African (TOP)")
  
  # regional AFR plink2
  regional_afr_plink2 <- tmp_afr_sumstats_plink2 %>% filter(chr==current_variant_chr) %>% filter(pos>current_variant_lowerbound) %>% filter(pos<current_variant_upperbound) %>% rename(P=`P-value`) %>% select(pos,P) %>% mutate(Ancestry="African (plink2)")
  
  # regional EUR 
  regional_eur <- tmp_eur_sumstats %>% filter(chr==current_variant_chr) %>% filter(pos>current_variant_lowerbound) %>% filter(pos<current_variant_upperbound) %>% select(pos,P) %>% mutate(Ancestry="European")
  
  # combined regional plot df
  combined_regional <- rbind(regional_afr_TOP,regional_afr_plink2,regional_eur)
  
  # Calculate the horizontal line threshold (-log10 of 1.25e-8)
  sig_threshold <- -log10(1.25e-8)
  p <- ggplot(combined_regional, aes(x = pos / 1e6, y = -log10(P))) +
    geom_point(alpha = 0.7, size = 1.5) +  
    geom_hline(yintercept = sig_threshold, color = "red", linetype = "dotted") +
    facet_wrap(~ Ancestry, ncol = 1) + 
    labs(x = paste0("Position on chromosome ", current_variant_chr, " (Mb)"),
         y = expression(-log[10](P))) +                               
    theme_bw() +                                                      
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))                        
  
  # Save the plot as a high quality PNG file
  ggsave(filename = paste0("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/QC/regional_plot_", subtype, "_", gsub("\\:",".",current_variant), ".png"),
         plot = p,
         width = 6,
         height = 5,    
         dpi = 300)      
}












