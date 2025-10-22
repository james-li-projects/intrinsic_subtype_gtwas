library(dplyr)
library(data.table)
library(tidyr)

jnci_1 <- fread("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/01_spredixcan_eqtl_sqtl/output/spredixcan_eqtl_mashr/spredixcan_igwas_gtexmashrv8_intrinsic_subtype_1__PM__Breast_Mammary_Tissue.csv") %>% filter(pvalue<2.6E-6) %>% separate(gene,into=c("ENSG","decimal")) %>% select(ENSG,zscore)

jnci_2 <- fread("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/01_spredixcan_eqtl_sqtl/output/spredixcan_eqtl_mashr/spredixcan_igwas_gtexmashrv8_intrinsic_subtype_2__PM__Breast_Mammary_Tissue.csv") %>% filter(pvalue<2.6E-6) %>% separate(gene,into=c("ENSG","decimal")) %>% select(ENSG,zscore)

jnci_5 <- fread("/gpfs/data/huo-lab/jmcclellan/acat_brca_combined/01_spredixcan_eqtl_sqtl/output/spredixcan_eqtl_mashr/spredixcan_igwas_gtexmashrv8_intrinsic_subtype_5__PM__Breast_Mammary_Tissue.csv") %>% filter(pvalue<2.6E-6) %>% separate(gene,into=c("ENSG","decimal")) %>% select(ENSG,zscore)

ajhg_1 <- data.frame(
  ENSG = c(
    "ENSG00000179750", "ENSG00000152348", "ENSG00000205464", 
    "ENSG00000123106", "ENSG00000071967", "ENSG00000145451", 
    "ENSG00000186468"
  ),
  zscore = c(
    -5.41, -5.53, -6.08, -5.40, -5.59, -8.96, 6.17
  )
)

ncomm_5 <- data.frame(
  ENSG=c("ENSG00000163064"), zscore=c(-5.8)
)

# assembling preivous zscore dfs
previous_HRPOS_HER2NEG <- rbind(jnci_1,ajhg_1) %>% mutate(subtype="HRPOS_HER2NEG") %>%
  group_by(ENSG, subtype) %>%
  summarise(previous_mean_zscore = mean(zscore, na.rm = TRUE), .groups = "drop") %>% select(subtype,ENSG,previous_mean_zscore)
previous_HRPOS_HER2POS <- jnci_2 %>% mutate(subtype="HRPOS_HER2POS") %>%
  group_by(ENSG, subtype) %>%
  summarise(previous_mean_zscore = mean(zscore, na.rm = TRUE), .groups = "drop") %>% select(subtype,ENSG,previous_mean_zscore)
previous_HRNEG_HER2NEG <- rbind(jnci_5,ncomm_5) %>% mutate(subtype="HRNEG_HER2NEG") %>%
  group_by(ENSG, subtype) %>%
  summarise(previous_mean_zscore = mean(zscore, na.rm = TRUE), .groups = "drop") %>% select(subtype,ENSG,previous_mean_zscore)
previous_zscore_df <- rbind(
  previous_HRNEG_HER2NEG,
  previous_HRPOS_HER2POS,
  previous_HRPOS_HER2NEG
)

# obtaining currrent results
current_zscore_df<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/parsed_output/annotated_bonf_significant_twas_LOCI_SUBTYPE_COJO.tsv") %>% select(subtype,ENSG,zscore) %>%
  group_by(ENSG, subtype) %>%
  summarise(current_mean_zscore = mean(zscore, na.rm = TRUE), .groups = "drop")

# joined zscore df
joined_zscore_df <- inner_join(previous_zscore_df,current_zscore_df,by=c("subtype","ENSG"))

# plotting results
library(ggplot2)
library(ggpubr)

# Scatter plot with swapped x and y axes
p <- ggplot(joined_zscore_df, aes(x = current_mean_zscore, y = previous_mean_zscore)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  stat_cor(
    aes(label = paste("italic(r)^2~`=`~", signif(..r..^2, digits = 2), 
                      "*`,`~italic(P)~`=`~", signif(..p.., digits = 2))),
    method = "pearson",
    label.x = min(joined_zscore_df$current_mean_zscore),
    label.y = max(joined_zscore_df$previous_mean_zscore),
    size = 5,
    parse = TRUE
  ) +
  labs(x = "Mean Z-score of each gene across cell types and ancestries in the current study", 
       y = "Mean Z-score of each gene across previous studies") +
  theme_classic(base_size = 12)

# Save the plot as high-quality PNG
output_path <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/validation_previous_intrinsic_twas/zscore_scatter_plot.png"
ggsave(output_path, plot = p, width = 8, height = 6, dpi = 300)
