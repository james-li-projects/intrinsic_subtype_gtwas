library(DBI)
library(RSQLite)
library(dplyr)

# set working directory
setwd("/gpfs/data/huo-lab/Vanderbilt/Julian/02_merge_in_gtex_and_prep_data/output/predixcan_dbs")

# List all gene model .db files
db_files <- list.files(pattern = "\\.gene_models\\..*\\.db$")

# Initialize empty list to hold data frames
all_extra_dfs <- list()

# Loop over each db file
for (file in db_files) {
  # Connect to the SQLite database
  con <- dbConnect(SQLite(), file)
  
  # Try to read the "extra" table
  if ("extra" %in% dbListTables(con)) {
    df <- dbReadTable(con, "extra") %>%
      select(pred.perf.R2) %>%
      mutate(
        ancestry = ifelse(grepl("BLACK", file), "African-ancestry", "European-ancestry"),
        celltype = sub(".*\\.gene_models\\.(.*)\\.db", "\\1", file)
      )
    all_extra_dfs[[file]] <- df
  } else {
    warning(paste("No 'extra' table in", file))
  }
  
  # Disconnect after reading
  dbDisconnect(con)
}

# Combine all data frames into one
combined_df <- bind_rows(all_extra_dfs)

# View first few rows (optional)
head(combined_df)

# parsing this R2 table 
combined_df <- combined_df %>% filter(celltype!="Stromal_and_Immune_cells") %>% mutate(celltype=gsub("_"," ",celltype)) %>% mutate(celltype=ifelse(celltype=="Breast tissue","Bulk tissue",celltype))
# Set desired order of celltypes
desired_order <- c("Bulk tissue", "Epithelial cells", "Adipocytes", "Endothelial cells", "Fibroblasts")
# Apply factor level ordering to combined_df
combined_df$celltype <- factor(combined_df$celltype, levels = desired_order)
print(levels(combined_df$celltype))
combined_df$pred.perf.R2 <- as.numeric(combined_df$pred.perf.R2)


# plotting R2 values for all successfully modeled genes
library(ggplot2)
library(dplyr)
library(scales)  # for comma formatting

# Corrected number of genes tested (based on revised description)
tested_genes_df <- tibble::tribble(
  ~ancestry,           ~celltype,           ~n_tested,
  "European-ancestry", "Bulk tissue",        18083,
  "European-ancestry", "Epithelial cells",    8900,
  "European-ancestry", "Adipocytes",          8965,
  "European-ancestry", "Endothelial cells",   9089,
  "European-ancestry", "Fibroblasts",         9251,
  "African-ancestry",  "Bulk tissue",        20203,
  "African-ancestry",  "Epithelial cells",    9819,
  "African-ancestry",  "Adipocytes",          9402,
  "African-ancestry",  "Endothelial cells",   9621,
  "African-ancestry",  "Fibroblasts",        10148
)

# Compute number of genes modeled (R² ≥ 0.01)
label_data <- combined_df %>%
  filter(!is.na(pred.perf.R2)) %>%
  group_by(ancestry, celltype) %>%
  summarise(n_modeled = sum(pred.perf.R2 >= 0.01), .groups = "drop") %>%
  left_join(tested_genes_df, by = c("ancestry", "celltype")) %>%
  mutate(
    label = paste0(
      "Genes modeled: ", comma(n_modeled),"\n",
      "Genes tested: ", comma(n_tested)
    )
  )
# Apply factor level ordering to label_data
label_data$celltype  <- factor(label_data$celltype,  levels = desired_order)

# Set output path
output_path <- "/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/performance_metrics/rsq_histogram_by_ancestry_celltype.png"

# Plot histogram with text labels
p <- ggplot(combined_df, aes(x = sqrt(pred.perf.R2))) +
  geom_histogram(bins = 50, fill = "grey60", color = "black") +
  geom_text(
    data = label_data,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.05, vjust = 1.5,
    inherit.aes = FALSE,
    size = 3.5
  ) +
  facet_grid(rows = vars(ancestry), cols = vars(celltype), scales = "free_x") +
  theme_classic() +
  labs(
    # x = expression("Prediction Model " * R^2),
    x = expression("Prediction Model Correlation Coefficient (r)"),
    y = "Count"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) 

# Save to file
ggsave(filename = output_path, plot = p, width = 15, height = 6, dpi = 300)
