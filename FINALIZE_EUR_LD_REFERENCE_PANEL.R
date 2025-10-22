library(data.table)
library(dplyr)
library(stringr)

# reading in variant annotations for eur sumstats 
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_eur_meta_sumstats")
eur_variant_df <- data.frame()
for (file in list.files()) {
  print(file)
  tmp_eur_variant_df <- fread(file) %>% select(ID)
  eur_variant_df <- rbind(eur_variant_df,tmp_eur_variant_df)
}
eur_variant_df <- unique(eur_variant_df)

# reading in variant annotations for afr sumstats 
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/GWAS/processed_afr_meta_sumstats")
afr_variant_df <- data.frame()
for (file in list.files()) {
  print(file)
  tmp_afr_variant_df <- fread(file) %>% select(ID)
  afr_variant_df <- rbind(afr_variant_df,tmp_afr_variant_df)
}
afr_variant_df <- unique(afr_variant_df)

# combining both ancestry variant ID
combined_variant_df <- unique(rbind(afr_variant_df,eur_variant_df))

# importing bim
bim<-fread("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/1kg/chr_prefix_EUR_1kg.multi.norm.bim")

# assembling conversion dictionary
bim_orient_1 <- bim %>% mutate(ID=paste(V1,V4,V5,V6,sep=":")) %>% select(ID,V2)
bim_orient_2 <- bim %>% mutate(ID=paste(V1,V4,V6,V5,sep=":")) %>% select(ID,V2)
filtered_bim_orient_1 <- bim_orient_1 %>% filter(ID %in% combined_variant_df$ID)
filtered_bim_orient_2 <- bim_orient_2 %>% filter(ID %in% combined_variant_df$ID)
filtered_bim_all_orient <- rbind(
  filtered_bim_orient_1,
  filtered_bim_orient_2)

#############################################
# identify duplicated rows based on V2        
duplicated_rows <- filtered_bim_all_orient[filtered_bim_all_orient$V2 %in% filtered_bim_all_orient$V2[duplicated(filtered_bim_all_orient$V2)], ] %>% arrange(V2)
allele_match_same_order <- function(id, v2) {
  id_parts <- unlist(strsplit(id, ":"))
  v2_parts <- unlist(strsplit(v2, ":"))
  return(id_parts[2] == v2_parts[2] && id_parts[3] == v2_parts[3] && id_parts[4] == v2_parts[4])
}
duplicated_rows$keep <- mapply(allele_match_same_order, duplicated_rows$ID, duplicated_rows$V2)
correct_order_rows <- duplicated_rows[duplicated_rows$keep == TRUE, ]
correct_order_rows$keep <- NULL

# identify rows with unique records based on V2
unique_rows <- filtered_bim_all_orient[, .N, by = V2][N == 1]
unique_rows <- filtered_bim_all_orient[V2 %in% unique_rows$V2]

# combining all these rows based on V2 duplicate filtering
filtered_bim_all_orient <- rbind(unique_rows,correct_order_rows)

#############################################
# identify duplicated rows based on ID        
duplicated_rows <- filtered_bim_all_orient[filtered_bim_all_orient$ID %in% filtered_bim_all_orient$ID[duplicated(filtered_bim_all_orient$ID)], ] %>% arrange(ID)
allele_match_same_order <- function(v2, id) {
  v2_parts <- unlist(strsplit(v2, ":"))
  id_parts <- unlist(strsplit(id, ":"))
  return(v2_parts[2] == id_parts[2] && v2_parts[3] == id_parts[3] && v2_parts[4] == id_parts[4])
}
duplicated_rows$keep <- mapply(allele_match_same_order, duplicated_rows$V2, duplicated_rows$ID)
correct_order_rows <- duplicated_rows[duplicated_rows$keep == TRUE, ]
correct_order_rows$keep <- NULL

# identify rows with unique records based on ID
unique_rows <- filtered_bim_all_orient[, .N, by = ID][N == 1]
unique_rows <- filtered_bim_all_orient[ID %in% unique_rows$ID]

# assembling final conversion dictionary
final_conversion_dict <- rbind(unique_rows,correct_order_rows) %>% select(V2,ID)

# extract variants from bfile that can be converted
setwd("/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/1kg/")
write.table(final_conversion_dict %>% select(V2),file="extract_variant.list",quote=F,row.names=F,col.names=F)
system("plink2 -bfile chr_prefix_EUR_1kg.multi.norm --extract extract_variant.list --make-bed --out chr_prefix_EUR_1kg.multi.norm.filtered")

# assemble finalized and renamed EUR reference panel
write.table(final_conversion_dict,file="variant_id_conversion.dict",quote=F,row.names=F,col.names=F,sep="\t")
system("plink2 -bfile chr_prefix_EUR_1kg.multi.norm.filtered --update-name variant_id_conversion.dict --make-bed --out FINAL_EUR_1kg")

# making a copy in a designated LD reference panel folder
system("plink2 -bfile chr_prefix_EUR_1kg.multi.norm.filtered --update-name variant_id_conversion.dict --make-bed --out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/eur")

