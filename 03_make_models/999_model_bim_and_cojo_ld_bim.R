main <- function(bim_path1, bim_path2, output_name, remove_b1_chr=T) {
	library(tidyverse)

  # Step 1: Load the two BIM files
  bim1 <- read_delim(bim_path1, delim = "\t", col_names = FALSE)
  bim2 <- read_delim(bim_path2, delim = "\t", col_names = FALSE)
  
  # Step 2: Assign proper column names to the BIM files (based on standard BIM format)
  colnames(bim1) <- c("chr", "snp", "cm", "pos", "allele1", "allele2")
  colnames(bim2) <- c("chr", "snp", "cm", "pos", "allele1", "allele2")

	if (isTRUE(remove_b1_chr)){
		bim1 <- bim1 %>%
			mutate(snp = str_replace(snp, "chr", ""))
	}
  
  # Step 3: Merge both BIM files on chromosome and position
  merged_bim <- left_join(bim1, bim2, by = c("chr", "pos"), suffix = c("_bim1", "_bim2"))
  
  # Step 4: Filter based on the conditions
  # Condition 1: Same chromosome, position, and alleles
  condition1 <- merged_bim %>%
    filter(allele1_bim1 == allele1_bim2 & allele2_bim1 == allele2_bim2)
  
  # Condition 2: Same chromosome, position, but flipped alleles
  condition2 <- merged_bim %>%
    filter(allele1_bim1 == allele2_bim2 & allele2_bim1 == allele1_bim2)
  
  # Step 5: Combine results from both conditions and remove duplicates
  result <- bind_rows(condition1, condition2) %>%
    select(snp_bim1, snp_bim2) %>%
    distinct()  # Remove duplicate SNP pairs
  
  # Step 6: Write the result to a tab-delimited file
  write_delim(result, output_name, delim = "\t")
  
  # Return the result
  # return(result)
}

if (interactive()){
# 	main("../input/vander_gtex_genotypes/WHITE.chrALL_vander_gtex.bim",
# 			 "../input/cojo_ld_ref_panels/eur.bim", 
# 			 "../output/WHITE.vander_gtex_model_ids_to_james_cojo_ld_ref.tsv")
# 	main("../input/vander_gtex_genotypes/BLACK.chrALL_vander_gtex.bim",
# 			 "../input/cojo_ld_ref_panels/afr.bim", 
# 			 "../output/BLACK.vander_gtex_model_ids_to_james_cojo_ld_ref.tsv")

	main(
			 "../input/cojo_ld_ref_panels/eur.bim", 
			 "../input/vander_gtex_genotypes/WHITE.chrALL_vander_gtex.bim",

			 "../output/WHITE.james_cojo_ld_ref_ids_to_vander_gtex.tsv")
	main(
			 "../input/cojo_ld_ref_panels/afr.bim", 
			 "../input/vander_gtex_genotypes/BLACK.chrALL_vander_gtex.bim",

			 "../output/BLACK.james_cojo_ld_ref_ids_to_vander_gtex.tsv")
}
