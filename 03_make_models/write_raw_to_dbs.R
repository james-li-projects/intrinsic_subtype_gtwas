CELL_TYPES <- c("Stromal_and_Immune_cells", "Adipocytes", "Breast_tissue", "Epithelial_cells", "Endothelial_cells")
SUBTYPES <- c("HRNEG_HER2NEG", "HRNEG_HER2POS", "HRPOS_HER2NEG", "HRPOS_HER2POS")
TABLE_NAMES <- c("extra", "weights")
JM_TO_JL_MAPPING <- c("BLACK" = "afr", "WHITE" = "eur") # We referred to these groups differently in our files/code

main <- function(){
  library(tidyverse)
  library(glue)
  library(DBI)
  library(argparse)

  # Set up argument parser
  parser <- ArgumentParser()
  parser$add_argument("--ancestry", required=FALSE, default="BLACK", help="Ancestry to process")
  parser$add_argument("--ctype", required=FALSE, default="Breast_tissue",  help="Cell type to process")
  args <- parser$parse_args()
  ancestry <- args$ancestry
  ctype <- args$ctype

  wipe_dbs(ancestry, ctype)

  for (mtype in c("gene")){
    the_db <- glue("../output/predixcan_dbs/{ancestry}.{mtype}_models.{ctype}.db")
    con <- dbConnect(RSQLite::SQLite(), dbname = the_db)

    pat <- glue("{ancestry}\\.ENSG\\d{{1,}}\\.\\d{{1,}}\\.gene\\.{ctype}\\.extra")
    relevant_files <- list.files("../output/raw/", pat, full.names=T)
    write_file_list_to_table(relevant_files, con, the_db, ctype, ancestry, mtype)
  }
}


write_file_list_to_table <- function(file_list, con, db_name, cell_type, ancestry, mtype, test=FALSE){
	library(data.table)
	library(glue)
	library(DBI)

	cojo_ref <- fread(glue("../output/{ancestry}.vander_gtex_model_ids_to_james_cojo_ld_ref.tsv"))
	master_extra <- tibble()
	master_weights <- tibble()

	if (isTRUE(test)){
		# file_list <- file_list[1:2]
		ensembl_ids <- c(
										 "ENSG00000108784.9"
		)

		# Filter file_list to only those containing at least one exact match
		file_list <- file_list[str_detect(file_list, paste(ensembl_ids, collapse="|"))]
	}

	for (cur_path in file_list){
		gene_name <- str_extract(cur_path, "(ENSG\\d{1,}\\.\\d{1,})")

		if (str_detect(cur_path, "extra")){
			ins_table <- "extra"
		} else if (str_detect(cur_path, "weight")){
			ins_table <- "weights"
		} else {
			next # Don't insert model performance metrics
		}


		cur_df <- fread(cur_path)
		if (ins_table == "weights"){
			cur_df <- cur_df %>%
				filter(weight != 0)

		} else if (ins_table == "extra"){
			cur_df <- cur_df %>%
				filter(`n.snps.in.model` != 0) %>%
				filter(pred.perf.R2 > .01)
		}

		# Don't put something in the extra table that won't have any 
		# corresponding rows in the weights table and vice-versa
		if (nrow(cur_df) == 0){
			# message(glue("Removed {cur_path} (all 0 weights)"))
			# file.remove(cur_path)
			next
		} else {
			weight_df <- fread(glue("../output/raw/{ancestry}.{gene_name}.{mtype}.{cell_type}.weight.tsv")) %>%
				filter(weight != 0) %>%
				mutate(varID = str_replace(varID, "chr", ""))
			master_weights <- bind_rows(master_weights, weight_df)

			# append_res <- dbAppendTable(con, "weights", weight_df)
			# if (append_res >= 1 & append_res == nrow(weight_df)){
			# } else {
			# 	message(glue("Weight insert failed for {cur_path}."))
			# }
		}

		master_extra <- bind_rows(master_extra, cur_df)
		# append_res <- dbAppendTable(con, ins_table, cur_df)
		# # Check that the inserts worked
		# if (append_res >= 1 & append_res == nrow(cur_df)){
		# 	num_pass <- num_pass + 1
		# } else {
		# 	message(glue("Insert failed for {cur_path}."))
		# }
	}

	# Create a unique dataframe to avoid redundant joins
  master_weights_uniq_df <- master_weights %>%
    distinct(varID) %>%
    left_join(cojo_ref %>% select(varID=snp_bim1, cojo_snp_id=snp_bim2), by="varID")

  # Assign ancestry based ld region
  # Extract columns to compare to
  master_weights_uniq_df <- master_weights_uniq_df %>%
    mutate(
      chromosome = as.numeric(str_extract(varID, "^[^:]+")),
      snp_position = as.numeric(str_extract(varID, "(?<=:)[^:]+"))
    )

	# Has columns chr, start, stop
	# e.g.
	#      chr    start     stop
	#   <int>    <int>    <int>
	#1:    16 87636152 88487630
	ld_ances_name <- JM_TO_JL_MAPPING[ancestry]
  # ances_ld_df <- fread(glue("../../040_parse_through_initial_twas/input/ld_blocks/{ancestry}_ld.bed"))
  ances_ld_df <- fread(glue("../../040_parse_through_initial_twas/input/ld_blocks/ma_focus_lds/grch38.{ld_ances_name}.loci.bed"))

	# Assign ancestry-based LD region using a loop
	master_weights_uniq_df$ancestry.ld.start <- NA_real_
	master_weights_uniq_df$ancestry.ld.stop <- NA_real_

	for (i in seq_len(nrow(master_weights_uniq_df))) {
		chrom <- master_weights_uniq_df$chromosome[i]
		pos <- master_weights_uniq_df$snp_position[i]

		matched_rows <- ances_ld_df[ances_ld_df$chr == chrom & ances_ld_df$start <= pos & ances_ld_df$stop >= pos, ]
		
		if (nrow(matched_rows) > 0) {
			master_weights_uniq_df$ancestry.ld.start[i] <- matched_rows$start[1]
			master_weights_uniq_df$ancestry.ld.stop[i] <- matched_rows$stop[1]
		}
	}

	# Rename and clean up columns
	master_weights_uniq_df <- master_weights_uniq_df %>%
		rename(ancestry.ld.chr = chromosome) %>%
		select(-snp_position)

  # Join in columns that indicate whether the SNP is in the original GWAS
  # Useful later when determining which list of SNPs we want conditional effects for
  for (stype in SUBTYPES){
    col_name <- glue("in.GWAS.{stype}")
    gwas_path <- glue("../../040_parse_through_initial_twas/output/metaxcan_gwas_parse/COJO.{ancestry}.{stype}.txt.gz")
    gwas_df <- fread(gwas_path) %>%
      select(cojo_snp_id=panel_variant_id) %>%
      mutate(!!col_name := 1)

    master_weights_uniq_df <- master_weights_uniq_df %>%
      left_join(gwas_df, by = "cojo_snp_id") %>%
      mutate(!!col_name := replace_na(.[[col_name]], 0))
  }


  # Join the processed unique data back to the main dataframe
  master_weights <- master_weights %>%
    left_join(master_weights_uniq_df, by = c("varID"))

	# Print the difference in lengths
  print(glue("Length of master_weights: {nrow(master_weights)}"))
  print(glue("Length of master_weights_uniq_df: {nrow(master_weights_uniq_df)}"))


	if (length(unique(master_weights$gene)) != nrow(master_extra)){
		browser()
	} else {
		append_res_w <- dbAppendTable(con, "weights", master_weights)
		append_res_e <- dbAppendTable(con, "extra", master_extra)
		if (append_res_w >= 1 & append_res_w == nrow(master_weights)){
			message("Inserted weights")
		} else {
			message(glue("Insert failed for {cur_path}."))
		}

		if (append_res_e >= 1 & append_res_e == nrow(master_extra)){
			message("Inserted extra")
		} else {
			message(glue("Insert failed for {cur_path}."))
		}
	}

	message(glue("Finished writing {db_name}"))
  dbDisconnect(con)
}

	# Insert into DB (copypasta)
	# 	sqlite> SELECT * FROM extra LIMIT 5;
	# gene|genename|gene_type|n.snps.in.model|pred.perf.R2|pred.perf.pval|pred.perf.qval
	# ENSG00000000457.13|SCYL3|protein_coding|1|||
	# ENSG00000000460.16|C1orf112|protein_coding|1|||
	# ENSG00000001167.14|NFYA|protein_coding|1|||
	# ENSG00000001460.17|STPG1|protein_coding|2|||
	# ENSG00000001461.16|NIPAL3|protein_coding|2|||
	# sqlite> SELECT * FROM weights LIMIT 5;
	# gene|rsid|varID|ref_allele|eff_allele|weight
	# ENSG00000180549.7|rs4880192|chr9_137032610_A_G_b38|A|G|-0.137037138259738
	# ENSG00000107281.9|rs2891950|chr9_137045741_C_G_b38|C|G|0.134760315076842
	# ENSG00000054179.11|chr9_137054225_G_A_b38|chr9_137054225_G_A_b38|G|A|0.0563672816193278
	# ENSG00000054179.11|rs12551623|chr9_137054480_G_T_b38|G|T|0.120511488623684
	# ENSG00000186193.8|rs140069694|chr9_137065212_T_C_b38|T|C|0.412661833832297

wipe_dbs <- function(ances, ctype){
	library(DBI)
	library(glue)

	for (model_type in c("gene"))
		the_db <- glue("../output/predixcan_dbs/{ances}.{model_type}_models.{ctype}.db")
		con <- dbConnect(RSQLite::SQLite(), dbname = the_db)


		# Add columns to say whether the SNP is present in the GWAS.
		# (Also if it's present in GWAS it's present in James's COJO reference LD)
		for (stype in c("cojo_snp_id", "ancestry.ld.chr", "ancestry.ld.start",
										"ancestry.ld.stop", paste0("in.GWAS.", SUBTYPES))){
			table_info <- dbGetQuery(con, "PRAGMA table_info(weights);")
			column_exists <- stype %in% table_info$name
			if (!column_exists) {
				dbExecute(con, glue("ALTER TABLE weights ADD COLUMN '{stype}' INTEGER;"))  # Adjust the column type as needed
				print(glue("Column '{stype}' has been added to the 'weights' table.\n"))
			} else {
				print(glue("Column '{stype}' already exists in the 'weights' table.\n"))
			}
		}

		dbExecute(con, "DELETE FROM weights; VACUUM;")
		dbExecute(con, "DELETE FROM extra; VACUUM;")
		dbExecute(con, "VACUUM;")

		dbDisconnect(con)
		message(glue("{the_db} wiped"))
}

if (interactive()){
	main()
} else {
	main()
}
