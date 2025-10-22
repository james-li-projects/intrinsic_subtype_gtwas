arg_parser <- function(){
  require(argparse)
  parser <- ArgumentParser(description=paste0("",
                                              "",
                                              "",
                                              "",
                                              ""))

  parser$add_argument("--chr_num", type="integer", required=T, default=NULL,
                      help="LD block chromosome"
  )


  parser$add_argument("--group", type="character", required=T, default=NULL,
                      help="WHITE or BLACK"
  )
  parser$add_argument("--y_type", type="character", required=T, default=NULL,
                      help="gene or intron"
  )

  parser_args  <- parser$parse_args()
  return(parser_args)
}

main <- function(group="WHITE", y_target=NULL, y_type=NULL,

								 # Arguments to support running by chr (reduce I/O)
								 y_source,
								 y_target_info,
								 geno_obj,
								 prior_df,
								 cov_and_pc_df,
								 sample_ids.train,
								 sample_ids.test,
								 chr_num,

								 window_size=1000000){
	require(tidyverse)
	require(data.table)
	require(glue)
	require(genio)
	require(MiXcan)
	require(purrr)

	print(glue("{group} | {y_target} | {y_type}"))

	if (y_type == "gene"){
		mixcan_y <- y_source %>%
			filter(Name == y_target) %>%
			select(-Name, -genename) %>%
			t()
		mixcan_y <- mixcan_y[sample_ids.train$id_clean_gtex_per,] %>% as.double()
		y_target_info <- y_target_info %>%
			filter(gene_id == y_target)

		if (nrow(y_target_info) == 0){
			message(glue("No gencode 26 info for {y_target} | {group}"))
		}

		y_start <- y_target_info$start_location
		y_end <- y_target_info$end_location

		geno_obj$bim_train <- geno_obj$bim %>%
			filter(pos >= max(0, y_start - window_size) &
						 pos <= y_end + window_size)

		geno_obj$train_X <- geno_obj$X[geno_obj$bim_train$id, sample_ids.train$id_dirty] %>% t()

		if (dim(geno_obj$train_X)[1] != length(sample_ids.train$id_dirty)){
			geno_obj$train_X <- geno_obj$train_X %>% t()
		}

		if (nrow(geno_obj$bim_train) == 0){
			message(glue("{y_target} found no variants in the window."))
			return(tibble(y_target=y_target))
		}

		y_target_label <- y_target
	} else if (y_type == "intron"){
		mixcan_y <- y_source %>%
			filter(Name == y_target) %>%
			select(-Name, -gene) %>%
			distinct() %>% # We have intron id and gene, one intron id might be mapped to more than one gene, (with numerical info repeated)
			t()

		mixcan_y <- mixcan_y[sample_ids.train$id_clean_gtex_dash,] %>% as.double()
		y_target_info <- y_target_info %>%
			filter(id == y_target) %>%
			separate(id, into = c("intron_chr", "intron_start", "intron_end"), remove=F, sep=":", convert=T)

		y_chr <- y_target_info$intron_chr
		y_start <- y_target_info$intron_start
		y_end <- y_target_info$intron_end

		y_target_label <- glue("intron_{y_chr}_{y_start}_{y_end}")

		geno_obj$bim_train <- geno_obj$bim %>%
			filter(pos >= max(0, y_start - window_size) &
						 pos <= y_end + window_size)

		geno_obj$train_X <- geno_obj$X[geno_obj$bim_train$id, sample_ids.train$id_dirty] %>% t()
	}

	tibble(want_prior=c("Epithelial cells", "Endothelial cells", "Fibroblasts", "Adipocytes")) %>%
		pmap(~run_mixcan_for_single_prior(id_df=sample_ids.train, group=group, mixcan_y=mixcan_y, geno_obj=geno_obj, cov=cov_and_pc_df,
																					want_prior=.,
																					y_target=y_target,
																					y_target_label=y_target_label,
																					y_type=y_type, y_source=y_source,
																					y_target_info=y_target_info,
																					prior_df=prior_df,
																					id_df.test=sample_ids.test,
																					chr_num=chr_num
																					))
	return()
}

run_mixcan_for_single_prior <- function(id_df, group, mixcan_y, geno_obj, cov,
																				want_prior, y_target, y_target_label, y_type, y_source, y_target_info, prior_df,

																				id_df.test,
																				chr_num,

																				num_retries = 20,
																				sqlite_retry_ms = 172800000,
																				write_to_db=F,
																				nested_cv=T
																				){
	want_prior_db_name <- want_prior %>% str_replace(" ", "_")
	perf_info_fname <- glue("../output/raw/{group}.{y_target}.{y_type}.{want_prior_db_name}.model_perf.tsv")

	# DEBUG
	# If we have performance info already then we shouldn't have to retry it, skip to next
	breast_check_condition <- file.exists(glue("../output/raw/{group}.{y_target}.{y_type}.Breast_tissue.extra.tsv"))

	if (file.exists(perf_info_fname) & breast_check_condition){ # Can retry if there is no breast tissue thing for this
		message(glue("{perf_info_fname} and breast tissue file found already, skipping. . ."))
		return()
	}

	if (!breast_check_condition){
		message(glue("{group}.{y_target}.{y_type} with {want_prior} prior doesn't have breast tissue file so let's try to get a breast file"))
	}

	# 		write_tsv(extra_insert_tib, glue("../output/raw/{group}.{y_target}.{y_type}.{want_prior_db_name}.extra.tsv"))
	# 		write_tsv(weight_df, glue("../output/raw/{group}.{y_target}.{y_type}.{want_prior_db_name}.weight.tsv"))
	# 		write_tsv(extra_insert_tib, glue("../output/raw/{group}.{y_target}.{y_type}.Breast_tissue.extra.tsv"))
	# 		write_tsv(weight_df, glue("../output/raw/{group}.{y_target}.{y_type}.Breast_tissue.weight.tsv"))

	
	# Get prior
	scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
	prior_df <- scale_values(prior_df[c(want_prior),])

	# Need 0 < prior < 1
	prior_df[prior_df == 0] <- 0 + 1e-6
	prior_df[prior_df == 1] <- 1 - 1e-6

	tries <- 0
	num_tissue_test_X_all_0 <- 0
	num_no_pred <- 0
	num_tissue_0_weights <- 0
	while (tries < num_retries){
		# start.time <- Sys.time()
		MiXcan_result <- MiXcan(y=mixcan_y, x=geno_obj$train_X, 
														cov = cov %>%
															select(-id) %>%
															as.matrix(),
														pi= prior_df[id_df$`id_clean_gtex_per`],
														xNameMatrix = geno_obj$bim_train$id,
														target=y_target
														)
		# end.time <- Sys.time()
		# time.taken <- end.time - start.time
		# browser()
		mixcan_result_type <- MiXcan_result$type
		tries <- tries + 1

		# Can sometimes get NoPredictor results by chance, maybe insert a retry here?
		if (mixcan_result_type != "NoPredictor"){
			weights <- MiXcan_extract_weight_with_tiss_also(MiXcan_result)

			if (isTRUE(nested_cv)){
				boot_res <- list(weights=weights)
				tissue_test_rsq <- MiXcan_result$tiss_lambda_info$rsq_mean
				tissue_test_rho <- MiXcan_result$tiss_lambda_info$cor_rho_mean
				tissue_test_cor_pval <- MiXcan_result$tiss_lambda_info$cor_test_pval
				return_tib <- tibble(chr_num=chr_num, target=y_target, mixcan_result=mixcan_result_type,
														 prior_used=want_prior, num_attempts=tries,
														 num_tissue_test_X_all_0=num_tissue_test_X_all_0,
														 num_no_pred=num_no_pred,
														 num_tissue_0_weights=num_tissue_0_weights,
														 
														 enm_rsq=tissue_test_rsq,
														 enm_rho=tissue_test_rho,
														 enm_test_cor_pval=tissue_test_cor_pval
														 )
			} else {
				boot_res <- Mixcan_bootstrap_stats_in_test(group, weights, geno_obj$X, y_target, y_type, y_source, mixcan_result_type, id_df.test)
				tissue_test_rsq <- boot_res$performance %>%
					filter(model == "tissue") %>%
					pull(Rsquared)
				tissue_test_X_all_0 <- boot_res$performance %>%
					filter(model == "tissue") %>%
					pull(boot_X_all_0)

				if (isTRUE(tissue_test_X_all_0)){
					num_tissue_test_X_all_0 <- num_tissue_test_X_all_0 + 1
				}

				if (all(boot_res$weights$weight_tissue == 0)){
					num_tissue_0_weights <- num_tissue_0_weights + 1
				}
			}

			if (is.na(tissue_test_rsq)){ # Retry when we get a null model
				if (tries == num_retries){
					message(glue("\t\tNo Rsquared retrievable after max tries ({num_retries}) | prior={want_prior} | num_tissue_test_X_all_0 {num_tissue_test_X_all_0} | num_no_pred {num_no_pred} | 0 weight instances {num_tissue_0_weights}"))
					return_tib <- tibble(chr_num=chr_num, target=y_target, mixcan_result=mixcan_result_type, prior_used=want_prior, num_attempts=tries,
															 num_tissue_test_X_all_0=num_tissue_test_X_all_0, num_no_pred=num_no_pred, num_tissue_0_weights=num_tissue_0_weights)

					if (!file.exists(perf_info_fname)){
						write_tsv(return_tib, perf_info_fname)
					}
					return()
					# return(return_tib)
				}
				next
			}

			message(glue("\t\tMiXcan successful after {tries} attempts | prior={want_prior} | num_tissue_test_X_all_0 {num_tissue_test_X_all_0} | num_no_pred {num_no_pred} | 0 weight instances {num_tissue_0_weights}"))
			if (isTRUE(nested_cv)){
			} else {
				return_tib <- tibble(chr_num=chr_num, target=y_target, mixcan_result=mixcan_result_type, prior_used=want_prior, num_attempts=tries,
																 num_tissue_test_X_all_0=num_tissue_test_X_all_0, num_no_pred=num_no_pred, num_tissue_0_weights=num_tissue_0_weights)
			}

			tries <- as.double("inf")
		} else {
			num_no_pred <- num_no_pred + 1
		}

	}

	if (mixcan_result_type == "NoPredictor"){
		message(glue("\t\tNo predictor after max tries ({num_retries}) | prior={want_prior} | num_tissue_test_X_all_0 {num_tissue_test_X_all_0} | num_no_pred {num_no_pred} | 0 weight instances {num_tissue_0_weights}"))
		return_tib <- tibble(chr_num=chr_num, target=y_target, mixcan_result=mixcan_result_type, prior_used=want_prior, num_attempts=tries,
															 num_tissue_test_X_all_0=num_tissue_test_X_all_0, num_no_pred=num_no_pred, num_tissue_0_weights=num_tissue_0_weights)

		if (!file.exists(perf_info_fname)){
			write_tsv(return_tib, perf_info_fname)
		}
		return()
		# return(return_tib)
	}
	if (!file.exists(perf_info_fname)){
		write_tsv(return_tib, perf_info_fname)
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

	# Might have multiple tissue models, so put one in as soon as possible, if
	# one is there already, leave it.
	# if (want_prior == expected_last_want_prior)

	if (isTRUE(write_to_db)){
		tiss_db <- glue("../output/predixcan_dbs/{group}.{y_type}_models.Breast_tissue.db")
		out_db_con <- DBI::dbConnect(RSQLite::SQLite(), dbname = tiss_db)
		RSQLite::sqliteSetBusyHandler(out_db_con, sqlite_retry_ms)
		check_df <- DBI::dbReadTable(out_db_con, "extra") %>%
			filter(gene == y_target_label)
		DBI::dbDisconnect(out_db_con)

		check_condition <- (nrow(check_df) == 0)
	} else {
		check_condition <- !file.exists(glue("../output/raw/{group}.{y_target}.{y_type}.Breast_tissue.extra.tsv"))
	}

	if (check_condition){
		if ("weight_tissue" %in% names(boot_res$weights)){
		} else {
			if (mixcan_result_type == "NonSpecific"){
				boot_res$weights <- boot_res$weights %>%
					mutate(weight_tissue = weight_cell_1)
			} else {
				message("wut")
			}
		}

		n.snps.in.model <- boot_res$weights %>% dplyr::filter(weight_tissue != 0) %>% nrow()
		# n.snps.in.model <- sum(boot_res$weights$weight_tissue != 0) # %>% dplyr::filter(weight_tissue != 0) %>% nrow()
		extra_insert_tib <- tibble(gene=y_target_label,
															 genename=y_target_info$gene_name,
															 gene_type=y_target_info$gene_type,
															 n.snps.in.model=n.snps.in.model,
															 pred.perf.R2=tissue_test_rsq,
															 pred.perf.pval=tissue_test_cor_pval,
															 pred.perf.qval=NA)
		weight_df <- boot_res$weights %>%
								  mutate(gene=y_target_label, rsid=NA) %>%
									select(gene, rsid, varID=xNameMatrix, weight=weight_tissue) %>%
									inner_join(geno_obj$bim %>% select(varID=id, ref_allele=ref, eff_allele=alt), by="varID") %>%
									select(gene, rsid, varID, ref_allele, eff_allele, weight) %>%
									filter(weight != 0)

		if (isTRUE(write_to_db)){
			tiss_db <- glue("../output/predixcan_dbs/{group}.{y_type}_models.Breast_tissue.db")
			out_db_con <- DBI::dbConnect(RSQLite::SQLite(), dbname = tiss_db)
			RSQLite::sqliteSetBusyHandler(out_db_con, sqlite_retry_ms)
			append_extra_res <- DBI::dbAppendTable(out_db_con, "extra", extra_insert_tib)
			append_weight_res <- DBI::dbAppendTable(out_db_con, "weights", weight_df)
			# Check that the inserts worked
			if (append_extra_res >= 1 & append_weight_res == nrow(weight_df)){
				print(glue("\t\tInsert into tissue database worked for {y_target}. | prior={want_prior}"))
			} else {
				message(glue("One or both inserts failed for {y_target}."))
			}

			DBI::dbDisconnect(out_db_con)
		} else {

			# "Double write" performance info. Whichever corresponding cell prior
			# "model_perf" was used for Breast tissue will be replicated as a Breast
			# tissue "model_perf" file, this way when writing into the tissue db
			# later we can properly attribute the additional performance metrics
			# 
			# ^^^ This will avoid a case where the tissue model was created using
			# prior A, (of A,B,C), but we might not know which of A, B, C
			# "model_perf" files is appropriate without this quick double write
			if (n.snps.in.model > 0){
				tiss_perf_info_fname <- glue("../output/raw/{group}.{y_target}.{y_type}.Breast_tissue.model_perf.tsv")

				write_tsv(return_tib, tiss_perf_info_fname)
				write_tsv(extra_insert_tib, glue("../output/raw/{group}.{y_target}.{y_type}.Breast_tissue.extra.tsv"))
				write_tsv(weight_df, glue("../output/raw/{group}.{y_target}.{y_type}.Breast_tissue.weight.tsv"))
			}
		}
	}

	if (mixcan_result_type == "CellTypeSpecific"){
		extra_insert_tib <- tibble(
		  gene=y_target_label,
			genename=y_target_info$gene_name,
			gene_type=y_target_info$gene_type,
			n.snps.in.model=boot_res$weights %>%
				filter(weight_cell_1 != 0) %>%
				nrow(),
			pred.perf.R2=tissue_test_rsq,
			pred.perf.pval=tissue_test_cor_pval,
			pred.perf.qval=NA)
		weight_df <- boot_res$weights %>%
								  mutate(gene=y_target_label, rsid=NA) %>%
									select(gene, rsid, varID=xNameMatrix, weight=weight_cell_1) %>%
									inner_join(geno_obj$bim %>% select(varID=id, ref_allele=ref, eff_allele=alt), by="varID") %>%
									select(gene, rsid, varID, ref_allele, eff_allele, weight) %>%
									filter(weight != 0)

		if (isTRUE(write_to_db)){
			cell1_db <- glue("../output/predixcan_dbs/{group}.{y_type}_models.{want_prior_db_name}.db")
			out_db_con <- DBI::dbConnect(RSQLite::SQLite(), dbname = cell1_db)
			RSQLite::sqliteSetBusyHandler(out_db_con, sqlite_retry_ms)
			append_extra_res <- DBI::dbAppendTable(out_db_con, "extra", extra_insert_tib) 
			append_weight_res <- DBI::dbAppendTable(out_db_con, "weights", weight_df)
			# Check that the inserts worked
			if (append_extra_res >= 1 & append_weight_res == nrow(weight_df)){
				print(glue("\t\tInsert into {want_prior} database worked for {y_target}. | prior={want_prior}"))
			} else {
				message(glue("One or both inserts failed for {y_target}."))
			}
			DBI::dbDisconnect(out_db_con)
		} else {
			extra_name <- glue("../output/raw/{group}.{y_target}.{y_type}.{want_prior_db_name}.extra.tsv")
			weight_name <- glue("../output/raw/{group}.{y_target}.{y_type}.{want_prior_db_name}.weight.tsv")

			if (!file.exists(extra_name)){
				write_tsv(extra_insert_tib, extra_name)
			} 

			if (!file.exists(weight_name)){
				write_tsv(weight_df, weight_name)
			}
		}


		if (want_prior == "Epithelial cells"){ # In this case we can extract cell type 2 as "Stromal_and_Immune_cells"
			extra_insert_tib <- tibble( gene=y_target_label,
																 genename=y_target_info$gene_name,
																 gene_type=y_target_info$gene_type,
																 n.snps.in.model=boot_res$weights %>%
																	 filter(weight_cell_2 != 0) %>% nrow(),
																pred.perf.R2=tissue_test_rsq,
																pred.perf.pval=tissue_test_cor_pval,
																 pred.perf.qval=NA)
			weight_df <- boot_res$weights %>%
										mutate(gene=y_target_label, rsid=NA) %>%
										select(gene, rsid, varID=xNameMatrix, weight=weight_cell_2) %>%
										inner_join(geno_obj$bim %>% select(varID=id, ref_allele=ref, eff_allele=alt), by="varID") %>%
										select(gene, rsid, varID, ref_allele, eff_allele, weight) %>%
										filter(weight != 0)
			if (isTRUE(write_to_db)){
				cellOther_db <- glue("../output/predixcan_dbs/{group}.{y_type}_models.Stromal_and_Immune_cells.db")
				out_db_con <- DBI::dbConnect(RSQLite::SQLite(), dbname = cellOther_db)
				RSQLite::sqliteSetBusyHandler(out_db_con, sqlite_retry_ms)
				append_extra_res <- DBI::dbAppendTable(out_db_con, "extra", extra_insert_tib)
				append_weight_res <- DBI::dbAppendTable(out_db_con, "weights", weight_df)
				# Check that the inserts worked
				if (append_extra_res >= 1 & append_weight_res == nrow(weight_df)){
					print(glue("\t\tInsert into Stromal_and_Immune_cells database worked for {y_target}. | prior={want_prior}"))
				} else {
					message(glue("One or both inserts failed for {y_target}."))
				}
				DBI::dbDisconnect(out_db_con)
			} else {
				extra_name <- glue("../output/raw/{group}.{y_target}.{y_type}.Stromal_and_Immune_cells.extra.tsv")
				weight_name <- glue("../output/raw/{group}.{y_target}.{y_type}.Stromal_and_Immune_cells.weight.tsv")

				if (!file.exists(extra_name)){
					write_tsv(extra_insert_tib, extra_name)
				}

				if (!file.exists(weight_name)){
					write_tsv(weight_df, weight_name)
				}
			}

		}
	} else if (mixcan_result_type == "NonSpecific"){
		message(glue("\t\tNo cell type db insertions for {y_target}. Model is NonSpecific | prior={want_prior}"))
	}


	# return_tib$rsquared <- boot_res$performance %>%
	# 					filter(model == "tissue") %>%
	# 					pull(Rsquared) %>%
	# 					mean()
	# return(return_tib)
}

get_resample_df <- function(results_obj, outer_or_inner, analysis_or_assess, outer_num=1,
														inner_num=1){
	if (outer_or_inner == "outer"){
		return(results_obj$splits[[outer_num]] %>% analysis_or_assess())
	} else if (outer_or_inner == "inner"){
		cur_obj <- results_obj %>%
			filter(id == glue("Fold{outer_num}"))
		return(cur_obj$inner_resamples[[1]]$splits[[inner_num]] %>% analysis_or_assess())
	}

}

double_check_splits <- function(results_obj){
	for (outer_id in 1:5){
		outer_test <- get_resample_df(results_obj, "outer", assessment, outer_num=outer_id)
		for (inner_id in 1:10){
			inner_train <- get_resample_df(results_obj, "inner", analysis, outer_id, inner_id)
			inner_test <- get_resample_df(results_obj, "inner", assessment, outer_id, inner_id)

			# The test set of outer loop should not overlap with any part of the inner loop
			outer_inner_intersect <- intersect(c(inner_train$id_col, inner_test$id_col), outer_test$id_col)
			if (length(outer_inner_intersect) >= 1){
				print(outer_inner_intersect)
			}
		}
	}
}

get_lambda <- function(x, y, xNameMatrix, tiss=F, n_train_test_folds=5, n_k_folds=10, target=NULL){
	require(tidymodels)
	colnames(y) <- "y"

	train_dat <- bind_cols(y, x,
												 .name_repair="unique_quiet") %>% # STFU cell model colname
	  mutate(id_col = row_number())

	results <- nested_cv(train_dat,
                     outside = vfold_cv(v=n_train_test_folds),
                     inside = vfold_cv(v=n_k_folds))

	require(furrr)
	plan(multicore) # set --ntasks=1 and --cpus-per-task=X sensibly: this plan will use X

	# Uncomment to check/debug inner/outer loops
	# inner_best_df <- map(results$inner_resamples, summarize_tune_results, perf_metrics=F) %>%
	# 	bind_rows() %>%
	# 	rename(lambda=penalty)
	inner_best_df <- future_map(results$inner_resamples, summarize_tune_results, perf_metrics=F, .options = furrr_options(seed = NULL)) %>%
		bind_rows() %>%
		rename(lambda=penalty)

	outer_loop_input <- results %>%
		bind_cols(inner_best_df)

	outer_loop_output <- outer_loop_input %>%
		future_pmap_dfr(~glmnet_perf_metrics(..1, ..4, perf_metrics=tiss), .options = furrr_options(seed = NULL))
	# Uncomment to check/debug inner/outer loops
	# outer_loop_output <- outer_loop_input %>%
	# 	pmap_dfr(~glmnet_perf_metrics(..1, ..4, perf_metrics=T))

	# https://rdrr.io/cran/nestedcv/src/R/nestedcv.R
	# Following the above, we return the median lambda value, and output the final mean metrics
	if (isTRUE(tiss)){
		return_list <- list(lambda=median(outer_loop_output$penalty),
												cor_rho_mean = mean(outer_loop_output$rho, na.rm=T),
												cor_test_pval=pchisq(-2 * sum(log(outer_loop_output$cor.test.p)), 2 * n_train_test_folds, lower.tail=F),
												rsq_mean = mean(outer_loop_output$rsq, na.rm=T))
	} else {
		# If it's a cell model we don't got no metrics anyway cause they would suck
		return_list <- list(lambda=median(outer_loop_output$penalty))
	}
	return(return_list)
}

glmnet_perf_metrics <- function(object, lambda=NULL, perf_metrics=F, ...){
	y_col <- ncol(object$data)

	if (isTRUE(perf_metrics)){
		penalty.factor <- rep(1, ncol(analysis(object)) - 1)
	} else {
		penalty.factor <- c(0, rep(1, ncol(analysis(object)) - 2))
	}

	if (is.null(lambda)){
		mod <- linear_reg(mixture=0.5, penalty=tune()) %>%
			set_engine("glmnet", penalty.factor=penalty.factor) %>%
			fit(y ~ ., data=analysis(object))
		pred_penalty <- mod$fit$lambda
	} else {
		mod <- linear_reg(mixture=0.5, penalty=lambda) %>%
			set_engine("glmnet", penalty.factor=penalty.factor) %>%
			fit(y ~ ., data=analysis(object))
		pred_penalty <- lambda
	}

	if (!is.null(lambda)){
		# If we're specifying lambda, we're trying to get perf metrics, so we only
		# calculate predictions based on SNP weights
		coefs <- tidy(mod) %>%
			filter(estimate != 0) %>%
			filter(str_detect(term, "chr")) %>%
			mutate(term = str_replace_all(term, "`", ""))
		perf_pred_x <- assessment(object) %>%
			select(all_of(coefs$term))

		coef_matrix <- as.matrix(coefs$estimate) %>% t()
		perf_pred_x_matrix <- perf_pred_x %>% as.matrix() %>% t()
		pred_vec <- coef_matrix%*%perf_pred_x_matrix %>% as.vector()

		preds <- tibble(.pred=pred_vec) %>%
			bind_cols(assessment(object) %>% dplyr::select(y)) %>%
			add_rowindex() %>%
			mutate(penalty=lambda)
	} else {
		preds <- multi_predict(mod, assessment(object) %>% dplyr::select(-y), penalty=pred_penalty) %>% 
			bind_cols(assessment(object) %>% dplyr::select(y)) %>%
			add_rowindex() %>%
			unnest(cols = ".pred")
	}

	mse_vec <- preds %>%
		group_by(penalty) %>%
		rmse(y, .pred) %>%
		mutate(mse = .estimate ^ 2)

	# Matches some of the elastic net criteria in https://github.com/hakyimlab/summary-gwas-imputation/blob/206dac587824a6f207e137ce8c2d7b15d81d5869/src/elastic_net.R#L24
	if (isTRUE(perf_metrics)){
		enm_vec <- preds %>% 
			group_by(penalty) %>%
			group_map(~enm_metrics_func(.x), .keep=T) %>%
			bind_rows()
	} else {
		enm_vec <- tibble(penalty=NA)
	}


	return_tib <- tibble(bind_cols(mse_vec %>% select(penalty, mse), enm_vec %>% select(-penalty)))

	return(return_tib)
}

calc_R2 <- function(y, y_pred) {
	tss <- sum(y**2)
	rss <- sum((y - y_pred)**2)
	1 - rss/tss
}

enm_metrics_func <- function(data, na.rm=T){
  # perf_estimates <- caret::postResample(data$.pred, data$y )
	rho <- ifelse(((sd(data$.pred) != 0) & (sd(data$y) != 0)), cor(data$.pred, data$y), 0)
  cor.test.p <- ifelse((sd(data$.pred) != 0) & (sd(data$y) != 0), cor.test(data$.pred, data$y)$p.value, runif(1)) %>% unlist()
  rsq <- rsq(data, y, .pred)$.estimate

	return(tibble(penalty=unique(data$penalty), rho=rho, cor.test.p=cor.test.p, rsq=rsq))
}

tune_over_perf <- function(object, perf_metrics){
	metrics <- glmnet_perf_metrics(object=object, perf_metrics=perf_metrics)
	return(metrics)
}

summarize_tune_results <- function(object, perf_metrics=F){
	one_inner_loop_res <- map_dfr(object$splits, tune_over_perf, perf_metrics=perf_metrics)

	cv.lambda.min.df <- one_inner_loop_res %>%
		filter(mse == min(mse, na.rm=T))

	# Lambda min or lambda 1se?
	mse_se <- sd(one_inner_loop_res$mse) / 10 # num_folds
	cv.lambda.1se.df <- one_inner_loop_res %>%
		filter(mse <= min(mse, na.rm=T) + mse_se) %>%
		arrange(desc(penalty)) %>% head(1)

	return_df <- cv.lambda.1se.df
	
	if (nrow(return_df) > 1){
		return_tib <- return_df %>%
			filter(penalty = max(penalty, na.rm=T))
	} else {
		return_tib <- return_df
	}

	return(return_tib)
}

MiXcan=function(y, x, cov=NULL, pi, xNameMatrix=NULL, target, yName=NULL,
                foldid=NULL) {
  # create new predictors
  x=as.matrix(x); y=as.matrix(y)
  n=nrow(x); p=ncol(x)

  if(is.null(xNameMatrix)) {xNameMatrix=paste0("SNP", 1:p)}
  if (is.null(foldid)) {foldid= sample(1:10, n, replace=T)}

  if (is.null(cov)) {
    pcov=0; xcov=x;
    ci=pi-0.5; z=ci*x; z=as.matrix(z); xx=as.matrix(cbind(ci, x, z))
  }

  if (is.null(cov)==F) {
    cov=as.matrix(cov)
    pcov=ncol(cov); xcov=as.matrix(cbind(x, cov))
    ci=pi-0.5; z=ci*x;
    xx=as.matrix(cbind(ci, x, z, cov))
  }

  # tissue model
	tiss_lambda_info <- get_lambda(xcov,y, xNameMatrix, tiss=T, target=target)
  #ft00=glmnet::cv.glmnet(x=xcov, y=y,family="gaussian",  foldid=foldid, alpha=0.5)
  ft0=glmnet::glmnet(x=xcov, y=y,  family="gaussian", lambda = tiss_lambda_info$lambda, alpha=0.5)
  est.tissue=c(ft0$a0,as.numeric(ft0$beta))

  # cell type specific model
	cell_lambda_info <- get_lambda(xx, y, xNameMatrix, tiss=F, target=target)
  # ft11=glmnet::cv.glmnet(x=xx, y=y,
  #                        penalty.factor=c(0, rep(1, ncol(xx)-1)),
  #                        family="gaussian", foldid=foldid, alpha=0.5)

  ft=glmnet::glmnet(x=xx, y=y, penalty.factor=c(0, rep(1, ncol(xx)-1)), family="gaussian",
                    lambda = cell_lambda_info$lambda, alpha=0.5)
  est=c(ft$a0,as.numeric(ft$beta))
  beta10=est[1]+est[2]/2
  beta20=est[1]-est[2]/2
  beta11=est[3: (p+2)] + est[(p+3): (2*p+2)]/2
  beta21=est[3: (p+2)] - est[(p+3): (2*p+2)]/2


  ## add inference for difference > 0
  Type ="NonSpecific"
  idx.diff=((p+3): (2*p+2)) -1
  idx.nonzero=which(est[-1]!=0)
  idx.nonzero.diff=intersect(idx.diff, idx.nonzero)

  if (length(idx.nonzero.diff)!=0) {
    xx.select=xx[,idx.nonzero]
    beta.ols.boot=NULL; # beta.en.boot=NULL
    for(boot in 1:200) {
      id=sample(n, n, replace =T)
      gfit = lm(y[id,]~xx.select[id,])
      beta.ols.boot =rbind(beta.ols.boot, coef(gfit)[-1])
    }
    beta.range=apply(beta.ols.boot, 2, function(f) quantile(f, prob=c(0.025, 0.975), na.rm=T))
    beta.diff.range=beta.range[, match(idx.nonzero.diff, idx.nonzero)]
    # print(beta.diff.range)
    if (is.null(dim(beta.diff.range))) {
      any.nonzero= beta.diff.range[1] * beta.diff.range[2]>0} else {
        print(apply(beta.diff.range, 2, function(f) f[1] * f[2]>0))
        any.nonzero= any(apply(beta.diff.range, 2, function(f) f[1] * f[2]>0), na.rm=T)
      }

    if (is.na(any.nonzero)==F & any.nonzero==T) {Type ="CellTypeSpecific"}
  }


  if (Type!="CellTypeSpecific") {
    beta1=beta2=est.tissue
  } else {
    if (is.null(cov)) {
      beta1=c(beta10, beta11)
      beta2=c(beta20, beta21)
    }
    if (is.null(cov)==F) {
      beta_cov=est[ (2*p+3): (2*p+2+pcov)]
      beta1=c(beta10, beta11, beta_cov)
      beta2=c(beta20, beta21, beta_cov)
    }
  }
  beta.all.models=cbind(est.tissue, beta1, beta2)
  colnames(beta.all.models)=c("Tissue", "Cell1", "Cell2")
  beta.SNP.cell1=data.frame(xNameMatrix, weight=beta1[2:(p+1)])
  beta.SNP.cell2=data.frame(xNameMatrix, weight=beta2[2:(p+1)])


  if (suppressWarnings(
    all(c(beta.SNP.cell1$weight, beta.SNP.cell2$weight)==0) )) {
    Type ="NoPredictor"}


  return(list(type=Type,
              beta.SNP.cell1=beta.SNP.cell1,
              beta.SNP.cell2=beta.SNP.cell2,
              beta.all.models=beta.all.models,
              glmnet.cell=ft,
              glmnet.tissue=ft0,
              yName=yName,
              xNameMatrix=xNameMatrix,
							tiss_lambda_info=tiss_lambda_info,
							cell_lambda_info=cell_lambda_info
							))
}

get_cov_and_pcs_df <- function(group, y_type, ids, num_geno_PC=5){
	vander_cov <- fread("../input/Komen.cov.tsv") %>%
		filter(Barcode %in% ids) %>%
		select(id=Barcode, Age)
	gtex_cov <- fread("../input/gtex_cov.tsv") %>%
		filter(SUBJID %in% ids) %>%
		select(id=SUBJID, Age=AGE)


	genoPC_col_select <- c("id", paste0("PC", 1:num_geno_PC))
	genoPC_df <- fread(glue("../input/pcs/{group}.vander_gtex_chrALL_pruned_genotypes_PCs.eigenvec")) %>%
		select(-`#FID`) %>%
		rename(id=IID) %>%
		mutate(id = case_when(str_detect(id, "Komen") ~ str_replace(id, "Komen-BC-", ""),
													TRUE ~ id)) %>%
		select(all_of(genoPC_col_select))

	names(genoPC_df) <- str_replace(names(genoPC_df), "PC", "genoPC")
	if (y_type == "gene"){
		pc_name <- "expression"
		featurePC_df <- fread(glue("../input/pcs/{group}.expression_vander_gtex_top_pcs.tsv")) %>%
			mutate(id = case_when(str_detect(id, "GTEX") ~ str_replace(id, "\\.", "-"),
													str_detect(id, "Komen") ~ str_replace(id, "Komen-BC-", ""),
													TRUE ~ id))
	} else if (y_type == "intron"){
		featurePC_df <- fread(glue("../input/pcs/{group}.splicing_vander_gtex_top_pcs.tsv"))
		pc_name <- "intron"
	}
	names(featurePC_df) <- str_replace(names(featurePC_df), "PC", glue("{pc_name}PC"))

	return_cov <- bind_rows(vander_cov, gtex_cov) %>%
		arrange(match(id, ids)) %>%
		left_join(genoPC_df, by="id") %>%
		left_join(featurePC_df, by="id")

	rownames(return_cov) <- ids
	return(return_cov)
}

Mixcan_bootstrap_stats_in_test <- function(group, weights, full_X, y_target, y_type, y_source, mixcan_res_type, id_df.test, num_reps=100,
																					 screw_boot=T){
	if (isTRUE(screw_boot)){
		num_reps <- 1
	}
	# The IDs here really kind of suck
	do_one_boot_run <- function(bait){

		if (isTRUE(screw_boot)){
		  boot_ids <- id_df.test$`id_clean_gtex_per`
		} else {
		  boot_ids <- sample(id_df.test$`id_clean_gtex_per`, replace=T, size=length(id_df.test$`id_clean_gtex_per`))
		}
		# print(rlang::hash(boot_ids))
		if (y_type == "gene"){
			mixcan_y <- y_source %>%
				filter(Name == y_target) %>%
				t()

			mixcan_y <- mixcan_y[boot_ids,] %>% as.double()
		} else if (y_type == "intron"){
			mixcan_y <- y_source %>%
				filter(Name == y_target) %>%
				select(-all_of(c("Name", "gene"))) %>%
				distinct() %>%
				t()

			mixcan_y <- mixcan_y[boot_ids %>% str_replace("GTEX\\.", "GTEX-"),] %>% as.vector()
		}
	  boot_X <- full_X[weights$xNameMatrix, str_replace(str_replace(boot_ids, "K1", "Komen-BC-K1"), "GTEX\\.", "GTEX-")] %>% t()

		if (all(boot_X == 0)){
			boot_X_all_0 <- TRUE
		} else {
			boot_X_all_0 <- FALSE
		}

		pred_vals <- MiXcan_prediction_with_tiss_also(weights, boot_X)

		# If the MiXcan result is cell type specific, then there is a difference
		# between the weights, if not, there is no difference, so we clean things
		# up a bit
		pred_vals <- as_tibble(pred_vals)
		if (mixcan_res_type == "CellTypeSpecific"){
			new_pred_vals <- pred_vals
		} else if (mixcan_res_type == "NonSpecific"){
			colnames(pred_vals) <- c("tissue", "throwaway")

			new_pred_vals <- pred_vals[,c("tissue")]
		}

		return_df <- tibble()
		for (cname in colnames(new_pred_vals)){
		  perf_estimates <- caret::postResample(new_pred_vals[,cname], mixcan_y)
			return_df <- bind_rows(return_df, perf_estimates %>% as.list() %>% as_tibble() %>% mutate(model=cname))
		}
		return(return_df %>% mutate(boot_X_all_0 = boot_X_all_0))
	}

	boot_results <- 1:num_reps %>%
		map_dfr(~do_one_boot_run(.))
	if (mixcan_res_type == "CellTypeSpecific"){
		new_weights <- weights
	} else if (mixcan_res_type == "NonSpecific"){
		new_weights <- weights[,c("xNameMatrix", "weight_cell_1")]
		names(new_weights) <- c("xNameMatrix", "weight_tissue")
		boot_results <- boot_results %>%
			distinct()
	}

	return(list(performance=boot_results, weights=new_weights))
}

MiXcan_prediction_with_tiss_also <- function(weights, new_x){
	# Need this, sometimes we get single SNP models
	if (dim(weights)[1] == 1){
		new_xx <- t(new_x)
	} else {
		new_xx <- new_x
	}

  yhat_MiXcan_cell_1 <- new_xx %*% as.matrix(weights[,"weight_cell_1"])
  yhat_MiXcan_cell_2 <- new_xx %*% as.matrix(weights[,"weight_cell_2"])
  yhat_MiXcan_prediction <- cbind(yhat_MiXcan_cell_1, yhat_MiXcan_cell_2)
  colnames(yhat_MiXcan_prediction) <- c("cell_1", "cell_2")

	if ("weight_tissue" %in% names(weights)){
		yhat_MiXcan_tissue <- new_xx %*% as.matrix(weights[,"weight_tissue"])
		yhat_MiXcan_prediction <- cbind(yhat_MiXcan_prediction, yhat_MiXcan_tissue)
    colnames(yhat_MiXcan_prediction) <- c("cell_1", "cell_2", "tissue")
	}
  return(yhat_MiXcan_prediction)
}

MiXcan_extract_weight_with_tiss_also <- function(model, keepZeroWeight=F){
	require(MiXcan)
	if (model$type == "CellTypeSpecific"){
		initial_df <- MiXcan_extract_weight(model, keepZeroWeight=T)

		beta.SNP.tissue <- data.frame(xNameMatrix=model$xNameMatrix, weight=model$beta.all.models[,1][2:(length(model$xNameMatrix)+1)])

		result_weight_tissue <-
			beta.SNP.tissue %>%
			dplyr::rename(weight_tissue = weight)

		combined_df <- initial_df %>%
			left_join(result_weight_tissue, by = "xNameMatrix")

	  if (keepZeroWeight==F) {
      combined_df <- combined_df %>%
        dplyr::filter(!(weight_tissue == 0 & weight_cell_1 == 0 & weight_cell_2 == 0))
    }
	} else {
		return(MiXcan_extract_weight(model, keepZeroWeight=keepZeroWeight))
	}
}

make_test_lists <- function(group="WHITE"){
	require(tidyverse)
	require(data.table)
	require(glue)
	require(purrr)
	gcode_26_df <- fread("../input/gencode_v26_all.txt") %>%
		mutate(chr_num = str_extract(chromosome, "\\d{1,}"))

	group_df <- fread(glue("../input/model_target_lists/{group}.full_non-repeat_gene_list.txt"))

	got_data_df <- gcode_26_df %>%
		inner_join(group_df, by = c("gene_id" = "id")) %>%
		filter(chr_num %in% 1:22)

	intron_df <- fread(glue("../input/model_target_lists/{group}.full_non-repeat_intron_list.txt")) %>%
		mutate(chr = str_extract(id, "^\\d{1,2}") %>% as.integer()) %>%
		inner_join(gcode_26_df, by = c("gene" = "gene_id"))

	for (cur_chr in 1:22){
		write_df <- got_data_df %>%
			filter(chr_num == cur_chr) %>%
			select(id=gene_id)

		write_intron_df <- intron_df %>%
			filter(chr == cur_chr)

		write_tsv(write_df, glue("../input/model_target_lists/by_chr/{group}.chr{cur_chr}_gene_list.txt"))
		write_tsv(write_intron_df, glue("../input/model_target_lists/by_chr/{group}.chr{cur_chr}_intron_list.txt"))
	}
}

do_by_chr <- function(group="WHITE", chr_num, y_type){
	require(genio)
	require(tidyverse)
	require(data.table)
	require(glue)

	targets_df <- fread(glue("../input/model_target_lists/by_chr/{group}.chr{chr_num}_{y_type}_list.txt"))
	sample_ids.train <- fread(glue("../input/splits/{group}.all.tsv"), col.names=c("FID", "id_dirty", "group")) %>%
		select(id_dirty) %>%
	  mutate(
					 id_clean_gtex_per = case_when(str_detect(id_dirty, "GTEX") ~ str_replace(id_dirty, "-", "\\."),
													str_detect(id_dirty, "Komen") ~ str_replace(id_dirty, "Komen-BC-", ""),
													TRUE ~ id_dirty),
					 id_clean_gtex_dash = case_when(str_detect(id_dirty, "Komen") ~ str_replace(id_dirty, "Komen-BC-", ""),
													TRUE ~ id_dirty),
					 ) %>%
		mutate(id=id_clean_gtex_dash)
	sample_ids.test <- fread(glue("../input/splits/{group}.test.tsv"), col.names=c("FID", "id_dirty", "group")) %>%
		select(id_dirty) %>%
	  mutate(
					 id_clean_gtex_per = case_when(str_detect(id_dirty, "GTEX") ~ str_replace(id_dirty, "-", "\\."),
													str_detect(id_dirty, "Komen") ~ str_replace(id_dirty, "Komen-BC-", ""),
													TRUE ~ id_dirty),
					 id_clean_gtex_dash = case_when(str_detect(id_dirty, "Komen") ~ str_replace(id_dirty, "Komen-BC-", ""),
													TRUE ~ id_dirty),
					 ) %>%
		mutate(id=id_clean_gtex_dash)

	cov_and_pc_df <- get_cov_and_pcs_df(group, y_type, sample_ids.train$id)

	if (y_type == "gene"){
		# exprMatrix <- read.table("../input/Komen.RNAseQC.TPM.GENCODEv26.hg38.txt", header=TRUE, row.names=1, as.is=TRUE)
		exprMatrix <- fread(glue("../input/{group}.Komen.GTEx.Breast.RBINT.TMM.TMP.Expr.GENCODEv26.hg38.txt"))
		gcode_26_df <- fread("../input/gencode_v26_all.txt")

		y_source <- exprMatrix %>%
			filter(Name %in% targets_df$id)

		y_target_info <- gcode_26_df %>%
			filter(gene_id %in% targets_df$id)


		y_start <- y_target_info$start_location
		y_end <- y_target_info$end_location

		y_source <- exprMatrix
	} else if (y_type == "intron"){
		# intron_info <- strsplit(y_target, ":")

		# y_chr <- as.integer(intron_info[[1]][1])
		# y_start <- intron_info[[1]][2] %>% as.double()
		# y_end <- intron_info[[1]][3] %>% as.double()
		# y_cluster <- intron_info[[1]][4] %>%
		# 	str_replace_all("_", "")

		# y_target_label <- glue("intron_{y_chr}_{y_start}_{y_end}_{y_cluster}")

		intronMatrix <- fread(glue("../input/{group}.Komen.GTEx.Breast.leafcutter.Introns.txt"))

		# mixcan_y <- intronMatrix %>%
		# 	filter(ID == y_target) %>%
		# 	select(-all_of(c("#Chr", "start", "end", "ID"))) %>%
		# 	t()
		# mixcan_y <- mixcan_y[paste0(sample_ids.train$id, ".markdup.sorted.bam"),] %>% as.vector()
		y_source <- intronMatrix
	}

	y_chr <- chr_num
	geno_obj <- read_plink(glue("../input/vander_gtex_genotypes/{group}.chr{y_chr}_vander_gtex"))

	# Some missing data here, do median imputation
	indx <- which(is.na(geno_obj$X), arr.ind = TRUE)
	geno_obj$X[indx] <- matrixStats::colMedians(geno_obj$X, na.rm = TRUE)[indx[, 2]]

	prior_df <- readRDS(glue("../input/xCell/{group}.Komen.GTEx.Breast.cell_scores.rds"))

	#write_df <- tibble()
	# Introns might be repeated, but we want to put the gene they go with in the
	# extra table
	for (cur_target in unique(targets_df$id)){ 
		# DEBUG
		# if (!(cur_target %in% c("ENSG00000223508.5","ENSG00000144619.14", "ENSG00000144619.14"))){
		# 	next
		# }
		if (y_type == "gene"){
		} else if (y_type == "intron"){
			y_target_info <- targets_df
		}
		main(group=group, y_target=cur_target, y_type=y_type,
				 y_source=y_source, y_target_info=y_target_info, geno_obj=geno_obj, prior_df=prior_df,
				 cov_and_pc_df=cov_and_pc_df, sample_ids.train=sample_ids.train,
				 sample_ids.test=sample_ids.test,
				 chr_num=chr_num
				 )
		# write_df <- bind_rows(write_df, cur_target_df)
	}

	# write_tsv(write_df, glue("../output/info/{group}.chr{chr_num}_{y_type}_info.tsv"))
}

if (interactive()){
	options(expressions = 5e5)
	do_by_chr(group="WHITE", chr_num=22, y_type="gene")
	# make_test_lists()
	# make_test_lists("BLACK")
	# do_by_chr(chr_num=22, y_type="intron")
	# main(y_target="ENSG00000000971.15", y_type="gene")
	# main(y_target="1:16765:16854:clu_7615_-", y_type="intron")
} else {
	options(error = function() traceback(2))
	options(expressions = 5e5)
	sink(stdout(), type="message")
	cli_args <- arg_parser()
	print(cli_args)
	do.call(do_by_chr, cli_args)
}
