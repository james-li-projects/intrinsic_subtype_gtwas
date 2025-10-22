main <- function(group="WHITE", k_fold=8){
	require(tidyverse)
	require(data.table)
	require(glue)
	require(R.utils)

	target_df <- fread(glue("../input/model_target_lists/{group}.full_non-repeat_gene_list.txt"))
	gcode_df <- fread("../input/gencode_v26_all.txt") %>%
		select(gene_id, gene_name, start_location, end_location, gene_type)

	info_df <- tibble()
	# WHITE.ENSG00000172215.5.gene.Fibroblasts.model_perf.tsv
	relevant_files <- list.files("../output/raw/",
															 pattern=glue("{group}.*model_perf.tsv"),
															 full.names=T)
	for (cur_file in relevant_files){
		cur_df <- fread(cur_file)

		weight_fname <- str_replace(cur_file, "model_perf\\.tsv", "weight\\.tsv")
		if (file.exists(weight_fname)){
			num_weights <- countLines(weight_fname)[[1]] - 1
			cur_df <- cur_df %>%
				mutate(num_weights = num_weights)
		}
		if (str_detect(cur_file, "Breast_tissue")){
			cur_df <- cur_df %>%
				mutate(breast_tiss = T)
		} else {
			cur_df <- cur_df %>%
				mutate(breast_tiss = F)
		}
		
		info_df <- bind_rows(info_df, cur_df)
	}

	out_df <- target_df %>%
		left_join(gcode_df, by = c("id" = "gene_id")) %>%
		left_join(info_df, by = c("id" = "target"))

	write_tsv(out_df, glue("../output/info_gwide/{group}.gene_gwide_info.tsv"))
}

if (interactive()){
	main()
	main("BLACK")
} else {
	main()
	main("BLACK")
}
