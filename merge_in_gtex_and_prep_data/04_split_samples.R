main <- function(){
	split_samples()
}

split_samples <- function(){
	require(tidyverse)
	require(data.table)
	require(glue)

		# select(id = V1) %>%
		# mutate(id = str_extract(id, "K\\d{1,}"))
	vander_fam <- fread("../output/vander_rsq_filtered_genotypes/chr1_rsq0.7.fam") %>%
		select(fid = V1, id = V2) %>%
		mutate(old_id = id) %>%
		mutate(id = case_when(str_detect(id, "GTEX") ~ id,
													str_detect(id, "Komen") ~ str_extract(id, "K\\d{1,}")))
	gtex_fam <- fread("../input/gtex_v8_whole_chr_bbf/processed_gtex.fam") %>%
		select(fid = V1, id = V2) %>%
		mutate(fid = as.character(fid)) %>%
		mutate(old_id = id) %>%
		mutate(id = case_when(str_detect(id, "GTEX") ~ id,
													str_detect(id, "Komen") ~ str_extract(id, "K\\d{1,}")))
	eg_fam_df <- bind_rows(vander_fam, gtex_fam)

	vander_patient_cov_df <- fread("../input/Komen.cov.tsv") %>%
		filter(Race %in% c("WHITE", "AFRNAMER")) %>%
		rename(id = Barcode)

	gtex_black <- fread("../input/gtex_AFR.list", col.names="id") %>% mutate(Race = "AFR")
	gtex_white <- fread("../input/gtex_EUR.list", col.names="id") %>% mutate(Race = "EUR")
	vander_gtex_patient_cov_df <- vander_patient_cov_df %>%
		bind_rows(gtex_black) %>%
		bind_rows(gtex_white) %>% as_tibble()

	rna_patient_df <- fread("../../00_nf_rnaseq/input/full_rnaseq_vander_sample_sheet.csv")
	w.gencode_patient_df <- tibble(id = names(fread(glue("../output/WHITE.Komen.GTEx.Breast.RBINT.TMM.TMP.Expr.GENCODEv26.hg38.txt")))) %>%
		mutate(id = case_when(str_detect(id, "GTEX") ~ str_replace(id, "\\.", "-"),
													TRUE ~ id)) %>%
		filter(id %in% vander_gtex_patient_cov_df$id)
	b.gencode_patient_df <- tibble(id = names(fread(glue("../output/BLACK.Komen.GTEx.Breast.RBINT.TMM.TMP.Expr.GENCODEv26.hg38.txt")))) %>%
		mutate(id = case_when(str_detect(id, "GTEX") ~ str_replace(id, "\\.", "-"),
													TRUE ~ id)) %>%
		filter(id %in% vander_gtex_patient_cov_df$id)

	gencode_patient_df <- bind_rows(w.gencode_patient_df, b.gencode_patient_df)

	patients_with_model_dat_df <- inner_join(eg_fam_df, vander_gtex_patient_cov_df, by = c("id")) %>%
		filter(#id %in% rna_patient_df$sample,
					 id %in% gencode_patient_df$id) %>%
	  mutate(id = case_when(str_detect(fid, "Komen") ~ fid,
													TRUE ~ id))

	# For now, just doing two groups, white and non-white
	white_pat_with_geno_df <- patients_with_model_dat_df %>%
		filter(Race %in% c("WHITE", "EUR")) %>%
		mutate(clust = "WHITE") %>%
		select(fid, id, clust)

	black_pat_with_geno_df <- patients_with_model_dat_df %>%
		filter(Race %in% c("AFR", "AFRNAMER")) %>%
		mutate(clust = "BLACK") %>%
		select(fid, id, clust)

	set.seed(209)
	white_train <- white_pat_with_geno_df %>%
		sample_frac(.8)
	set.seed(415)
	black_train <- black_pat_with_geno_df %>%
		sample_frac(.8)

	white_test <- anti_join(white_pat_with_geno_df, white_train, by = "id")
	black_test <- anti_join(black_pat_with_geno_df, black_train, by = "id")

	write_tsv(white_pat_with_geno_df, glue("../output/splits/WHITE.all.tsv"), col_names=F)
	write_tsv(black_pat_with_geno_df, glue("../output/splits/BLACK.all.tsv"), col_names=F)
	write_tsv(white_train, glue("../output/splits/WHITE.train.tsv"), col_names=F)
	write_tsv(white_test, glue("../output/splits/WHITE.test.tsv"), col_names=F)
	write_tsv(black_train, glue("../output/splits/BLACK.train.tsv"), col_names=F)
	write_tsv(black_test, glue("../output/splits/BLACK.test.tsv"), col_names=F)
	message("Samples written to ../output/splits/")
}

if (interactive()){
	main()
} else {
	main()
}
