main <- function(group="WHITE"){
	require(tidyverse)
	require(data.table)
	require(glue)
	shared_df <- fread("../output/vander_rsq_filt_and_gtex_overlap_vars.txt", col.names=c("SNP"))

	master_keep <- tibble()
	master_exclude <- tibble()
	for (chrom in 1:22){
		gtex_freq <- fread(glue("../output/freq/gtex.{group}.chr{chrom}.frq.strat")) %>%
			filter(SNP %in% shared_df$SNP)

		vander_freq <- fread(glue("../output/freq/vander.{group}.chr{chrom}_rsq0.7.frq.strat")) %>%
			filter(SNP %in% shared_df$SNP)

		gtex_vander_freq <- gtex_freq %>%
			select(SNP, g.MAF=MAF) %>%
			inner_join(vander_freq %>% select(SNP, v.MAF=MAF), by = "SNP") %>%
			mutate(dif.MAF = abs(g.MAF - v.MAF))

		og_shared <- nrow(gtex_vander_freq)

		gtex_vander_keep <- gtex_vander_freq %>%
			filter(dif.MAF < 0.15)
	  gtex_vander_exclude <- anti_join(gtex_vander_freq, gtex_vander_keep, by = "SNP")

		maf_filt_shared <- nrow(gtex_vander_keep)

		master_keep <- bind_rows(master_keep, gtex_vander_keep %>% select(SNP))
		master_exclude <- bind_rows(master_exclude, gtex_vander_exclude %>% select(SNP))

		message(glue("{group} | Chr{chrom}: {og_shared} to {maf_filt_shared} variants ({og_shared-maf_filt_shared} removed)"))

		write_tsv(gtex_vander_keep %>% select(SNP), glue("../output/race_maf_dif_filt_and_vander_rsq_filt_lists/keep.{group}.chr{chrom}.txt"), col_names=F)
		# write_tsv(gtex_vander_exclude %>% select(SNP), glue("../output/race_maf_dif_filt_and_vander_rsq_filt_lists/exclude.{group}.chr{chrom}.txt"), col_names=F)
	}
	write_tsv(master_keep, glue("../output/race_maf_dif_filt_and_vander_rsq_filt_lists/keep.{group}.chrALL.txt"), col_names=F)
	# write_tsv(master_exclude, glue("../output/race_maf_dif_filt_and_vander_rsq_filt_lists/exclude.{group}.chrALL.txt"), col_names=F)
}


if (interactive()){
	main()
	main("BLACK")
} else {
	main("BLACK")
}
