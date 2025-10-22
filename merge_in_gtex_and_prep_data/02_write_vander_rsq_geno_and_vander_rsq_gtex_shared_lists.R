main <- function(){
	write_geno_lists()
}

write_geno_lists <- function(maf=.01, rsq=.7){
	require(tidyverse)
	require(data.table)
	require(glue)

	gtex_vars <- fread("../input/gtex_v8_whole_chr_bbf/processed_gtex.bim")
	master_df <- tibble()
	for (cur_chr in 1:22){
		vars_df <- fread(glue("../input/ImputedGenotyping_TopMed/chr{cur_chr}.info.gz"))
		filtered_df <- vars_df %>%
			filter(Rsq > rsq) %>%
		  select(SNP, MAF, `REF(0)`,	`ALT(1)`)
		message(glue("Chr{cur_chr} {nrow(vars_df)} to {nrow(filtered_df)} rows"))

		master_df <- bind_rows(master_df, filtered_df)
		write_delim(filtered_df, glue("../output/vander_rsq_filtered_var_lists/chr{cur_chr}_variants_rsq{rsq}.txt"),
								col_names=F)
	}
	write_delim(master_df %>% select(SNP), glue("../output/vander_rsq_filtered_var_lists/chrALL_variants_rsq{rsq}.txt"),
							col_names=F)

	# in bim file, a1 is alt (minor), a2 is ref (major)
	# https://www.cog-genomics.org/plink2/formats#bim

	# Write vars shared between rsq filtered vander and gtex
	shared_vars <- tibble(id = intersect(gtex_vars$V2, master_df$SNP))
	write_delim(shared_vars, glue("../output/vander_rsq_filt_and_gtex_overlap_vars.txt"),
							col_names=F)
}

if (interactive()){
	main()
} else {
	main()
}
