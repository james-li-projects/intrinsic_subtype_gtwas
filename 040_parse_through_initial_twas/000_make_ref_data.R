# Using bbf files used in gene expression model creation, create files for each
# ancestry to list variants 

# Also creates necessary files for metaxcan/alvaro parsing for COJO LD references
arg_parser <- function(){
  require(argparse)
  parser <- ArgumentParser(description=paste0("",
                                              "",
                                              "",
                                              "",
                                              ""))
  parser$add_argument("--ancestry", type="character", required=T, default=NULL,
                      help="black or white"
  )
	parser$add_argument("--cojo_ld", action = "store_true",
										help = "Enable COJO LD stuff", default = FALSE)

  # Can usually use function defaults here

  # Put in {no_ext} to put the name of the db with no extension in the out ld,
  # put {part_num} to put the part number in (for big lds that need to be split.)
  parser$add_argument("--scratch_path", type="character", required=F, default="/scratch/jmcclellan/tmp/vander/",
                      help="Scratch directory to store LD bigsnpr backingfiles."
  )
  parser_args  <- parser$parse_args()
  return(parser_args)
}
main <- function(ancestry,
								 scratch_path = "/scratch/jmcclellan/tmp/vander/",
								 cojo_ld=F
								 ){
	library(tidyverse)
	library(data.table)
	library(glue)
	library(arrow)
	library(bigsnpr)
	library(tools)

	if (isTRUE(cojo_ld)){
		vars_fname <- glue("../input/cojo_ld_ref_jcm/COJO.{ancestry}.variants.txt.gz") 
		if (file.exists(vars_fname)){
			return()
		}
	} else {
		vars_fname <- glue("../input/model_geno_ref/{ancestry}.variants.txt.gz") 
	}

	alt_ancestry <- switch(ancestry, "BLACK"="afr", "WHITE"="eur")

	if (isTRUE(cojo_ld)){
		bed_path <- glue("../input/LD_REFERENCE_PANELS/{alt_ancestry}.bed")
	} else {
		bed_path <- glue("../input/vander_gtex_genotypes/{ancestry}.chrALL_vander_gtex.bed")
	}

	fbase <- glue("{scratch_path}{thing}",
								thing=file_path_sans_ext(basename(bed_path)))
	frds <-  glue("{fbase}.rds")

	if (file.exists(frds)){
		rds <- frds
	} else {
		rds <- snp_readBed(bed_path, backingfile=fbase)
	}
	ref_obj <- snp_attach(rds)

	ref_obj$map <- ref_obj$map %>%
		mutate(idx = row_number())

	if (isTRUE(cojo_ld)){
		variant_df <- ref_obj$map %>%
			as_tibble() %>%
			mutate(allele_1_frequency=0.15) %>% # placeholder freq, unused going forward
			mutate(rsid=marker.ID) %>%
			select(chromosome, position=physical.pos, id=marker.ID, allele_0=allele2, allele_1=allele1, allele_1_frequency, rsid)
	} else {
		variant_df <- ref_obj$map %>%
			as_tibble() %>%
			mutate(new_id=str_replace(`marker.ID`, "chr", "")) %>%
			mutate(allele_1_frequency=0.15) %>% # placeholder freq, unused going forward
			mutate(rsid=marker.ID) %>%
			select(chromosome, position=physical.pos, id=new_id, allele_0=allele2, allele_1=allele1, allele_1_frequency, rsid)
	}

	metadata_fname <- glue("../input/model_geno_ref/{ancestry}.metadata.parquet")

	if (!file.exists(vars_fname)){
		write_tsv(variant_df, vars_fname)
	}

	if (isTRUE(cojo_ld)){
		return()
	}

  if (!file.exists(metadata_fname)){
		# For alvaro's parquet reading, need to have metadata in chromosome sorted row groups
	  write_dataset(variant_df %>% mutate(chrom=chromosome), metadata_fname, partitioning="chrom", format="parquet")
	}

	indv_names <- ref_obj$fam %>% pull(sample.ID)

	# Median impute
	missing_counts <- big_counts(ref_obj$genotype)
	missing_indx <- which(missing_counts[4,] > 0, arr.ind = TRUE)

	if (length(missing_indx) > 0){
		for (cmi in missing_indx){
			cur_col_med <- median(ref_obj$genotype[,cmi], na.rm = TRUE)
			na_rows <- which(is.na(ref_obj$genotype[,cmi]))
			for (crow in na_rows){
				ref_obj$genotype[crow,cmi] <- cur_col_med
			}
		}
	}

	for (chr_num in unique(ref_obj$map$chromosome)){
		chr_idx <- ref_obj$map %>%
		 filter(chromosome == chr_num) %>%
		 pull(idx)

		var_names <- ref_obj$map %>%
		 filter(chromosome == chr_num) %>%
		 mutate(new_id=str_replace(`marker.ID`, "chr", "")) %>%
		 pull(new_id)

		write_df <- ref_obj$genotype[,chr_idx] %>%
			as_tibble() %>%
			mutate(individual=indv_names, .before=V1)
		names(write_df) <- c("individual", var_names)

		chr_parq_name <- glue("../input/model_geno_ref/{ancestry}.chr{chr_num}.geno.parquet")
		if (!file.exists(chr_parq_name)){
			write_parquet(write_df, chr_parq_name)
		}
	}
	message(glue("Finished {ancestry}"))
}

if (interactive()){
	# main("BLACK")
	# main("WHITE")
	# main("BLACK", cojo_ld=T)
	# main("WHITE", cojo_ld=T)
} else {
  options(error = function() traceback(2))
  sink(stdout(), type="message")
  cli_args <- arg_parser()
  print(cli_args)
  do.call(main, cli_args)
}
