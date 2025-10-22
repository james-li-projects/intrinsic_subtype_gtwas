main <- function(test=T){
	require(tidyverse)
	require(data.table)
	require(purrr)
	require(glue)

	cov_df <- fread("../input/Komen.cov.tsv")

	make_one_fastq_row <- function(cur_bc){
		fastq_files <- list.files("../input/Komen_RNASeq", pattern=cur_bc, full.names=T)
		if (length(fastq_files) == 0){
			print(glue("No fastq files found for {cur_bc}"))
			return(tibble()) # No fastq files for these patients
		}
		return_tib <- tibble(Barcode=cur_bc, fastq_1=fastq_files[1],
												 fastq_2=fastq_files[2], strandedness="auto")
	}

	barcodes <- cov_df %>% pull(Barcode)
	fastq_df <- barcodes %>% map_dfr(~make_one_fastq_row(.))

	samp_sheet_df <- inner_join(cov_df, fastq_df) %>%
		select(sample=Barcode, fastq_1, fastq_2, strandedness)

	if (isTRUE(test)){
		prefix <- "test"
		samp_sheet_df <- samp_sheet_df %>%
			head(10)
	} else {
		prefix <- "full"
	}
	write_csv(samp_sheet_df, glue("../input/{prefix}_rnaseq_vander_sample_sheet.csv"))

}

if (interactive()){
	main()
	main(F)
} else {
	main()
	main(F)
}
