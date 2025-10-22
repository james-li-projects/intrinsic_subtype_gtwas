# Converts hg19 AFR/EUR ld blocks to hg38.
#
# In the case of NAs that occur in the middle of a chr, split the difference
# with the next row's block (assuming next row's block has no NA) (do manually to file after)
#
# In the case NAs occur at the end of a chr manually add a sufficiently long
# enough distance to end of the chromosome, a few MB or so

ANCESTIRES <- c("WHITE", "BLACK")

main <- function(ancestry, chr_num){
	library(tidyverse)
	library(data.table)
	library(bigsnpr)
	library(glue)


	prefix <- switch(ancestry, "BLACK"="AFR",
									 "WHITE"="EUR")

	hg19_blocks_df <- fread(glue("../input/ld_blocks/{prefix}/fourier_ls-chr{chr_num}.bed")) %>%
		mutate(chr_old = chr,
					 start=as.numeric(start),
					 stop=as.numeric(stop)
					 ) %>%
		mutate(chr = str_extract(chr_old, "\\d{1,}")) 

	hg19_starts <- hg19_blocks_df %>%
		select(chr, pos=start) %>%
		mutate(row_num = row_number())

	hg19_stops <- hg19_blocks_df %>%
		select(chr, pos=stop) %>%
		mutate(row_num=row_number())

	hg38_starts <- snp_modifyBuild(hg19_starts,
																 "/gpfs/data/gao-lab/Julian/software/liftOver",
																 #check_reverse=F,
																 from="hg19", to="hg38")
																 #local_chain="../input/ld_blocks/hg19ToHg38.over.chain.gz")
	hg38_stops <- snp_modifyBuild(hg19_stops,
																 "/gpfs/data/gao-lab/Julian/software/liftOver",
																 #check_reverse=F,
																 from="hg19", to="hg38")
																 # local_chain="../input/ld_blocks/hg19ToHg38.over.chain.gz")

	hg38_starts <- hg38_starts %>%
		rename(start=pos)

	hg38_stops <- hg38_stops %>%
		rename(stop=pos)

	combined_df <- inner_join(hg38_starts, hg38_stops, by = c("chr", "row_num"))

	write_tsv(combined_df %>% select(-row_num), glue("../input/ld_blocks/hg38_ld/{ancestry}_chr{chr_num}_ld.bed"))
}

split_block_diff <- function(row1_start, row2_end){
	dist <- row2_end - row1_start
	row1_end <- row1_start + dist/2

	return(row1_end)
}

if (interactive()){
	for (chr_num in 1:22){
		main("WHITE", chr_num)
		main("BLACK", chr_num)
	}
} else {
}
