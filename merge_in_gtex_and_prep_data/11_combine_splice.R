main <- function(group="WHITE"){
	require(edgeR)
	require(data.table)
	require(tidyverse)
	require(RNOmni)
	require(PCAForQTL)
	require(glue)
	print(glue("{group}"))

	# cluster IDs from vander data have the form: clu_<number>_{+ or -}, so we
	# remove the extraneous bit to make it just clu_<number> to try and join with
	# gtex cluster ID
	vander_splice <- read_vander_splice() %>%
		separate(ID, into = c("chr", "start", "end", "cluster"), sep=":") %>%
		mutate(cluster = str_replace(cluster, "_-", "")) %>%
		mutate(cluster = str_replace(cluster, "_\\+", "")) %>%
		mutate(chr = as.integer(chr)) %>%
		select(-`#Chr`)

	# ID:  "chr1:14829:14970:clu_52500:ENSG00000227232.5"
	gtex_counts <- fread("/gpfs/data/huo-lab/BCAC/james.li/GTEx/GTEx_Analysis_v8_sQTL_phenotype_matrices/Breast_Mammary_Tissue.v8.leafcutter_phenotypes.bed.gz") %>%
		separate(ID, into = c("chr", "start", "end", "cluster", "gene"), sep=":") %>%
		mutate(chr = str_extract(chr, "\\d{1,}")) %>%
		mutate(chr = as.integer(chr)) %>%
		select(-`#Chr`) %>%
		distinct()

	joint_df <- vander_splice %>%
		inner_join(gtex_counts, by = c("chr", "start", "end", "cluster"))
	joint_df_no_cluster <- vander_splice %>%
		select(-cluster) %>%
		inner_join(gtex_counts %>% select(-cluster), by = c("chr", "start", "end")) %>%
		relocate(gene, .after = end)
	
	group_file <- glue("../output/splits/{group}.all.tsv")
	if (!file.exists(group_file)){
		stop(glue("{group} membership file doesn't exist yet, go ahead and finish subsequent steps then return here."))
	}
	group_df <- fread(group_file, header=F) %>%
		mutate(V2 = case_when(#str_detect(V2, "GTEX") ~ str_replace(V2, "-", "\\."),
													str_detect(V2, "Komen") ~ str_replace(V2, "Komen-BC-", ""),
													TRUE ~ V2))

	pc_file <- glue("../output/pcs/{group}.vander_gtex_chrALL_pruned_genotypes_PCs.eigenvec")

	if (!file.exists(pc_file)){
		stop("PC file doesn't yet exist for the group, go ahead and get that finished up first in subsequent steps then return here.")
	}

	genoPCs <- fread(pc_file) %>%
		select(-`#FID`)
	names(genoPCs) <- c("id", paste0("genoPC", 1:10))

	genoPCs <- genoPCs %>%
		mutate(gtex_or_vander = case_when(str_detect(id, "GTEX") ~ "GTEX",
																			TRUE ~ "vander")) %>%
	  mutate(id = case_when(str_detect(id, "GTEX") ~ str_replace(id, "-", "\\."),
													str_detect(id, "Komen") ~ str_replace(id, "Komen-BC-", ""),
													TRUE ~ id))
	# Calculate splicing group PCs

	  # mutate(id = case_when(str_detect(id, "GTEX") ~ str_replace(id, "-", "\\."),
		# 											str_detect(id, "Komen") ~ str_replace(id, "Komen-BC-", ""),
		# 											TRUE ~ id)) %>%
	gtex_afr <- fread("../input/gtex_AFR.list", col.names="id") %>% mutate(Race = "AFR")
	gtex_eur <- fread("../input/gtex_EUR.list", col.names="id") %>% mutate(Race = "EUR")
	vander_patient_cov_df <- fread("../input/Komen.cov.tsv") %>%
		filter(Race %in% c("WHITE", "AFRNAMER")) %>%
		rename(id = Barcode)

	if (group == "WHITE"){
		group_ids <- c(vander_patient_cov_df %>%
									 filter(Race == "WHITE") %>%
									 pull(id), gtex_eur$id)
	} else if (group == "BLACK"){
		group_ids <- c(vander_patient_cov_df %>%
									 filter(Race == "AFRNAMER") %>%
									 pull(id), gtex_afr$id)
	}

	joint_df_write <- joint_df_no_cluster %>%
		mutate(Name = glue("{chr}:{start}:{end}"), .before=chr) %>%
		select(-all_of(c("chr", "start", "end")))
	names(joint_df_write) <- names(joint_df_write) %>% str_replace("\\.markdup\\.sorted\\.bam", "")
	# joint_df_write <- joint_df_write[,c("Name", "gene", intersect(group_ids, names(joint_df_write)))] %>%
	# 	distinct(Name, .keep_all=T)
	joint_df_write <- joint_df_write[,c("Name", "gene", intersect(group_ids, names(joint_df_write)))]

	# Probably shouldn't have to pivot longer then wider, but it works
	group_splice <- joint_df_no_cluster %>%
		distinct(chr, start, end, .keep_all=T) %>%
		pivot_longer(5:dim(joint_df_no_cluster)[2], names_to="id", values_to="qqnorm") %>%
		pivot_wider(names_from = c("chr", "start", "end", "gene"), names_sep=":", values_from="qqnorm") %>%
		mutate(id = case_when(str_detect(id, "markdup") ~ str_replace(id, "\\.markdup\\.sorted\\.bam", ""),
													TRUE ~ id
													)) %>%
		filter(id %in% group_ids) %>%
		arrange(match(id, group_ids))

	group_ids_with_splice <- group_splice$id

	prcompResult <- prcomp(group_splice %>% select(-id), center=TRUE, scale.=TRUE)
	PCs <- prcompResult$x
	rownames(PCs) <- group_splice$id

	resultRunElbow <- PCAForQTL::runElbow(prcompResult=prcompResult)
  print(resultRunElbow)
	PCsTop <- PCs[,1:resultRunElbow]

	# Write splice
	write_tsv(joint_df_write, glue("../output/{group}.Komen.GTEx.Breast.leafcutter.Introns.txt"))
	write_tsv(joint_df_write %>% select(id=Name, gene), glue("../output/model_target_lists/{group}.full_non-repeat_intron_list.txt"))
	# Only write PCs for patients that also have genotype data
	write_tsv(PCsTop %>%
					 	  as_tibble() %>%
							mutate(id=group_ids_with_splice, .before="PC1") %>%
							filter(id %in% group_df$V2) %>%
							arrange(match(id, group_df$V2)),
						glue("../output/pcs/{group}.splicing_vander_gtex_top_pcs.tsv"))

	cov_df <- get_covs(group_df$V2)
	pc_combined <- bind_cols(genoPCs %>% arrange(match(id, group_df$V2)) %>% select(-id, -gtex_or_vander), PCsTop[group_df$V2,])
	rownames(pc_combined) <- group_df$V2
	cutoff <- 0.9

	while (cutoff > 0){
		cov_df <- cov_df %>% mutate(gtex_or_not = case_when(str_detect(id, "GTEX") ~ TRUE, TRUE ~ FALSE)) %>%
			mutate(gtex_or_not = as.factor(gtex_or_not))
		knownCovariatesFiltered <- PCAForQTL::filterKnownCovariates(cov_df %>%
	 	select(-id), pc_combined, unadjustedR2_cutoff=cutoff)
		print(cutoff)
		print(names(knownCovariatesFiltered))
		cutoff <- cutoff - 0.1
	}
	expr <- group_splice %>% select(-id) %>% as.matrix()

	RNGkind("L'Ecuyer-CMRG")
	set.seed(1)
	resultRunBE<-PCAForQTL::runBE(expr,B=20,alpha=0.05)
	K_elbow<-resultRunElbow #12.
	K_BE<-resultRunBE$numOfPCsChosen #29.
	K_GTEx<-60 #GTEx uses 60 PEER factors, and they are almost identical to the top 60 PCs.
	PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE","GTEx"),values=c(K_elbow,K_BE,K_GTEx),
													 titleText=glue("{group}"))
	ggplot2::ggsave(glue("../output/pc_fig/{group}_intron.jpg"),width=16,height=11,unit="cm")
	return()
}

get_covs <- function(ids){
	vander_cov <- fread("../input/Komen.cov.tsv") %>%
		filter(Barcode %in% ids) %>%
		select(id=Barcode, Age, BMI)
	gtex_cov <- fread("../input/gtex_cov.tsv") %>%
		filter(SUBJID %in% ids) %>%
		select(id=SUBJID, Age=AGE, BMI)


	return_cov <- bind_rows(vander_cov, gtex_cov) %>%
		arrange(match(id, ids))
	rownames(return_cov) <- ids
	return(return_cov)
}

read_vander_splice <- function(){
	master_df <- tibble()
	for (chr_num in 1:22){
		cur_df <- fread(glue("../input/leafcutter_output/Vander417_perind.counts.gz.qqnorm_{chr_num}"))
		master_df <- bind_rows(master_df, cur_df)
	}

	return(master_df)
}

if (interactive()){
	main()
	main("BLACK")
} else {
	main()
	main("BLACK")
}
