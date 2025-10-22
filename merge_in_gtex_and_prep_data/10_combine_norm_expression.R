main <- function(group="WHITE"){
	require(edgeR)
	require(data.table)
	require(tidyverse)
	require(RNOmni)
	require(PCAForQTL)
	require(glue)

	if (!file.exists(glue("../output/{group}.Komen.GTEx.Breast.RBINT.TMM.TMP.Expr.GENCODEv26.hg38.txt"))){
		message(glue("{group} | Merged vander and gtex expression does not exist, creating . . "))
		merge_and_rb_int(group)
	} else {
		message(glue("{group} | Merged vander and gtex expression exists, reading. . . "))
	}

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

	dataGeneExpressionFP <- fread(glue("../output/{group}.Komen.GTEx.Breast.RBINT.TMM.TMP.Expr.GENCODEv26.hg38.txt"))

	# Get intersection of genotype + gene expression IDs for current group
	expr <- dataGeneExpressionFP %>%
		select(-Name, -genename) %>%
		t()

	genotyped_and_expr <- intersect(genoPCs$id, rownames(expr))
	cov_df <- get_covs(genotyped_and_expr)

	# expr <- expr[genotyped_and_expr,]
	genoPCs <- genoPCs %>%
		filter(id %in% genotyped_and_expr) %>%
		arrange(match(id, genotyped_and_expr))
	rownames(genoPCs) <- genoPCs %>% pull(id)

	prcompResult <- prcomp(expr, center=TRUE, scale.=TRUE)
	PCs <- prcompResult$x

	resultRunElbow <- PCAForQTL::runElbow(prcompResult=prcompResult)
  print(resultRunElbow)

	PCsTop <- PCs[,1:resultRunElbow]
	geno_pc_names <- rownames(PCsTop)

	write_tsv(PCsTop %>%
					 	  as_tibble() %>%
							mutate(id=geno_pc_names, .before="PC1") %>%
							filter(id %in% genotyped_and_expr),
						glue("../output/pcs/{group}.expression_vander_gtex_top_pcs.tsv"))

	# Can't get filterKnownCovariates to work with filtering only genoPCs
	#
	# # Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
	# #   NA/NaN/Inf in 'y'
	# # In addition: Warning message:
	# # In storage.mode(v) <- "double" : NAs introduced by coercion
	pc_combined <- bind_cols(genoPCs %>% select(-id, -gtex_or_vander), PCsTop[genotyped_and_expr,])
	rownames(pc_combined) <- genotyped_and_expr
	# knownCovariatesFiltered <- PCAForQTL::filterKnownCovariates(cov_df %>% select(-id), pc_combined, unadjustedR2_cutoff=0.9)
	# print(names(knownCovariatesFiltered))
	cutoff <- 0.9

	gtex_member_and_pcs(cov_df=cov_df, PCs=PCsTop[genotyped_and_expr,], group=group)
	while (cutoff > 0){
		cov_df <- cov_df %>% mutate(gtex_or_not = case_when(str_detect(id, "GTEX") ~ TRUE, TRUE ~ FALSE)) %>%
			mutate(gtex_or_not = as.factor(gtex_or_not))
		knownCovariatesFiltered <- PCAForQTL::filterKnownCovariates(cov_df %>%
	 	select(-id), pc_combined, unadjustedR2_cutoff=cutoff)
		print(cutoff)
		print(names(knownCovariatesFiltered))
		cutoff <- cutoff - 0.1
	}
	# knownCovariatesFiltered <- PCAForQTL::filterKnownCovariates(cov_df %>% select(-id), PCsTop, unadjustedR2_cutoff=0.9)
	# knownCovariatesFiltered <- PCAForQTL::filterKnownCovariates(genoPCs %>% select(-id), PCsTop, unadjustedR2_cutoff=0.9)

	RNGkind("L'Ecuyer-CMRG")
	set.seed(1)
	resultRunBE<-PCAForQTL::runBE(expr,B=20,alpha=0.05)
	K_elbow<-resultRunElbow #12.
	K_BE<-resultRunBE$numOfPCsChosen #29.
	K_GTEx<-60 #GTEx uses 60 PEER factors, and they are almost identical to the top 60 PCs.
	PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE","GTEx"),values=c(K_elbow,K_BE,K_GTEx),
													 titleText=glue("{group}"))
	ggplot2::ggsave(glue("../output/pc_fig/{group}_gene.jpg"),width=16,height=11,unit="cm")
	return()
}

get_covs <- function(ids){
	vander_cov <- fread("../input/Komen.cov.tsv") %>%
		filter(Barcode %in% ids) %>%
		select(id=Barcode, Age, BMI)
	gtex_cov <- fread("../input/gtex_cov.tsv") %>%
		mutate(SUBJID = str_replace(SUBJID, "-", "\\.")) %>%
		filter(SUBJID %in% ids) %>%
		select(id=SUBJID, Age=AGE, BMI)


	return_cov <- bind_rows(vander_cov, gtex_cov) %>%
		arrange(match(id, ids))
	rownames(return_cov) <- ids
	return(return_cov)
}

merge_and_rb_int <- function(group="WHITE"){
	require(glue)
	# tutorial link for reference: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf
	#######################################################
	# Importing and joining vandy and GTEx count matrices #
	#######################################################

	gtex_afr <- fread("../input/gtex_AFR.list", col.names="id") %>% mutate(Race = "AFR") %>%
		mutate(id = str_replace(id, "-", "\\."))
	gtex_eur <- fread("../input/gtex_EUR.list", col.names="id") %>% mutate(Race = "EUR") %>%
		mutate(id = str_replace(id, "-", "\\."))
	vander_patient_cov_df <- fread("../input/Komen.cov.tsv") %>%
		filter(Race %in% c("WHITE", "AFRNAMER")) %>%
		rename(id = Barcode)

	# importing vanderbilt count matrix 
	vandy_counts <- fread("../input/Komen.RNAseQC.counts.GENCODEv26.hg38.txt") 
	vandy_counts <- data.frame(vandy_counts)

	# importing GTEx count matrix
	gtex_counts <- fread("/gpfs/data/huo-lab/BCAC/james.li/GTEx/GTEx_counts/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct.gz")
	gtex_counts <- data.frame(gtex_counts)
	gtex_counts <- gtex_counts %>% select(-id)
	names(gtex_counts)[3:dim(gtex_counts)[2]] <- names(gtex_counts)[3:dim(gtex_counts)[2]] %>% str_extract("(GTEX\\.[\\w\\d]{1,}).", group = 1) # Extract IDs that are similar to ones we expect later on: e.g. "GTEX.XY12"

	if (group == "WHITE"){
		gtex_counts <- gtex_counts[,c("Name", "Description", intersect(gtex_eur$id, colnames(gtex_counts)))]
		vandy_counts <- vandy_counts[,c("Name", "genename", intersect(vander_patient_cov_df %>% filter(Race == "WHITE") %>% pull(id), colnames(vandy_counts)))]
	} else if (group == "BLACK"){
		gtex_counts <- gtex_counts[,c("Name", "Description", intersect(gtex_afr$id, colnames(gtex_counts)))]
		vandy_counts <- vandy_counts[,c("Name", "genename", intersect(vander_patient_cov_df %>% filter(Race == "AFRNAMER") %>% pull(id), colnames(vandy_counts)))]
	}

	# joining the two count matrices and formatting them: making rownames the ENSG IDs and the column names the sample IDs
	Counts <- inner_join(vandy_counts,gtex_counts,by=c("Name" = "Name", "genename" = "Description"))

# 	Counts_other <- inner_join(vandy_counts %>%
# 														   mutate(no_dec = str_extract(Name, "ENSG\\d{1,}")),
# 													   gtex_counts %>%
# 														   mutate(no_dec = str_extract(Name, "ENSG\\d{1,}")),
# 														 by= "Name")
# 	Counts_strict <- inner_join(vandy_counts %>%
# 														   mutate(no_dec = str_extract(Name, "ENSG\\d{1,}")),
# 													   gtex_counts %>%
# 														   mutate(no_dec = str_extract(Name, "ENSG\\d{1,}")),
# 														 by= c("Name" = "Name", "genename" = "Description"))

	ensg_to_name_df <- Counts %>%
		select(Name, genename) %>%
		distinct()

	rownames(Counts) <- Counts$Name
	Counts <- Counts %>% select(-Name,-genename)
	count_mat <- data.matrix(Counts)

	############################################
	# Importing Joined Count Matrix into EdgeR #
	############################################

	# importing the count matrix into an edgeR object called DGEList
	dgList <- DGEList(counts=count_mat, genes=rownames(count_mat))

	#########################
	# Filtering the dataset #
	#########################

	total_num_samples <- nrow(dgList$samples)

	# filtering based on TPM>=0.1 in >=20% samples
	countsPerMillion <- cpm(dgList)
	countCheck <- countsPerMillion >= 0.1
	keep <- which(rowSums(countCheck) >= 0.2*total_num_samples)
	dgList <- dgList[keep,]

	# filtering based on Raw read count>=6 in >=20% samples
	rawReadCounts <- dgList$counts
	readCountCheck <- rawReadCounts >= 6
	keep <- which(rowSums(readCountCheck) >= 0.2*total_num_samples)
	dgList <- dgList[keep,]

	#######################################
	# Performing TMM Normalization of TPM #
	#######################################

	dgList <- calcNormFactors(dgList, method="TMM")
	TMM_TPM_Expr <- cpm(dgList)

	# Rank based inverse norm by row (gene expression)
	RB_INT_TMM_TPM_Expr <- apply(TMM_TPM_Expr, 1, RankNorm) %>% t() %>% as_tibble()
	first_name <- names(RB_INT_TMM_TPM_Expr)[1]

	RB_INT_TMM_TPM_Expr <- RB_INT_TMM_TPM_Expr %>%
		mutate(Name = rownames(TMM_TPM_Expr), .before = !!first_name)

	write_df <- ensg_to_name_df %>%
		inner_join(RB_INT_TMM_TPM_Expr, by = "Name")

	# outputting TMM normalized TPM 
	# save(TMM_TPM_Expr,file="/gpfs/data/huo-lab/Vanderbilt/james.li/TMM_TPM_Expr.RData")
	oname <- glue("../output/{group}.Komen.GTEx.Breast.RBINT.TMM.TMP.Expr.GENCODEv26.hg38.txt") 

	write_tsv(write_df, oname)
	message(glue("{oname} written"))
}

gtex_member_and_pcs <- function(cov_df, PCs, group){
	require(glmnet)
	require(glue)

	cov_df <- cov_df %>% mutate(gtex_or_not = case_when(str_detect(id, "GTEX") ~ TRUE, TRUE ~ FALSE))

	glmmod <- glmnet(PCs, y=as.factor(cov_df$gtex_or_not), alpha=1, family="binomial")

	# Do correlations
	high_cor <- 0
	high_cor_pc <- ""
	for (cname in colnames(PCs)){
		cur_cor <- cor(PCs[,cname], cov_df$gtex_or_not)
		print(glue("{group} | {cname} has corr: {cur_cor} with gtex_or_not"))
		if (abs(cur_cor) > high_cor){
			high_cor <- cur_cor
			high_cor_pc <- cname
		}
	}
	print(glue("Highest correlation ({high_cor}) PC for {group}: {high_cor_pc}"))
}

if (interactive()){
	main()
	main("BLACK")
} else {
	main()
	main("BLACK")
}
