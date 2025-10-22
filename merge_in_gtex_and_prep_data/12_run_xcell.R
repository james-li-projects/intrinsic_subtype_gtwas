main <- function(group="WHITE"){
	require(xCell)
	require(tidyverse)
	require(glue)
	require(data.table)
	print(group)

	res <- runonce::save_run({
		exprMatrix <- read.table(glue("../output/{group}.Komen.GTEx.Breast.RBINT.TMM.TMP.Expr.GENCODEv26.hg38.txt"), header=TRUE, row.names=1, as.is=TRUE)
		alt_names <- exprMatrix[,1]
		ensg_names <- rownames(exprMatrix)

		# Get rid of repeated names (< ~ 200 or so given ~50K gene rows)
		names_tib <- tibble(ensg = rownames(exprMatrix), symb = alt_names)
		repeat_symbols <- names_tib %>%
			group_by(symb) %>%
			summarize(n = n()) %>%
			filter(n > 1)

		bad_ensg <- names_tib %>%
			filter(symb %in% repeat_symbols$symb) %>%
			pull(ensg)

		exprMatrix <- exprMatrix[!(rownames(exprMatrix) %in% bad_ensg),]

		ensg_names <- rownames(exprMatrix)

		ensg_df <- tibble(ensg=ensg_names) %>%
			mutate(ensg_no_dec = str_extract(ensg, "ENSG\\d{1,}"))
		if (length(bad_ensg) > 0){
		  write_tsv(tibble(id=bad_ensg), "../output/misc/ensg_for_repeat_gene_names__Komen.RNAseQC.TPM.GENCODEv26.hg38.tsv")
		}

		write_tsv(tibble(id=rownames(exprMatrix)), glue("../output/model_target_lists/{group}.full_non-repeat_gene_list.txt"))
		rownames(exprMatrix) <- exprMatrix[,1]

		browser()
		xCellAnalysis(exprMatrix, genes=rownames(exprMatrix), 
									cell.types.use = c("Epithelial cells", "Endothelial cells",
																		 "Fibroblasts", "Adipocytes", "B-cells",
																		 "CD4+ T-cells", "CD8+ T-cells", "DC",
																		 "Eosinophils", "Macrophages", "Monocytes",
																		 "Mast cells", "Neutrophils", "NK cells"
																		  ))
		# xCellAnalysis(exprMatrix, genes=rownames(exprMatrix))
	}, file = glue("../output/xCell/{group}.Komen.GTEx.Breast.cell_scores.rds"))
	return()
}

if (interactive()){
	main()
	main("BLACK")
} else {
	main()
	main("BLACK")
}
