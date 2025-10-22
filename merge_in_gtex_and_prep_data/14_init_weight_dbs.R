main <- function(group="WHITE",
								 template_db_path="../input/template_predixcan.db" # Just a wiped version of gtex v8 db files
								 ){
	require(tidyverse)
	require(glue)
	require(data.table)

	model_names <- str_replace(c("Epithelial cells", "Endothelial cells", "Fibroblasts",
									 "Adipocytes", "Breast_tissue", "Stromal_and_Immune_cells"), " ", "_")

	for (cname in model_names){
		out_db_path <- glue("../output/predixcan_dbs/{group}.gene_models.{cname}.db")
		out_db_path2 <- glue("../output/predixcan_dbs/{group}.intron_models.{cname}.db")
		file.copy(template_db_path, out_db_path, overwrite=T)
		file.copy(template_db_path, out_db_path2, overwrite=T)
	}

}

if (interactive()){
	main()
	main("BLACK")
} else {
	main()
	main("BLACK")
}
