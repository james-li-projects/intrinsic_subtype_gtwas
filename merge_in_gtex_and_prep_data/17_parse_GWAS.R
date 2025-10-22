SUBTYPES <- c("HRNEG_HER2NEG", "HRNEG_HER2POS", "HRPOS_HER2NEG",
							"HRPOS_HER2POS")
ANCESTRIES <- c("eur", "afr")




# weights
# gene	rsid	varID	ref_allele	eff_allele	weight
# ENSG00000000419.12	NA	chr20:49960606:A:G	A	G	0.00212503299873433
# ENSG00000000419.12	NA	chr20:49962465:G:A	G	A	0.0047444638540068895
# ENSG00000000419.12	NA	chr20:49962811:G:T	G	T	0.001626264821874481
# ENSG00000000419.12	NA	chr20:49963473:G:C	G	C	5.968002020021716e-4
# ENSG00000000419.12	NA	chr20:49965204:T:C	T	C	1.6656478154722823e-4
# ENSG00000000419.12	NA	chr20:49979436:G:T	G	T	-0.02143471029117039
# ENSG00000000419.12	NA	chr20:50029855:C:T	C	T	-0.0296760168701782
# ENSG00000000419.12	NA	chr20:50353984:G:A	G	A	-0.13801125061690067
# ENSG00000000419.12	NA	chr20:50402316:C:T	C	T	-0.0045109016900066285

# variant_id	panel_variant_id	chromosome	position	effect_allele	non_effect_allele	frequency	pvalue	zscore	effect_size	standard_error	sample_size	n_cases
# rs888953847	NA	chr1	594445	T	C	0.0006	0.9626	-0.046891127	-0.0152	0.3236	472238	75024
# rs1040232850	NA	chr1	595762	CTG	C	0.9986	0.267	-1.1099977	-0.2494	0.2247	472238	75024
# rs1390538076	NA	chr1	630947	A	G	0.0003	0.8214	0.22574468	0.0975	0.43200000000000005	472238	75024
# rs1250812823	NA	chr1	646157	A	G	0.0013	0.3411	-0.95199406	-0.3028	0.318	472238	75024
# rs991450070	NA	chr1	726526	C	G	0.0011	0.3542	0.92647344	0.2938	0.3171	472238	75024
# rs151190501	chr1_727233_G_A_b38	chr1	727233	A	G	0.0191	0.6424	-0.4643458	-0.021	0.0453	472238	75024
# rs61769339	chr1_727242_G_A_b38	chr1	727242	A	G	0.1105	0.04962	1.9632254	0.0328	0.0167	472238	75024
# rs1362083662	NA	chr1	730012	T	G	0.0009	0.9441	-0.070117675	-0.0319	0.4544	472238	75024
# rs1039158266	NA	chr1	730804	T	G	0.0013	0.8176	0.2306329	0.054000000000000006	0.2343	472238	75024

main <- function(){
	library(tidyverse)
	library(data.table)
	library(glue)

	for (stype in SUBTYPES){
		for (ances in ANCESTRIES){
			input_df_path <- glue("../input/{ances}_meta_sumstats/{stype}.tsv")
			
			cur_df <- fread(input_df_path) %>%
				mutate(frequency=.1, variant_id=NA, sample_size=472238, n_cases=75024) %>%
				rename(effect_allele=EffectAllele, non_effect_allele=BaselineAllele, effect_size=BETA, standard_error=SE, pvalue=P) %>%
				separate(ID, c("chromosome", "position", NA, NA), sep=":") %>%
				mutate(chromosome_fixed = glue("chr{chromosome}"),
							 panel_variant_id=glue("{chromosome_fixed}:{position}:{non_effect_allele}:{effect_allele}"),
							 zscore = effect_size / standard_error
							 ) %>%
				mutate(chromosome=chromosome_fixed) %>%
				select(variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, frequency, pvalue, zscore, effect_size, standard_error, sample_size, n_cases)

			# For later on, since I named stuff WHITE and BLACK after
			if (ances == "eur"){
				ances_tag <- "WHITE"
			} else if (ances == "afr"){
				ances_tag <- "BLACK"
			}

			write_tsv(cur_df, glue("../output/james_meta_gwas_parsed/{ances_tag}.{stype}.txt.gz"))
			message(glue("Finished {input_df_path}"))
		}
	}
}

if (interactive()){
	main()
} else {
	main()
}
