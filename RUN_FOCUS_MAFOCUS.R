# FOCUS
subtypes=("HRPOS_HER2NEG" "HRPOS_HER2POS" "HRNEG_HER2POS" "HRNEG_HER2NEG")
races=("BLACK" "WHITE")
ancestries=("afr" "eur")
celltypes=("Adipocytes" "Breast_tissue" "Endothelial_cells" "Epithelial_cells" "Stromal_and_Immune_cells")
for subtype in "${subtypes[@]}"; do
for celltype in "${celltypes[@]}"; do
for i in "${!races[@]}"; do
race=${races[$i]}
ancestry=${ancestries[$i]}
upper_ancestry=$(echo "$ancestry" | tr '[:lower:]' '[:upper:]')
/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/software/ma-focus/bin/focus finemap \
/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/FOCUS/input_gwas/${subtype}.${race}.tsv.gz \
/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/${ancestry} \
/gpfs/data/huo-lab/Vanderbilt/Julian/05_ma_focus/input/focus_imported_dbs/${race}.gene_models.${celltype}.db \
--prior-prob "gencode38" --locations 38:${upper_ancestry} \
--out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/FOCUS/output/focus_${subtype}_${race}_${celltype}.tsv
done
done
done

# MA-FOCUS
subtype_list=("HRPOS_HER2NEG" "HRPOS_HER2POS" "HRNEG_HER2POS" "HRNEG_HER2NEG")
celltype_list=("Adipocytes" "Breast_tissue" "Endothelial_cells" "Epithelial_cells" "Stromal_and_Immune_cells")
for subtype in "${subtype_list[@]}"; do
for celltype in "${celltype_list[@]}"; do
/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/software/ma-focus/bin/focus finemap \
/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/FOCUS/input_gwas/${subtype}.WHITE.tsv.gz:/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/FOCUS/input_gwas/${subtype}.BLACK.tsv.gz \
/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/eur:/gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/data/LD_REFERENCE_PANELS/afr \
/gpfs/data/huo-lab/Vanderbilt/Julian/05_ma_focus/input/focus_imported_dbs/WHITE.gene_models.${celltype}.db:/gpfs/data/huo-lab/Vanderbilt/Julian/05_ma_focus/input/focus_imported_dbs/BLACK.gene_models.${celltype}.db \
--prior-prob "gencode38" --locations 38:EUR-AFR \
--out /gpfs/data/huo-lab/BCAC/james.li/XANCESTRY_GWAS_TWAS/output/TWAS/FOCUS/output/mafocus_${subtype}_${celltype}.tsv
done
done


