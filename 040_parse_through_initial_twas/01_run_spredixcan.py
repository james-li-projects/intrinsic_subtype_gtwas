#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import time
import glob
import re
import pandas as pd
import logging

import parsl
from parsl.dataflow.futures import AppFuture
from config import config, id_for_memo_File
from parsl.data_provider.files import File
from parsl.app.app import bash_app 
from parsl.app.errors import MissingOutputs, BashExitFailure
from functools import partial
from math import ceil
from datetime import datetime

from options import MetaxcanParse
from metaxcan_parse import metaxcan_parse_gwas
from metaxcan_impute import metaxcan_impute_gwas
from metaxcan_spredixcan import spredixcan
from cojo_prep import filter_gwas_by_regions, combine_spredixcans


config.retries = 2
parsl.load(config)
parsl.set_stream_logger(level=logging.DEBUG)

# CELL_TYPES = ["Stromal_and_Immune_cells", "Adipocytes", "Breast_tissue", "Epithelial_cells", "Endothelial_cells"]
CELL_TYPES = ["Fibroblasts", "Stromal_and_Immune_cells", "Adipocytes", "Breast_tissue", "Epithelial_cells", "Endothelial_cells"]
ANCESTRIES = ["BLACK", "WHITE"]
# LD_NAME_MAP = {
#    "BLACK": "eur.afr",
#    "WHITE": "eur"
# }
JM_TO_JL_KEY = {"BLACK": "afr",
                "WHITE": "eur"}

final_spredixcans = {
    "outputs": [],
    "ancestry": [],
    "celltype": [],
    "subtype": []
}

# Track futures
parse_gwas_futures = []
parse_gwas_outputs = {}  # key = (ances, subtype) → future
cojo_parse_futures = []
cojo_parse_outputs = {}  # key = (ances, subtype) → future
impute_futures = []

for ances in ANCESTRIES:
    jl_ances = JM_TO_JL_KEY[ances]
    gwas_dir = f"../input/{jl_ances}_meta_sumstats/"
    gwas_fnames = os.listdir(gwas_dir)

    for subtype_gwas in gwas_fnames:
        subtype_name = os.path.splitext(subtype_gwas)[0]

        # Launch main parse
        parse_future = metaxcan_parse_gwas(
            os.path.join(gwas_dir, subtype_gwas),
            os.path.join(f"../input/model_geno_ref/{ances}.variants.txt.gz"),
            "james_gwas_brca_stype_by_ances",
            outputs=[File(f"../output/metaxcan_gwas_parse/{ances}.{subtype_name}.txt.gz")]
        )
        parse_gwas_futures.append(parse_future)
        parse_gwas_outputs[(ances, subtype_name)] = parse_future

        # Launch COJO parse
        cojo_future = metaxcan_parse_gwas(
            os.path.join(gwas_dir, subtype_gwas),
            os.path.join(f"../input/cojo_ld_ref_jcm/COJO.{ances}.variants.txt.gz"),
            "james_gwas_brca_stype_by_ances",
            outputs=[File(f"../output/metaxcan_gwas_parse/COJO.{ances}.{subtype_name}.txt.gz")]
        )
        cojo_parse_futures.append(cojo_future)
        cojo_parse_outputs[(ances, subtype_name)] = cojo_future

# Wait for all parse futures (main and cojo)
for future in parse_gwas_futures:
    future.result()
for future in cojo_parse_futures:
    future.result()

# Now safe to run filter_gwas_by_regions and schedule imputation
for (ances, subtype_name), cojo_future in cojo_parse_outputs.items():
    jl_ances = JM_TO_JL_KEY[ances] # eur or afr
    filter_gwas_by_regions(
        gwas=cojo_future.outputs[0],
        id_convert_file=f"../input/id_convert_files/{ances}.vander_gtex_model_ids_to_james_cojo_ld_ref.tsv",
        #region_file=f"../input/ld_blocks/{ances}_ld.bed",
        region_file=f"../input/ld_blocks/ma_focus_lds/grch38.{jl_ances}.loci.bed",
        inputs=cojo_future.outputs,
        out_dir="../output/by_region_cojo_gwas/"
    )

# Now launch imputation
for (ances, subtype_name), parse_future in parse_gwas_outputs.items():
    jl_ances = JM_TO_JL_KEY[ances] # eur or afr
    postprocess_pattern = f"{ances}.{subtype_name}"
    impute_future = metaxcan_impute_gwas(
        parse_future.outputs[0].filepath,
        # f"../input/ld_blocks/hg38_ld/{ances}_chr{{chr_num}}_ld.bed",
        f"../input/ld_blocks/ma_focus_lds/chr{{chr_num}}_grch38.{jl_ances}.loci.bed",
        f"../input/model_geno_ref/{ances}.chr{{chr_num}}.geno.parquet",
        f"../input/model_geno_ref/{ances}.metadata.parquet/chrom={{chr_num}}/part-0.parquet",
        postprocess_pattern,
        output_format=f"../output/metaxcan_gwas_impute_regions/{postprocess_pattern}_chr{{chr_num}}_sb{{sub_batch}}_reg{{reg}}_ff{{ff}}_by_region.txt.gz",
        output=[File(f"../output/metaxcan_gwas_imputed/{postprocess_pattern}.txt.gz")]
    )
    impute_futures.append((impute_future, ances, subtype_name))

# Wait for imputation before S-PrediXcan
for impute_future, ances, subtype_name in impute_futures:
    impute_future.result()

    for ctype in CELL_TYPES:
        spredi_future = spredixcan(
            gwas_file=impute_future.outputs[0].filepath,
            db_file=f"../input/predixcan_dbs/{ances}.gene_models.{ctype}.db",
            covar_file=f"../input/predixcan_lds/{ances}.gene_models.{ctype}.txt.gz",
            outputs=[File(f"../output/spredixcan_results/{ances}.{ctype}.{subtype_name}.txt")]
        )

        final_spredixcans["outputs"].append(spredi_future)
        final_spredixcans["ancestry"].append(ances)
        final_spredixcans["celltype"].append(ctype)
        final_spredixcans["subtype"].append(subtype_name)

# Combine spredixcans only for successful TWAS
filtered_data = {"outputs": [], "ancestry": [], "celltype": [], "subtype": []}
for future, ancestry, celltype, subtype in zip(
    final_spredixcans["outputs"],
    final_spredixcans["ancestry"],
    final_spredixcans["celltype"],
    final_spredixcans["subtype"]
):
    try:
        future.result()  # Check if future completed successfully
        filtered_data["outputs"].append(future.outputs[0])
        filtered_data["ancestry"].append(ancestry)
        filtered_data["celltype"].append(celltype)
        filtered_data["subtype"].append(subtype)
    except Exception as e:
        print(f"Skipping failed future: {e}")


combine_spredixcans(inputs=filtered_data["outputs"],
        ancestry=filtered_data["ancestry"],
        celltype=filtered_data["celltype"],
        subtype=filtered_data["subtype"])

parsl.wait_for_current_tasks()
time.sleep(20)
