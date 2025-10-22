#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
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

from options import MetaxcanParse
from metaxcan_parse import metaxcan_parse_gwas, small_metaxcan_parse_gwas
from metaxcan_impute import metaxcan_impute_gwas
from metaxcan_spredixcan import spredixcan
from cojo_prep import filter_gwas_by_regions, combine_spredixcans
from gcta_cojo import cojo_create_ma, get_model_snp_lists, get_cond_snp_list, gcta_cojo, combine_condTWAS_results


config.retries = 1
parsl.load(config)
parsl.set_stream_logger(level=logging.DEBUG)

CELL_TYPES = ["Fibroblasts", "Stromal_and_Immune_cells", "Adipocytes", "Breast_tissue", "Epithelial_cells", "Endothelial_cells"]
ANCESTRIES = ["BLACK", "WHITE"]
JM_TO_JL_KEY = {"BLACK": "afr",
                "WHITE": "eur"}

def run_cojo_ctwas(file_path, lead_var_after_before=2500000,
        clump_r2=0.65,
        sample_size=None,
        tracking_opath="../output/tracking/ancestry.subtype.celltype.cojo.condTWAS.tracking.tsv"):
    """
    Reads a tab-delimited file and processes columns: ensg_id, gene_start, gene_end, ancestry, celltype, subtype.
    Optionally selects a random sample of specified size.

    :param file_path: Path to the tab-delimited file.
    :param sample_size: Number of rows to randomly sample (default: None for full dataset).
    :param tracking_opath: Path to save the tracking dataframe.
    """
    # Read the file
    df = pd.read_csv(file_path, sep="\t", dtype={"ucsc_cytoband": str})

    # Sample if requested
    if sample_size:
        df = df.sample(n=sample_size, random_state=42)

    tracking_rows = []
    cond_results = []
    model_snp_results = []
    cojo_results = []
    parse_results = []
    impute_results = []
    spredixcan_results = []
    # Loop through each row
    cond_snp_odir = "../output/cond_snps/"
    for _, row in df.iterrows():
        ucsc_cytoband = row["ucsc_cytoband"]

        chr_num = int(ucsc_cytoband.split(':')[0])
        #try:
        #    chr_num = int(ucsc_cytoband.split(':')[0])
        #except AttributeError:
        #    chr_num = 14

        ensg_id = row["ensg_id"]
        gene_start = row["gene_start"]
        gene_end = row["gene_end"]
        ancestry = row["ancestry"]
        celltype = row["celltype"]
        subtype = row["subtype"]

        tracking_row = {
            "ensg_id": ensg_id,
            "ancestry": ancestry,
            "celltype": celltype,
            "subtype": subtype,
            "status": "successful"
        }

        print(f"Processing {ensg_id} ({gene_start}-{gene_end}), {ancestry}, {celltype}, {subtype}")

        gwas_ma = File(f"../output/cojo_ma/{ancestry}.{subtype}.ma")
        ma_res = cojo_create_ma(f"../output/metaxcan_gwas_parse/COJO.{ancestry}.{subtype}.txt.gz",
            outputs=[gwas_ma]
        )
        ma_res.result()
        ma_file = ma_res.outputs[0]

        cond_snp_oput = [File(os.path.join(cond_snp_odir, f"{ancestry}.{subtype}.{celltype}_{ensg_id}.clumped.cond.snplist"))]
        cond_snp_res = get_cond_snp_list(
            ensg_id, ancestry, subtype, celltype, chr_num,
            gene_start, gene_end, lead_var_after_before, clump_r2=clump_r2,
            ma_file=ma_file.filepath,
            out_dir=cond_snp_odir,
            outputs=cond_snp_oput#,
            #stdout=f"logs/{ancestry}.{subtype}.{celltype}_{ensg_id}.get_cond_snp_list.log",
            #stderr=f"logs/{ancestry}.{subtype}.{celltype}_{ensg_id}.get_cond_snp_list.log"
        )

        cond_results.append({
            "future": cond_snp_res,
            "output_file": cond_snp_oput[0],
            "tracking_row": tracking_row,
            "chr_num": chr_num,
            "ensg_id": ensg_id,
            "ancestry": ancestry,
            "subtype": subtype,
            "celltype": celltype,
            "gene": ensg_id,
            "gwas_ma": gwas_ma,
            "ma_future": ma_res

        })

    # After all tasks are launched, wait and check outputs
    model_snp_odir = "../output/model_snps/"
    for result in cond_results:
        result_text = result["future"].result()  # wait for task to complete
        result["ma_future"].result() # Let MA file creation complete (quick)

        if not os.path.exists(result["output_file"].filepath):
            result["tracking_row"]["status"] = result_text
            tracking_rows.append(result["tracking_row"])
            continue

        model_snp_oputs = [
            File(os.path.join(model_snp_odir, f"{result['ancestry']}.{result['subtype']}.{result['celltype']}_{result['ensg_id']}_model.ids.snplist")),
            File(os.path.join(model_snp_odir, f"{result['ancestry']}.{result['subtype']}.{result['celltype']}_{result['ensg_id']}_cojo.ids.snplist")),
        ]
        model_res = get_model_snp_lists(
            result["ensg_id"], celltype=result["celltype"],
            subtype=result["subtype"], ancestry=result["ancestry"],
            out_dir=model_snp_odir,
            outputs=model_snp_oputs
        )

        model_snp_results.append({
            "future": model_res,
            "model_snp_files": model_snp_oputs,
            "tracking_row": result["tracking_row"],
            "chr_num": result["chr_num"],
            "ensg_id": result["ensg_id"],
            "ancestry": result["ancestry"],
            "subtype": result["subtype"],
            "celltype": result["celltype"],
            "gwas_ma": result["gwas_ma"],
            "cond_snp_file": result["output_file"]
        })

    # Wait for all model SNP results, run COJO where successful
    cojo_odir = "../output/cojo_results/"
    for result in model_snp_results:
        result["future"].result()  # wait for task to complete
        for output_file in result["model_snp_files"]:
            if not os.path.exists(output_file.filepath):
                result["tracking_row"]["status"] = f"{output_file.filepath} missing"
                tracking_rows.append(result["tracking_row"])
                break
            # File(os.path.join(model_snp_odir, f"{result['ancestry']}.{result['subtype']}.{result['celltype']}_{result['ensg_id']}_model.ids.snplist")),

        cojo_prefix = os.path.join(cojo_odir, f"{result['ancestry']}.{result['subtype']}.{result['celltype']}_{result['ensg_id']}")
        cojo_oputs = [File(f"{cojo_prefix}.cma.cojo"), File(f"{cojo_prefix}.cma.cojo.parse")]

        cojo_res = gcta_cojo(result["gwas_ma"].filepath,
            result["model_snp_files"][1], result["cond_snp_file"],
            f"../input/LD_REFERENCE_PANELS/{JM_TO_JL_KEY[result['ancestry']]}",
            cojo_odir, cojo_prefix, outputs=cojo_oputs)
        cojo_results.append({
            "chr_num": result["chr_num"],
            "ensg_id": result["ensg_id"],
            "ancestry": result["ancestry"],
            "subtype": result["subtype"],
            "celltype": result["celltype"],
            "model_snp_files": result["model_snp_files"],
            "cond_snp_file": result["cond_snp_file"],
            "cojo_fut": cojo_res,
            "cojo_for_parse": cojo_oputs[1],
            "tracking_row": result["tracking_row"]
            })


    # Wait for COJO results, parse results where possible.
    parse_odir = "../output/metaxcan_gwas_parse/cojo/"
    for res in cojo_results:
        res["cojo_fut"].result() # Wait for cojo to finish
        if not os.path.exists(res["cojo_for_parse"]):
            res["tracking_row"]["status"] = f"{res['cojo_for_parse'].filepath} missing"
            tracking_rows.append(res["tracking_row"])
            next

        prefix = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
        parse_gwas_res = small_metaxcan_parse_gwas(res["cojo_for_parse"],
            f"../input/model_geno_ref/{res['ancestry']}.variants.txt.gz",
            "james_gwas_brca_stype_by_ances_cojo",
            outputs=[File(os.path.join(parse_odir, f"{prefix}.txt.gz"))]
        )

        parse_results.append({
            "chr_num": res["chr_num"],
            "ensg_id": res["ensg_id"],
            "ancestry": res["ancestry"],
            "subtype": res["subtype"],
            "celltype": res["celltype"],
            "model_snp_files": res["model_snp_files"],
            "cond_snp_file": res["cond_snp_file"],
            "parse_fut": parse_gwas_res,
            "tracking_row": res["tracking_row"],
            })

    # Parse COJO results 
    for res in parse_results:
        if res["parse_fut"].exception():
            res["tracking_row"]["status"] = f"Parsing after COJO failed"
            tracking_rows.append(res["tracking_row"])
            next

        if res["parse_fut"].done():
            parsed_gwas = res["parse_fut"].outputs[0]
            pp_pattern = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
            imputed_gwas = metaxcan_impute_gwas(parsed_gwas,
                    #f"../input/ld_blocks/hg38_ld/{res['ancestry']}_chr{{chr_num}}_ld.bed",
                    f"../input/ld_blocks/ma_focus_lds/chr{{chr_num}}_grch38.{JM_TO_JL_KEY[res['ancestry']]}.loci.bed",
                    f"../input/model_geno_ref/{res['ancestry']}.chr{{chr_num}}.geno.parquet",
                    f"../input/model_geno_ref/{res['ancestry']}.metadata.parquet/chrom={{chr_num}}/part-0.parquet",
                    pp_pattern,
                    chromosome=res["chr_num"],
                    output_format=f"../output/metaxcan_gwas_impute_regions/cojo/{pp_pattern}_chr{{chr_num}}_sb{{sub_batch}}_reg{{reg}}_ff{{ff}}_by_region.txt.gz",
                    output=[File(f"../output/metaxcan_gwas_imputed/cojo/{pp_pattern}.txt.gz")])
            impute_results.append({
            "chr_num": res["chr_num"],
            "ensg_id": res["ensg_id"],
            "ancestry": res["ancestry"],
            "subtype": res["subtype"],
            "celltype": res["celltype"],
            "model_snp_files": res["model_snp_files"],
            "cond_snp_file": res["cond_snp_file"],
            "imputed_gwas_fut": imputed_gwas,
            "tracking_row": res["tracking_row"]})

    # Impute parsed results
    for res in impute_results:
        if res["imputed_gwas_fut"].exception():
            res["tracking_row"]["status"] = f"Imputation after COJO failed"
            tracking_rows.append(res["tracking_row"])
            next

        if res["imputed_gwas_fut"].done():
            pattern = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
            imputed_gwas = res["imputed_gwas_fut"].outputs[0]
            spredixcan_res = spredixcan(gwas_file=imputed_gwas,
                db_file=f"../input/predixcan_dbs/{res['ancestry']}.gene_models.{res['celltype']}.db",
                covar_file=f"../input/predixcan_lds/{res['ancestry']}.gene_models.{res['celltype']}.txt.gz",
                outputs=[File(f"../output/spredixcan_results/cojo/{pattern}.txt")]
                )

            spredixcan_results.append({
            "chr_num": res["chr_num"],
            "ensg_id": res["ensg_id"],
            "ancestry": res["ancestry"],
            "subtype": res["subtype"],
            "celltype": res["celltype"],
            "spredixcan_fut": spredixcan_res,
            "tracking_row": res["tracking_row"]})

    # Run spredixcan on imputed results
    for res in spredixcan_results:
        if res["spredixcan_fut"].exception():
            res["tracking_row"]["status"] = f"Spredixcan after COJO failed"
            tracking_rows.append(res["tracking_row"])
            next

        if res["spredixcan_fut"].done():
            tracking_rows.append(res["tracking_row"])
            pass

    # After loop: write tracking rows to file
    if tracking_rows:
        tracking_df = pd.DataFrame(tracking_rows)
        os.makedirs(os.path.dirname(tracking_opath), exist_ok=True)
        tracking_df.to_csv(tracking_opath, sep="\t", index=False)
        print(f"Tracking information written to {tracking_opath}")
    parsl.wait_for_current_tasks()

def old_run_cojo_ctwas_single_cond(file_path,
        tracking_opath="../output/tracking/ancestry.subtype.celltype.cojo.condTWAS.single_snp.tracking.tsv"):
    """
    Processes a simplified input file with columns: subtype, ensg_id, variant_id, chr_num.
    Runs COJO conditioning on a single variant per gene, across both ancestries.

    :param file_path: Path to the tab-delimited input file.
    :param tracking_opath: Path to save the tracking dataframe.
    """
    df = pd.read_csv(file_path, sep="\t")

    tracking_rows = []
    model_snp_results = []
    cojo_results = []
    parse_results = []
    impute_results = []
    spredixcan_results = []

    ancestries = ['BLACK', 'WHITE']
    model_snp_odir = "../output/model_snps_single_cond/"
    cojo_odir = "../output/cojo_results_single_cond/"
    parse_odir = "../output/metaxcan_gwas_parse/cojo_single_cond/"
    impute_odir = "../output/metaxcan_gwas_imputed/cojo_single_cond/"
    spredixcan_odir = "../output/spredixcan_results/cojo_single_cond/"

    for _, row in df.iterrows():
        subtype = row["subtype"]
        ensg_id = row["ensg_id"]
        variant_id = row["variant_id"]
        chr_num = int(row["chr_num"])

        for ancestry in ancestries:
            for celltype in CELL_TYPES:  # CELLTYPES assumed defined globally
                tracking_row = {
                    "ensg_id": ensg_id,
                    "ancestry": ancestry,
                    "celltype": celltype,
                    "subtype": subtype,
                    "status": "successful"
                }

                print(f"Processing {ensg_id}, variant {variant_id}, chr {chr_num}, {ancestry}, {celltype}, {subtype}")

                gwas_ma = File(f"../output/cojo_ma/{ancestry}.{subtype}.ma")
                ma_res = cojo_create_ma(f"../output/metaxcan_gwas_parse/COJO.{ancestry}.{subtype}.txt.gz",
                    outputs=[gwas_ma]
                )
                ma_res.result()

                model_snp_oputs = [
                    File(os.path.join(model_snp_odir, f"{ancestry}.{subtype}.{celltype}_{ensg_id}_model.ids.snplist")),
                    File(os.path.join(model_snp_odir, f"{ancestry}.{subtype}.{celltype}_{ensg_id}_cojo.ids.snplist")),
                ]
                model_res = get_model_snp_lists(
                    ensg_id, celltype=celltype,
                    subtype=subtype, ancestry=ancestry,
                    out_dir=model_snp_odir,
                    outputs=model_snp_oputs
                )

                model_snp_results.append({
                    "future": model_res,
                    "model_snp_files": model_snp_oputs,
                    "tracking_row": tracking_row,
                    "ensg_id": ensg_id,
                    "ancestry": ancestry,
                    "subtype": subtype,
                    "celltype": celltype,
                    "gwas_ma": gwas_ma,
                    "single_snp": variant_id,
                    "chr_num": chr_num
                })

    for result in model_snp_results:
        result["future"].result()

        # Skip if any output file is missing
        missing = False
        for output_file in result["model_snp_files"]:
            if not os.path.exists(output_file.filepath):
                result["tracking_row"]["status"] = f"{output_file.filepath} missing"
                tracking_rows.append(result["tracking_row"])
                missing = True
                break
        if missing:
            continue

        cojo_prefix = os.path.join(cojo_odir, f"{result['ancestry']}.{result['subtype']}.{result['celltype']}_{result['ensg_id']}")
        cojo_oputs = [File(f"{cojo_prefix}.cma.cojo"), File(f"{cojo_prefix}.cma.cojo.parse")]

        cojo_res = gcta_cojo(result["gwas_ma"].filepath,
            result["model_snp_files"][1], "",  # dummy cond_snp_list_path
            f"../input/LD_REFERENCE_PANELS/{JM_TO_JL_KEY[result['ancestry']]}",
            cojo_odir, cojo_prefix,
            single_snp=result["single_snp"],
            outputs=cojo_oputs)

        cojo_results.append({
            "ensg_id": result["ensg_id"],
            "ancestry": result["ancestry"],
            "subtype": result["subtype"],
            "celltype": result["celltype"],
            "chr_num": result["chr_num"],
            "model_snp_files": result["model_snp_files"],
            "cojo_fut": cojo_res,
            "cojo_for_parse": cojo_oputs[1],
            "tracking_row": result["tracking_row"]
        })

    for res in cojo_results:
        res["cojo_fut"].result()
        if not os.path.exists(res["cojo_for_parse"]):
            res["tracking_row"]["status"] = f"{res['cojo_for_parse'].filepath} missing"
            tracking_rows.append(res["tracking_row"])
            continue

        prefix = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
        parse_gwas_res = small_metaxcan_parse_gwas(res["cojo_for_parse"],
            f"../input/model_geno_ref/{res['ancestry']}.variants.txt.gz",
            "james_gwas_brca_stype_by_ances_cojo",
            outputs=[File(os.path.join(parse_odir, f"{prefix}.txt.gz"))]
        )

        parse_results.append({
            "ensg_id": res["ensg_id"],
            "ancestry": res["ancestry"],
            "subtype": res["subtype"],
            "celltype": res["celltype"],
            "chr_num": res["chr_num"],
            "model_snp_files": res["model_snp_files"],
            "parse_fut": parse_gwas_res,
            "tracking_row": res["tracking_row"]
        })

    for res in parse_results:
        if res["parse_fut"].exception():
            res["tracking_row"]["status"] = f"Parsing after COJO failed"
            tracking_rows.append(res["tracking_row"])
            continue

        if res["parse_fut"].done():
            parsed_gwas = res["parse_fut"].outputs[0]
            pp_pattern = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
            imputed_gwas = metaxcan_impute_gwas(parsed_gwas,
                    f"../input/ld_blocks/ma_focus_lds/chr{res['chr_num']}_grch38.{JM_TO_JL_KEY[res['ancestry']]}.loci.bed",
                    f"../input/model_geno_ref/{res['ancestry']}.chr{res['chr_num']}.geno.parquet",
                    f"../input/model_geno_ref/{res['ancestry']}.metadata.parquet/chrom={res['chr_num']}/part-0.parquet",
                    pp_pattern,
                    chromosome=res['chr_num'],
                    output_format=f"../output/metaxcan_gwas_impute_regions/cojo_single_cond/{pp_pattern}_chr{{chr_num}}_sb{{sub_batch}}_reg{{reg}}_ff{{ff}}_by_region.txt.gz",
                    output=[File(os.path.join(impute_odir, f"{pp_pattern}.txt.gz"))])

            impute_results.append({
                "ensg_id": res["ensg_id"], "ancestry": res["ancestry"],
                "subtype": res["subtype"],
                "celltype": res["celltype"],
                "chr_num": res["chr_num"],
                "imputed_gwas_fut": imputed_gwas,
                "tracking_row": res["tracking_row"]
            })

    for res in impute_results:
        if res["imputed_gwas_fut"].exception():
            res["tracking_row"]["status"] = f"Imputation after COJO failed"
            tracking_rows.append(res["tracking_row"])
            continue

        if res["imputed_gwas_fut"].done():
            pattern = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
            imputed_gwas = res["imputed_gwas_fut"].outputs[0]
            spredixcan_res = spredixcan(gwas_file=imputed_gwas,
                db_file=f"../input/predixcan_dbs/{res['ancestry']}.gene_models.{res['celltype']}.db",
                covar_file=f"../input/predixcan_lds/{res['ancestry']}.gene_models.{res['celltype']}.txt.gz",
                outputs=[File(os.path.join(spredixcan_odir, f"{pattern}.txt"))]
            )

            spredixcan_results.append({
                "ensg_id": res["ensg_id"],
                "ancestry": res["ancestry"],
                "subtype": res["subtype"],
                "celltype": res["celltype"],
                "chr_num": res["chr_num"],
                "spredixcan_fut": spredixcan_res,
                "tracking_row": res["tracking_row"]
            })

    for res in spredixcan_results:
        if res["spredixcan_fut"].exception():
            res["tracking_row"]["status"] = f"Spredixcan after COJO failed"
        tracking_rows.append(res["tracking_row"])

    if tracking_rows:
        tracking_df = pd.DataFrame(tracking_rows)
        os.makedirs(os.path.dirname(tracking_opath), exist_ok=True)
        tracking_df.to_csv(tracking_opath, sep="\t", index=False)
        print(f"Tracking information written to {tracking_opath}")
    parsl.wait_for_current_tasks()

def run_cojo_ctwas_single_cond(file_path,
        tracking_opath="../output/tracking/ancestry.subtype.celltype.cojo.condTWAS.single_snp.tracking.tsv"):
    """
    Processes a simplified input file with columns: subtype, ensg_id, variant_id, chr_num.
    Runs COJO conditioning on a single variant per gene, across both ancestries.

    :param file_path: Path to the tab-delimited input file.
    :param tracking_opath: Path to save the tracking dataframe.
    """
    df = pd.read_csv(file_path, sep="\t")

    tracking_rows = []
    model_snp_results = []
    cojo_results = []
    parse_results = []
    impute_results = []
    spredixcan_results = []

    ancestries = ['BLACK', 'WHITE']
    model_snp_odir = "../output/model_snps_single_cond/"
    cojo_odir = "../output/cojo_results_single_cond/"
    parse_odir = "../output/metaxcan_gwas_parse/cojo_single_cond/"
    impute_odir = "../output/metaxcan_gwas_imputed/cojo_single_cond/"
    spredixcan_odir = "../output/spredixcan_results/cojo_single_cond/"

    for _, row in df.iterrows():
        subtype = row["subtype"]
        ensg_id = row["ensg_id"]
        variant_id = row["variant_id"]
        chr_num = int(row["chr_num"])

        for ancestry in ancestries:
            for celltype in CELL_TYPES:  # CELLTYPES assumed defined globally
                tracking_row = {
                    "ensg_id": ensg_id,
                    "ancestry": ancestry,
                    "celltype": celltype,
                    "subtype": subtype,
                    "status": "successful"
                }

                print(f"Processing {ensg_id}, variant {variant_id}, chr {chr_num}, {ancestry}, {celltype}, {subtype}")

                gwas_ma = File(f"../output/cojo_ma/{ancestry}.{subtype}.ma")
                ma_res = cojo_create_ma(f"../output/metaxcan_gwas_parse/COJO.{ancestry}.{subtype}.txt.gz",
                    outputs=[gwas_ma]
                )
                ma_res.result()

                model_snp_oputs = [
                    File(os.path.join(model_snp_odir, f"{ancestry}.{subtype}.{celltype}_{ensg_id}_model.ids.snplist")),
                    File(os.path.join(model_snp_odir, f"{ancestry}.{subtype}.{celltype}_{ensg_id}_cojo.ids.snplist")),
                ]
                model_res = get_model_snp_lists(
                    ensg_id, celltype=celltype,
                    subtype=subtype, ancestry=ancestry,
                    out_dir=model_snp_odir,
                    outputs=model_snp_oputs
                )

                model_snp_results.append({
                    "future": model_res,
                    "model_snp_files": model_snp_oputs,
                    "tracking_row": tracking_row,
                    "ensg_id": ensg_id,
                    "ancestry": ancestry,
                    "subtype": subtype,
                    "celltype": celltype,
                    "gwas_ma": gwas_ma,
                    "single_snp": variant_id,
                    "chr_num": chr_num
                })

    # Wait for model SNP results
    for result in model_snp_results:
        try:
            result["future"].result()
        except Exception as e:
            result["tracking_row"]["status"] = f"Model SNP list generation failed: {str(e)}"
            tracking_rows.append(result["tracking_row"])
            continue

        # Skip if any output file is missing
        missing = False
        for output_file in result["model_snp_files"]:
            if not os.path.exists(output_file.filepath):
                result["tracking_row"]["status"] = f"{output_file.filepath} missing"
                tracking_rows.append(result["tracking_row"])
                missing = True
                break
        if missing:
            continue

        cojo_prefix = os.path.join(cojo_odir, f"{result['ancestry']}.{result['subtype']}.{result['celltype']}_{result['ensg_id']}")
        cojo_oputs = [File(f"{cojo_prefix}.cma.cojo"), File(f"{cojo_prefix}.cma.cojo.parse")]

        cojo_res = gcta_cojo(result["gwas_ma"].filepath,
            result["model_snp_files"][1], "",  # dummy cond_snp_list_path
            f"../input/LD_REFERENCE_PANELS/{JM_TO_JL_KEY[result['ancestry']]}",
            cojo_odir, cojo_prefix,
            single_snp=result["single_snp"],
            outputs=cojo_oputs)

        cojo_results.append({
            "ensg_id": result["ensg_id"],
            "ancestry": result["ancestry"],
            "subtype": result["subtype"],
            "celltype": result["celltype"],
            "chr_num": result["chr_num"],
            "model_snp_files": result["model_snp_files"],
            "cojo_fut": cojo_res,
            "cojo_for_parse": cojo_oputs[1],
            "tracking_row": result["tracking_row"]
        })

    # Wait for COJO results
    for res in cojo_results:
        try:
            res["cojo_fut"].result()
        except Exception as e:
            res["tracking_row"]["status"] = f"COJO failed: {str(e)}"
            tracking_rows.append(res["tracking_row"])
            continue
        
        if not os.path.exists(res["cojo_for_parse"].filepath):
            res["tracking_row"]["status"] = f"{res['cojo_for_parse'].filepath} missing"
            tracking_rows.append(res["tracking_row"])
            continue

        prefix = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
        parse_gwas_res = small_metaxcan_parse_gwas(res["cojo_for_parse"],
            f"../input/model_geno_ref/{res['ancestry']}.variants.txt.gz",
            "james_gwas_brca_stype_by_ances_cojo",
            outputs=[File(os.path.join(parse_odir, f"{prefix}.txt.gz"))]
        )

        parse_results.append({
            "ensg_id": res["ensg_id"],
            "ancestry": res["ancestry"],
            "subtype": res["subtype"],
            "celltype": res["celltype"],
            "chr_num": res["chr_num"],
            "model_snp_files": res["model_snp_files"],
            "parse_fut": parse_gwas_res,
            "tracking_row": res["tracking_row"]
        })

    # Parse COJO results
    for res in parse_results:
        try:
            if res["parse_fut"].exception():
                res["tracking_row"]["status"] = f"Parsing after COJO failed: {res['parse_fut'].exception()}"
                tracking_rows.append(res["tracking_row"])
                continue

            if res["parse_fut"].done():
                parsed_gwas = res["parse_fut"].outputs[0]
                pp_pattern = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
                imputed_gwas = metaxcan_impute_gwas(parsed_gwas,
                        f"../input/ld_blocks/ma_focus_lds/chr{res['chr_num']}_grch38.{JM_TO_JL_KEY[res['ancestry']]}.loci.bed",
                        f"../input/model_geno_ref/{res['ancestry']}.chr{res['chr_num']}.geno.parquet",
                        f"../input/model_geno_ref/{res['ancestry']}.metadata.parquet/chrom={res['chr_num']}/part-0.parquet",
                        pp_pattern,
                        chromosome=res['chr_num'],
                        output_format=f"../output/metaxcan_gwas_impute_regions/cojo_single_cond/{pp_pattern}_chr{{chr_num}}_sb{{sub_batch}}_reg{{reg}}_ff{{ff}}_by_region.txt.gz",
                        output=[File(os.path.join(impute_odir, f"{pp_pattern}.txt.gz"))])

                impute_results.append({
                    "ensg_id": res["ensg_id"], 
                    "ancestry": res["ancestry"],
                    "subtype": res["subtype"],
                    "celltype": res["celltype"],
                    "chr_num": res["chr_num"],
                    "imputed_gwas_fut": imputed_gwas,
                    "tracking_row": res["tracking_row"]
                })
        except Exception as e:
            res["tracking_row"]["status"] = f"Parsing after COJO failed with exception: {str(e)}"
            tracking_rows.append(res["tracking_row"])
            continue

    # Impute parsed results
    for res in impute_results:
        try:
            if res["imputed_gwas_fut"].exception():
                res["tracking_row"]["status"] = f"Imputation after COJO failed: {res['imputed_gwas_fut'].exception()}"
                tracking_rows.append(res["tracking_row"])
                continue

            if res["imputed_gwas_fut"].done():
                pattern = f"{res['ancestry']}.{res['subtype']}.{res['celltype']}_{res['ensg_id']}"
                imputed_gwas = res["imputed_gwas_fut"].outputs[0]
                spredixcan_res = spredixcan(gwas_file=imputed_gwas,
                    db_file=f"../input/predixcan_dbs/{res['ancestry']}.gene_models.{res['celltype']}.db",
                    covar_file=f"../input/predixcan_lds/{res['ancestry']}.gene_models.{res['celltype']}.txt.gz",
                    outputs=[File(os.path.join(spredixcan_odir, f"{pattern}.txt"))]
                )

                spredixcan_results.append({
                    "ensg_id": res["ensg_id"],
                    "ancestry": res["ancestry"],
                    "subtype": res["subtype"],
                    "celltype": res["celltype"],
                    "chr_num": res["chr_num"],
                    "spredixcan_fut": spredixcan_res,
                    "tracking_row": res["tracking_row"]
                })
        except Exception as e:
            res["tracking_row"]["status"] = f"Imputation after COJO failed with exception: {str(e)}"
            tracking_rows.append(res["tracking_row"])
            continue

    # Run spredixcan on imputed results
    for res in spredixcan_results:
        try:
            if res["spredixcan_fut"].exception():
                res["tracking_row"]["status"] = f"Spredixcan after COJO failed: {res['spredixcan_fut'].exception()}"
                tracking_rows.append(res["tracking_row"])
                continue
            
            if res["spredixcan_fut"].done():
                tracking_rows.append(res["tracking_row"])
        except Exception as e:
            res["tracking_row"]["status"] = f"Spredixcan after COJO failed with exception: {str(e)}"
            tracking_rows.append(res["tracking_row"])
            continue

    # After loop: write tracking rows to file
    if tracking_rows:
        tracking_df = pd.DataFrame(tracking_rows)
        os.makedirs(os.path.dirname(tracking_opath), exist_ok=True)
        tracking_df.to_csv(tracking_opath, sep="\t", index=False)
        print(f"Tracking information written to {tracking_opath}")
    else:
        print("No tracking rows to write - all tasks may have completed successfully")
    
    parsl.wait_for_current_tasks()



# Example usage:
try:
    # run_cojo_ctwas("../input/james_parsed_twas_output/james_cojo_input_table.tsv", sample_size=None)
    # combine_condTWAS_results("../output/spredixcan_results/cojo/", outputs=[File("../output/spredixcan_combined/cojo/james_vars_cojo_condTWAS_summary.tsv")])
    parsl.wait_for_current_tasks()
    # run_cojo_ctwas_single_cond("../input/jcm_mod_single_cond_input.tsv")
    # combine_condTWAS_results("../output/spredixcan_results/cojo_single_cond/",
    #         tracker_path="../output/tracking/ancestry.subtype.celltype.cojo.condTWAS.single_snp.tracking.tsv",
    #         outputs=[File("../output/spredixcan_combined/cojo/single_snp_james_vars_cojo_condTWAS_summary.tsv")])
    run_cojo_ctwas_single_cond("../input/jcm_10_15_2025_yijia_added_input.tsv",
                                tracking_opath="../output/tracking/10_15_2025.ancestry.subtype.celltype.cojo.condTWAS.single_snp.tracking.tsv")
    parsl.wait_for_current_tasks()
    combine_condTWAS_results("../output/spredixcan_results/cojo_single_cond/",
            tracker_path="../output/tracking/10_15_2025.ancestry.subtype.celltype.cojo.condTWAS.single_snp.tracking.tsv",
            outputs=[File("../output/spredixcan_combined/cojo/10_15_2025_single_snp_james_vars_cojo_condTWAS_summary.tsv")])
    parsl.wait_for_current_tasks()
finally:
    # parsl.wait_for_current_tasks()
    parsl.clear()

