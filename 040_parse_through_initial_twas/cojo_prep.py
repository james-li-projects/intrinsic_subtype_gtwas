#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import parsl
import subprocess
from subprocess import call
from parsl.app.app import bash_app, python_app, join_app
from parsl.data_provider.files import File
from parsl.dataflow.futures import AppFuture

from config import id_for_memo_File

@python_app(cache=True, executors=["parser"])
def filter_gwas_by_regions(
    gwas, region_file, id_convert_file, out_dir=".", inputs=[], outputs=[],
    stdout=parsl.AUTO_LOGNAME, stderr=parsl.AUTO_LOGNAME
):
    """
    Filters GWAS data based on genomic regions and writes separate files per
    region.
    
    Parameters:
        gwas (str): Path to the GWAS file.
        region_file (str): Path to the region file.
        id_convert_file (str): Path to the ID conversion file.
        out_dir (str, optional): Directory to save output files. Default is
            current directory.
    """
    try:
        import pandas as pd
        import os
        import re

        # Read the GWAS file
        gwas_df = pd.read_csv(gwas, sep="\t")

        # Extract numerical part from chromosome column
        gwas_df["chromosome"] = (
            gwas_df["chromosome"].astype(str).str.extract(r'(\d+)')[0]
        )
        gwas_df["chromosome"] = pd.to_numeric(
            gwas_df["chromosome"], errors='coerce'
        )

        # Read the ID conversion file
        id_convert_df = pd.read_csv(id_convert_file, sep="\t")

        # Join gwas_df to id_convert_df
        gwas_df = gwas_df.merge(id_convert_df, left_on="panel_variant_id", right_on="snp_bim2", how="left")

        # Filter to keep only non-NA values of snp_bim1
        gwas_df = gwas_df[~gwas_df["snp_bim1"].isna()]

        # Remove the snp_bim1 column
        gwas_df = gwas_df.drop(columns=["snp_bim1"])

        # Read the region file
        region_df = pd.read_table(region_file, delim_whitespace=True)

        # Ensure output directory exists
        os.makedirs(out_dir, exist_ok=True)

        # Extract filename without .txt.gz extension
        gwas_name = re.sub(r'\.txt\.gz$', '', os.path.basename(gwas))

        # Loop through each region and filter GWAS data
        for _, region in region_df.iterrows():
            chr_num, start, stop = region["chr"], region["start"], region["stop"]

            filtered_df = gwas_df[
                (gwas_df["chromosome"] == chr_num) &
                (gwas_df["position"] >= start) &
                (gwas_df["position"] <= stop)
            ]

            output_file = os.path.join(
                out_dir, f"{gwas_name}_chr{chr_num}_{start}_{stop}.txt.gz"
            )

            if not os.path.exists(output_file):
                filtered_df.to_csv(output_file, sep="\t", index=False)
                print(f"Saved {len(filtered_df)} entries to {output_file}")
            else:
                print(f"File already exists: {output_file} skipping write.")
    except Exception as e:
        print(f"Error: {e}")


@python_app(executors=["parser"])
def combine_spredixcans(ancestry, celltype, subtype, 
                         cytoband=File("../input/hg38_cytoBand.txt.gz"), 
                         gencode=File("../input/gencode_v26_all.txt"),
                         output_name="../output/spredixcan_combined/ancestry.subtype.celltype.TWAS.results.txt",
                         inputs=[], outputs=[],
                         stdout=parsl.AUTO_LOGNAME, stderr=parsl.AUTO_LOGNAME):

    import os
    import pandas as pd  

    outputs = inputs # leftover from when I wanted this function to fail if spredixcan input failed

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_name), exist_ok=True)
    
    # Load GENCODE file
    gencode_df = pd.read_csv(gencode.filepath, sep="\t")
    gencode_df["chromosome"] = gencode_df["chromosome"].str.extract(r'chr(\d+)')
    gencode_df["chromosome"] = pd.to_numeric(gencode_df["chromosome"], errors='coerce').astype('Int64')

    # Load Cytoband file
    cytoband_df = pd.read_csv(cytoband.filepath, sep="\t", header=None, 
                              names=["Chr", "cyto_start", "cyto_end", "ucsc_cytoband", "info"])
    cytoband_df["Chr"] = cytoband_df["Chr"].str.extract(r'chr(\d+)')
    cytoband_df["Chr"] = pd.to_numeric(cytoband_df["Chr"], errors='coerce').astype('Int64')
    cytoband_df = cytoband_df[cytoband_df["Chr"].notna()] # Don't want chromosomes without numbers in them, some weird stuff in cyto

    # Remake ucsc_cytoband column
    cytoband_df["ucsc_cytoband"] = cytoband_df.apply(lambda row: f"{row['Chr']}:{row['ucsc_cytoband']}", axis=1)

    all_data = []

    for output, anc, cell, sub in zip(outputs, ancestry, celltype, subtype):
        try:
            filepath = output.filepath
            df = pd.read_csv(filepath, usecols=["gene", "zscore", "pvalue"])
            df["ancestry"] = anc
            df["celltype"] = cell
            df["subtype"] = sub
            all_data.append(df)
        except Exception as e:
            print(f"Processing failed for output: Error: {e}")
    
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)

        final_df = final_df.merge(
            gencode_df[["chromosome", "gene_name", "start_location", "end_location", "gene_type", "gene_id"]], 
            left_on="gene", right_on="gene_id", how="left"
        ).drop(columns=["gene_id"])

        def find_cytoband(row): # See where the start location of the gene falls for cytobands
            match = cytoband_df[
                (cytoband_df["Chr"] == row["chromosome"]) &
                (cytoband_df["cyto_start"] <= row["start_location"]) &
                (cytoband_df["cyto_end"] >= row["start_location"]) 
            ]

            return match["ucsc_cytoband"].values[0] if not match.empty else None

        final_df["ucsc_cytoband"] = final_df.apply(find_cytoband, axis=1)

        final_df = final_df.rename(columns={
            "gene": "ensg_id",
            "start_location": "gene_start",
            "end_location": "gene_end"
        })

        final_df = final_df[
            ["ucsc_cytoband", "gene_name", "ensg_id", "gene_type", "gene_start", "gene_end", "zscore", "pvalue", "ancestry", "celltype", "subtype"]
        ]

        final_df.to_csv(output_name, sep="\t", index=False)
        print(f"Saved results to: {output_name}")

