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

@python_app(cache=True, executors=["parser"])
def cojo_create_ma(input_file, 
                   column_mappings={"panel_variant_id": "SNP", "effect_allele": "A1", 
                                    "non_effect_allele": "A2", "frequency": "freq", 
                                    "effect_size": "b", "standard_error": "se", 
                                    "pvalue": "p", "sample_size": "N"},
                   inputs=[], outputs=[],
                   stdout=parsl.AUTO_LOGNAME,
                   stderr=parsl.AUTO_LOGNAME):
    import os
    import pandas as pd
    
    try:
        # Check if output file already exists
        output_filepath = outputs[0].filepath
        if os.path.exists(output_filepath):
            print(f"Output file {output_filepath} already exists. Skipping processing.")
            return

        # Read the tab-separated file
        df = pd.read_csv(input_file, sep='\t')

        # Rename the columns based on user input
        df.rename(columns=column_mappings, inplace=True)

        # Ensure the required columns are present
        # A1 is the effect allele
        required_columns = ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")

        # Create output directory if it does not exist
        output_dir = os.path.dirname(output_filepath)
        os.makedirs(output_dir, exist_ok=True)

        # Save the output as a space-delimited file
        df.to_csv(output_filepath, sep=' ', index=False, columns=required_columns)
        print(f"File has been processed and saved to {output_filepath}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return


@python_app(cache=True, executors=["express"])
def get_model_snp_lists(gene,
        weight_db_id_col="varID",
        alt_id_col="cojo_snp_id",
        subtype=None,
        ancestry=None,
        celltype=None,
        db_dir="../input/predixcan_dbs/",
        out_dir=".",
        by_region_dir="../output/by_region_cojo_gwas",
        region_id_want="snp_bim2",
        inputs=[], outputs=[],
        stdout=parsl.AUTO_LOGNAME,
        stderr=parsl.AUTO_LOGNAME):
    import os
    import sqlite3
    import pandas as pd
    import glob

    # Construct the weight_db path
    weight_db = os.path.join(db_dir, f"{ancestry}.gene_models.{celltype}.db")

    # Connect to the SQLite database
    conn = sqlite3.connect(weight_db)

    # Define the query with additional columns
    gwas_column = f'"in.GWAS.{subtype}"' if subtype else 'NULL as "in.GWAS."'
    query = f"""
    SELECT {weight_db_id_col}, {alt_id_col},
           "ancestry.ld.chr", "ancestry.ld.start", "ancestry.ld.stop",
           {gwas_column}
    FROM weights
    WHERE gene = ?
    """

    # Execute the query and store the result in a DataFrame
    gene_weight_df = pd.read_sql_query(query, conn, params=(gene,))

    # Create a new "region_id" column
    gene_weight_df["region_id"] = (
        "chr" + gene_weight_df["ancestry.ld.chr"].astype(str) + "_" +
        gene_weight_df["ancestry.ld.start"].astype(str) + "_" +
        gene_weight_df["ancestry.ld.stop"].astype(str)
    )

    # Create a list of unique region IDs
    unique_region_ids = gene_weight_df["region_id"].unique().tolist()

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Create separate DataFrames for writing
    model_ids_df = gene_weight_df[[weight_db_id_col]]

    # Filter cojo_ids_df where GWAS column == 1 and drop NA
    if subtype:
        cojo_ids_df = gene_weight_df[gene_weight_df[f"in.GWAS.{subtype}"] == 1][[alt_id_col]].dropna()
    else:
        cojo_ids_df = pd.DataFrame(columns=[alt_id_col])

    # Store number of rows
    num_model_ids = model_ids_df.shape[0]
    num_cojo_ids = cojo_ids_df.shape[0]

    # If there are fewer cojo IDs than model IDs, bind from region files
    proxy_cojo_ids_df = pd.DataFrame()
    use_proxy = num_model_ids > num_cojo_ids and ancestry and subtype
    if use_proxy:
        dfs = []
        for region_id in unique_region_ids:
            file_path = os.path.join(by_region_dir, f"COJO.{ancestry}.{subtype}_{region_id}.txt.gz")
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep='\t')
                if region_id_want in df.columns:
                    dfs.append(df[[region_id_want]].dropna())
        if dfs:
            proxy_cojo_ids_df = pd.concat(dfs, ignore_index=True)
            proxy_cojo_ids_df.columns = [alt_id_col]  # Rename column to match cojo_ids_df
            cojo_ids_df = pd.concat([cojo_ids_df, proxy_cojo_ids_df], ignore_index=True)
            num_cojo_ids = cojo_ids_df.shape[0]  # Update count after binding

    # Define file paths
    base_name = f"{ancestry}.{subtype}.{celltype}_{gene}"
    file1 = os.path.join(out_dir, f"{base_name}_model.ids.snplist")
    file2 = os.path.join(out_dir, f"{base_name}_cojo.ids.snplist")

    # Write DataFrames to two headerless files only if they don't exist
    if not os.path.exists(file1) and num_model_ids > 0:
        model_ids_df.to_csv(file1, index=False, header=False, sep='\t')
    if not os.path.exists(file2) and num_cojo_ids > 0:
        cojo_ids_df.to_csv(file2, index=False, header=False, sep='\t')

    # Close the database connection
    conn.close()
    return

@python_app(cache=True, executors=["imputer"])
def get_cond_snp_list(
    gene, ancestry, subtype, celltype,
    chr_num, gene_start, gene_end, lead_var_after_before,
    clump_r2,
    ma_file,
    method="gcta",
    clump_kb=10000,
    clump_p1=1e0,
    clump_snp_field="P",
    clump_field="ID",  # Default clump field is "ID"
    lead_var_dir="../input/cojo_input_lead_variants/",
    out_dir=".",
    inputs=[], outputs=[],
    stdout=parsl.AUTO_LOGNAME,
    stderr=parsl.AUTO_LOGNAME
):

    import pandas as pd
    import os
    import parsl
    import subprocess

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Mapping ancestry to population label
    JM_TO_JL_KEY = {"BLACK": "afr", "WHITE": "eur"}
    ancestry_james = JM_TO_JL_KEY[ancestry]

    # Build path to the lead variants file
    lead_var_file = os.path.join(lead_var_dir, ancestry_james, f"{subtype}.tsv")

    # Read the lead variants dataframe
    lead_vars_df = pd.read_csv(lead_var_file, sep='\t')

    # Filter for rows where CHR == chr_num and POS is within adjusted gene range
    lead_vars_df = lead_vars_df[
        (lead_vars_df['CHR'] == chr_num) &
        (lead_vars_df['POS'] >= gene_start - lead_var_after_before) &
        (lead_vars_df['POS'] <= gene_end + lead_var_after_before)
    ]

    # If no rows remain, return a descriptive message
    if lead_vars_df.empty:
        return (
            f"no lead variants in {chr_num}_{gene_start}-{gene_end} gene_start/end -/+ {lead_var_after_before}"
        )

    # Write the file with header columns "ID" and "P"
    output_filename = f"{ancestry}.{subtype}.{celltype}_{gene}.initial.cond.snplist"
    output_path = os.path.join(out_dir, output_filename)
    lead_vars_df[['ID', 'P']].to_csv(output_path, index=False, header=True, sep='\t')

    # Variables for PLINK
    bfile = os.path.join("../input/LD_REFERENCE_PANELS", ancestry_james)
    cond_snp_list = output_path
    plink_clump_prefix = os.path.join(out_dir, f"{ancestry}.{subtype}.{celltype}_{gene}")
    plink_output = f"{plink_clump_prefix}.clumped"
    final_output = os.path.join(out_dir, f"{ancestry}.{subtype}.{celltype}_{gene}.clumped.cond.snplist")

    gcta_out = plink_clump_prefix + ".jma.cojo"

    # If there's only one lead variant, skip PLINK and just write it out
    if len(lead_vars_df) == 1:
        lead_vars_df[['ID']].to_csv(final_output, index=False, header=False)
        return f"{ancestry}.{subtype}.{celltype}_{gene} only one lead variant; written without PLINK."

    # Define the PLINK bash command with memory set to 4000
    plink_command = f"""
    plink \
     --bfile {bfile} \
     --clump {cond_snp_list} \
     --clump-r2 {clump_r2} \
     --clump-kb {clump_kb} \
     --clump-p1 {clump_p1} \
     --memory 12000 \
     --threads 1 \
     --clump-snp-field {clump_snp_field} \
     --clump-field {clump_field} \
     --out {plink_clump_prefix}
    """

    gcta_command = f"""
    cat {cond_snp_list} | cut -f 1 | grep -v 'ID' > {cond_snp_list}.slct
    gcta-1.94.1 \
     --bfile {bfile} \
     --cojo-file {ma_file} \
     --extract {cond_snp_list}.slct \
     --cojo-slct \
     --out {plink_clump_prefix}
    """

    if method == "gcta":
        subprocess.run(gcta_command, shell=True, executable="/bin/bash").check_returncode()
        want_output = gcta_out
    elif method == "plink":
        subprocess.run(plink_command, shell=True, executable="/bin/bash").check_returncode()
        want_output = plink_output

    # Check if plink_output exists using Python after running the command
    if os.path.exists(want_output):
        # If output exists, process it
        select_df = pd.read_csv(want_output, delim_whitespace=True)
        select_df[['SNP']].to_csv(final_output, index=False, header=False)
    else:
        # If output does not exist, write the SNP with the minimum P value from lead_vars_df
        min_p_snp = lead_vars_df.loc[lead_vars_df['P'].idxmin(), ['ID']]
        min_p_snp.to_csv(final_output, index=False, header=False)

    return(f"{ancestry}.{subtype}.{celltype}_{gene} condition list succesfully acquired.")


@python_app(cache=False, executors=["parser", "imputer"])
def combine_condTWAS_results(input_folder,
                   tracker_path="../output/tracking/ancestry.subtype.celltype.cojo.condTWAS.tracking.tsv",
                   inputs=[], outputs=[],
                   stdout=parsl.AUTO_LOGNAME,
                   stderr=parsl.AUTO_LOGNAME):
    import os
    import re
    import pandas as pd

    file_info_list = []
    master_df = pd.DataFrame()

    for filename in os.listdir(input_folder):
        if not filename.endswith(".txt"):
            continue

        # Remove the extension
        base = os.path.splitext(filename)[0]  # strip '.txt'

        try:
            ancestry, subtype, rest = base.split(".", 2)
            parts = rest.split("_")
            # Find the index where ensg_id starts
            ensg_start_idx = next(i for i, part in enumerate(parts) if part.startswith("ENSG"))
            celltype = "_".join(parts[:ensg_start_idx])
            ensg_id = "_".join(parts[ensg_start_idx:])

            file_info_list.append({
                "filename": filename,
                "ancestry": ancestry,
                "subtype": subtype,
                "celltype": celltype,
                "ensg_id": ensg_id
            })

            # Read and filter the file
            file_path = os.path.join(input_folder, filename)
            df = pd.read_csv(file_path)
            filtered_df = df[df["gene"] == ensg_id].copy()

            # Add metadata columns
            filtered_df["ancestry"] = ancestry
            filtered_df["subtype"] = subtype
            filtered_df["celltype"] = celltype

            # Rename columns
            filtered_df = filtered_df.rename(columns={
                "gene": "ensg_id",
                "zscore": "cond_zscore",
                "pvalue": "cond_pvalue"
            })

            master_df = pd.concat([master_df, filtered_df], ignore_index=True)

        except Exception as e:
            print(f"Failed to process file: {filename}. Error: {e}")

    # Select only relevant columns
    if not master_df.empty:
        master_df = master_df[["ensg_id", "ancestry", "subtype", "celltype", "cond_zscore", "cond_pvalue"]]

        # Load tracking file and join
        try:
            tracker_df = pd.read_csv(tracker_path, sep="\t")
            merged_df = tracker_df.merge(master_df,
                                         on=["ensg_id", "ancestry", "subtype", "celltype"],
                                         how="left")
            output_path = outputs[0].filepath
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            merged_df.to_csv(output_path, sep="	", index=False)
            print(f"Merged results written to {output_path}")
        except Exception as e:
            print(f"Failed to load or merge tracking file: {e}")

    return f"Parsed and combined {len(file_info_list)} files from {input_folder}"


@python_app(cache=True, executors=["imputer"])
def gcta_cojo(cojo_ma_file, snp_list_path, cond_snp_list_path, bfile, odir,
              out_prefix, single_snp=None, inputs=[], outputs=[], stdout=parsl.AUTO_LOGNAME,
              stderr=parsl.AUTO_LOGNAME):
    import os
    import subprocess
    import pandas as pd
    cojo_output, parse_ready_output = outputs
    os.makedirs(odir, exist_ok=True)

    extract_path = f"{out_prefix}.cojo.snplist"
    single_cond_path = f"{out_prefix}.single.cond" if single_snp else cond_snp_list_path

    # Write the single SNP condition file if applicable
    if single_snp:
        with open(single_cond_path, 'w') as f:
            f.write(f"{single_snp}\n")

    pre_cojo = f"""
    IN_DIR=`dirname {extract_path}`
    mkdir -p $IN_DIR
    cat {snp_list_path} {single_cond_path} | sort --unique > {extract_path}
    cat {snp_list_path} {single_cond_path} | sort | uniq -c | tr -s ' ' | grep -Pv '^ 1'
    OUT_DIR=`dirname {out_prefix}`
    mkdir -p $OUT_DIR
    """

    bash_command = f"""
    {pre_cojo}

    gcta-1.94.1 \
    --bfile {bfile} \
    --cojo-file {cojo_ma_file} \
    --cojo-cond {single_cond_path} \
    --extract {extract_path} \
    --out {out_prefix}
    """
    subprocess.run(bash_command, shell=True, executable="/bin/bash").check_returncode()

    if os.path.exists(cojo_output):
        cojo_df = pd.read_csv(cojo_output, delim_whitespace=True)
        snp_list_df = pd.read_csv(snp_list_path, header=None, names=["SNP"], delim_whitespace=True)

        missing_snps = set(snp_list_df['SNP']) - set(cojo_df['SNP'])
        if missing_snps:
            additional_rows = pd.DataFrame({
                'SNP': list(missing_snps),
                'bC': 0,
                'bC_se': 1,
                'pC': 1
            })
            cojo_df = pd.concat([cojo_df, additional_rows], ignore_index=True)

        cojo_df.loc[cojo_df['bC'].isna(), ['bC', 'bC_se', 'pC']] = [0, 1, 1]
        cojo_df.to_csv(parse_ready_output, sep='\t', index=False)
    else:
        raise FileNotFoundError(f"Expected output file not found: {cojo_output}")

    return


