#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import parsl
import subprocess
from subprocess import call
from parsl.app.app import bash_app, python_app, join_app
from parsl.dataflow.futures import AppFuture

from config import id_for_memo_File
from options import MetaxcanParse

METAXCAN_PARSE_OPTS = MetaxcanParse()

@bash_app(cache=True, executors=["parser"])
def metaxcan_parse_gwas(gwas_path, ref_metadata_path, parse_scheme_key, inputs=[], outputs=[],
                        stdout=parsl.AUTO_LOGNAME, 
                        stderr=parsl.AUTO_LOGNAME): 
    import os
    # if os.path.exists(outputs[0].filepath): 
    #     return("echo 'Output exists. Remove it or delete it.'")
    output = outputs[0].filepath

    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_parsing.py \
        -gwas_file {gwas_path} \
        {METAXCAN_PARSE_OPTS.parsing_info[parse_scheme_key]} \
        -snp_reference_metadata {ref_metadata_path} METADATA \
        -output {output}
    """
    return(bash_command) 

@bash_app(cache=True, executors=["express", "imputer"])
def small_metaxcan_parse_gwas(gwas_path, ref_metadata_path, parse_scheme_key, inputs=[], outputs=[],
                        stdout=parsl.AUTO_LOGNAME, 
                        stderr=parsl.AUTO_LOGNAME): 
    import os
    # if os.path.exists(outputs[0].filepath): 
    #     return("echo 'Output exists. Remove it or delete it.'")
    output = outputs[0].filepath

    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_parsing.py \
        -gwas_file {gwas_path} \
        {METAXCAN_PARSE_OPTS.parsing_info[parse_scheme_key]} \
        -snp_reference_metadata {ref_metadata_path} METADATA \
        -output {output}
    """
    return(bash_command) 

@python_app(cache=True)
def make_cojo_ma_file(study, outputs=[],
                      stdout=parsl.AUTO_LOGNAME, 
                      stderr=parsl.AUTO_LOGNAME): 
    import os
    import numpy as np
    import pandas as pd
    if not os.path.exists(os.path.dirname(outputs[0].filepath)):
        os.makedirs(os.path.dirname(outputs[0].filepath))
    if os.path.exists(outputs[0].filepath):
        return("echo 'Output exists. Remove it or delete it.'")
    output = outputs[0].filepath
    pre_output = os.path.join(os.path.dirname(output), "not_final_" + os.path.basename(output))
    bash_command = \
    f"""
    zcat {setup.GWAS_KEY[study]} | awk 'BEGIN {{FS="\t"}}; {{print $3"_"$4" "$5" "$6" "$8" "$12" "$13" "$11" "$10" {setup.SAMPLE_N_KEY[study]}"}}' > {pre_output}
    """
    proc = subprocess.run(bash_command, shell=True)
    if proc.returncode == 0:
        try:
            pre_cojo_ma_names = ["chromosome_position", "effect_allele", "non_effect_allele", "frequency", "effect_size", "standard_error", "pvalue", "zscore", "N"]
            pre_cojo_ma_df = pd.read_csv(pre_output, header=0, names=pre_cojo_ma_names, delim_whitespace=True)

            pre_cojo_ma_df['effect_size_alt'] = pre_cojo_ma_df['zscore'] / np.sqrt(setup.SAMPLE_N_KEY[study])
            pre_cojo_ma_df['N'][pre_cojo_ma_df['N'].isnull()] = setup.SAMPLE_N_KEY[study]

            pre_cojo_ma_df['effect_size'][pre_cojo_ma_df['effect_size'].isnull()] = pre_cojo_ma_df['effect_size_alt'][pre_cojo_ma_df['effect_size'].isnull()]
            pre_cojo_ma_df['standard_error'][pre_cojo_ma_df['standard_error'].isnull()] = 1 / np.sqrt(setup.SAMPLE_N_KEY[study]) 

            cojo_ma_names = ["chromosome_position", "effect_allele", "non_effect_allele", "frequency", "effect_size", "standard_error", "pvalue", "N"]
            pre_cojo_ma_df[cojo_ma_names].to_csv(output, sep=" ", header=True, index=False)
            os.remove(pre_output)
        except FileNotFoundError:
            print("Shouldn't have no file.")
    else:
        raise Exception(f"shell command failed? {proc}")
    return
