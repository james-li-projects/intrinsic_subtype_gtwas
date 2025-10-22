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
from options import MetaxcanParse

METAXCAN_PARSE_OPTS = MetaxcanParse()

@bash_app(cache=True, executors=["imputer", "parser"])
def spredixcan(gwas_file, db_file, covar_file,
        model_db_snp_key="varID",
        inputs=[], outputs=[],
                        stdout=parsl.AUTO_LOGNAME, 
                        stderr=parsl.AUTO_LOGNAME): 
    import os
    # if os.path.exists(outputs[0].filepath): 
    #     return("echo 'Output exists. Remove it or delete it.'")

    # Create output directory if it does not exist
    output_filepath = outputs[0].filepath
    output_dir = os.path.dirname(output_filepath)
    os.makedirs(output_dir, exist_ok=True)

    output = outputs[0].filepath

    # `bash_command` below does not need -chromosome argument
    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/MetaXcan/software/SPrediXcan.py \
        --gwas_file {gwas_file} \
        --snp_column panel_variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --zscore_column zscore \
        --model_db_path {db_file} \
        --covariance {covar_file} \
        --keep_non_rsid --additional_output \
        --model_db_snp_key {model_db_snp_key} \
        --throw \
        --output_file {output}
    """
    return(bash_command) 
