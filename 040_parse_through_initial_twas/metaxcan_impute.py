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

@bash_app(cache=True, executors=["imputer"])
def metaxcan_impute_gwas_part(gwas_file, by_region_file, pq_geno_file,
        pq_geno_meta_file, chromosome, cur_sub_batch,
        window=100000, parsimony=7, regularization=0.1,
        freq_filter=0.01, total_sub_batches=10, 
        inputs=[], outputs=[],
                        stdout=parsl.AUTO_LOGNAME, 
                        stderr=parsl.AUTO_LOGNAME): 
    import os
    # if os.path.exists(outputs[0].filepath): 
    #     return("echo 'Output exists. Remove it or delete it.'")
    output = outputs[0].filepath



    # `bash_command` below does not need -chromosome argument
    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_summary_imputation.py \
        -by_region_file {by_region_file} \
        -gwas_file {gwas_file} \
        -parquet_genotype {pq_geno_file} \
        -parquet_genotype_metadata {pq_geno_meta_file} \
        -window {window} \
        -parsimony {parsimony} \
        -regularization {regularization} \
        -frequency_filter {freq_filter} \
        -sub_batches {total_sub_batches} \
        -sub_batch {cur_sub_batch} \
        --standardise_dosages \
        -output {output}
    """
    return(bash_command) 

@bash_app(cache=True, executors=["parser"])
def metaxcan_postprocess_gwas(gwas_file, parts_folder, pattern, parsimony=9, 
                        inputs=[], outputs=[],
                        stdout=parsl.AUTO_LOGNAME, 
                        stderr=parsl.AUTO_LOGNAME): 
    import os
    # if os.path.exists(outputs[0].filepath):
    #     return("echo 'Output exists. Remove it or delete it.'")
    output = outputs[0].filepath
    # gwas_file.result()

    bash_command = \
    f"""
    python3 /gpfs/data/gao-lab/Julian/software/summary-gwas-imputation/src/gwas_summary_imputation_postprocess.py \
        -gwas_file {gwas_file} \
        -folder {parts_folder} \
        -pattern {pattern} \
        -parsimony {parsimony} \
        -output {output}
    """
    return(bash_command) 


# @join_app(cache=True)
def metaxcan_impute_gwas(gwas_file, ld_block_format, pq_geno_file_format, pq_geno_meta_fformat,
        postproc_pattern,
        output_format,
        chromosome=None,
        window=100000, parsimony=7, regularization=0.1,
        freq_filter=0.01, total_sub_batches=10, 
        output=None):

    import os
    region_folder, oput_fpath = os.path.split(output_format)
    region_folder = region_folder + "/" # prlly not necessary

    # Create output directory if it does not exist
    os.makedirs(region_folder, exist_ok=True)

    # Create output directory if it does not exist
    output_dir = os.path.dirname(output[0])
    os.makedirs(output_dir, exist_ok=True)

    impute_part_res = []
    if chromosome is not None:
        chr_range = [chromosome]
    else:
        chr_range = range(1,23)

    for chr_num in chr_range:
        for batch_num in range(total_sub_batches):
            ofile_name = output_format.format(chr_num=chr_num, sub_batch=batch_num,
                    reg=regularization, ff=freq_filter)
            impute_part_res.append(
                    metaxcan_impute_gwas_part(gwas_file,
                        ld_block_format.format(chr_num=chr_num),
                        File(pq_geno_file_format.format(chr_num=chr_num)),
                        File(pq_geno_meta_fformat.format(chr_num=chr_num)), # pq_geno_meta_file,
                        chromosome=chr_num,
                        cur_sub_batch=batch_num,
                        outputs=[File(ofile_name)]
                ))

    # Wait for imputation individual parts to finish
    # finished_impute_parts = [res.result() for res in impute_part_res]
    
    return(metaxcan_postprocess_gwas(gwas_file, region_folder, pattern=postproc_pattern, inputs=[res.outputs[0] for res in impute_part_res], outputs=output))

# def metaxcan_impute_gwas_part(gwas_file, by_region_file, pq_geno_file,
#         pq_geno_meta_file, chromosome, cur_sub_batch,
#         window=100000, parsimony=7, regularization=0.1,
#         freq_filter=0.01, total_sub_batches=10, 
    return(bash_command) 
