from parsl.data_provider.files import File

#--keep_all_original_entries \ # Don't need this for the COJO parse or the standard parse thankfully
class MetaxcanParse:
    def __init__(self):
        self.parsing_info = {
                "james_gwas_brca_stype_by_ances": """\
                    -output_column_map ID variant_id \
                    -output_column_map BaselineAllele non_effect_allele \
                    -output_column_map EffectAllele effect_allele \
                    -output_column_map BETA effect_size \
                    -output_column_map SE standard_error \
                    -output_column_map EAF frequency \
                    -output_column_map N sample_size \
                    -output_column_map N_case n_cases \
                    -split_column variant_id ':' chromosome position nea ea \
                    --chromosome_format \
                     -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
                """,
                "james_gwas_brca_stype_by_ances_cojo": """\
                    -output_column_map SNP variant_id \
                    -output_column_map bC effect_size \
                    -output_column_map bC_se standard_error \
                    -output_column_map pC pvalue \
                    -output_column_map freq frequency \
                    -output_column_map n sample_size \
                    -output_column_map n n_cases \
                    -split_column variant_id ':' chromosome position non_effect_allele effect_allele \
                    --chromosome_format \
                     -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
                     """
                }

# Misc. Helper functions
def filter_res(chr_num, app_fut):
    """
    Function to filter a list of DataFuture Objects. Combine with partial from functools.
    """
    import re
    has_match = re.search("chr{}_".format(chr_num), app_fut.outputs[0].filepath)
    return has_match

def filter_fp(chr_num, file_obj):
    """
    Function to filter a list of DataFuture Objects. Combine with partial from functools.
    """
    import re
    has_match = re.search("chr{}_".format(chr_num), file_obj.filepath)
    return has_match

def make_tmp_list_file(files):
    from uuid import uuid4
    import pandas as pd
    try:
        df = pd.DataFrame(data=[i.filepath for i in files])
    except AttributeError: 
        df = pd.DataFrame(data=[file_str for file_str in files])
    tmp_file = "tmp/" + str(uuid4()) + ".tmp"
    df.to_csv(tmp_file, index=False, header=False)
    return(tmp_file)

def extract_fp_strings(files):
    r_strings = []
    for fobj in files:
        r_strings.append(str(fobj.filepath))
    return(r_strings)
