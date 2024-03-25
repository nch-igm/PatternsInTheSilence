#!/usr/bin/env python

""" 

Usage: python score_filter_tcga_variants.py \
    -vf <scored variant file> -tf <tcga file> [-tf ...] \
        -tag <tcga tag> [-tag ...] -ac <AC max> -m <merge type> \
            -o <output folder> -t <output tag>
Reads: scored variant file, tcga file(s)
Writes: <output folder>context_scored_variant_frq_cohorts<tag>[_<tag>]_AC<AC max>_<merge type>_<output tag>.tsv

"""

import argparse
import pandas as pd

def read_in_merge_freq (filename, cohort_tag) :
    """
    read in variant frequency file of cohort
    """
    
    freq_df = pd.read_csv(filename,
                          sep="\t",
                           dtype={"CHR":str}).\
    drop(columns=["A1","A2"]).\
    rename(columns={"MAF":"MAF_"+cohort_tag,
                    "NCHROBS":"NCHROBS_"+cohort_tag})
    
    freq_count_df = freq_df[freq_df["NCHROBS_"+cohort_tag] > 0].copy()
    freq_count_df["ALT_CT_"+cohort_tag] = \
        round(freq_count_df["MAF_"+cohort_tag]*
              freq_count_df["NCHROBS_"+cohort_tag]).astype(int)
    freq_count_df[cohort_tag] = "True"
    
    print("--",cohort_tag+":")
    print(freq_count_df.head())
    
    return freq_count_df

def merge_freq_dfs (freq_df_dict, key_list) :
    """
    Merge freq_dfs across cohorts
    """
    
    num_keys = len(key_list)
    
    merge_df = freq_df_dict[key_list[0]].merge(freq_df_dict[key_list[1]],
                                               on=["CHR", "SNP", "POS", "REF", "ALT"],
                                               how="outer")
    if num_keys > 2 :
        for fi in range(2,num_keys) :
            merge_df = merge_df.merge(freq_df_dict[key_list[fi]],
                                      on=["CHR", "SNP", "POS", "REF", "ALT"],
                                      how="outer")
    else :
        pass
    
    merge_fill_df = merge_df.fillna(value={k:"False" for k in key_list})
    merge_df = merge_fill_df.fillna(0)
    
    return merge_df

def filter_scored_table (scored_df, ac_max) :
    """
    filter variant score table for max allele count and gnomAD covg (y column)
    """
    
    scored_covg_df = scored_df.query("y > -1")
    
    if scored_covg_df[["gnomAD2_EX_AC",
                       "gnomAD3_WG_AC"]].isnull().values.any() :
        print("Do you need to fill in NA values?")
    else :
        pass

    scored_covg_filtered_df = scored_covg_df[(scored_covg_df["gnomAD2_EX_AC"]+
                                              scored_covg_df["gnomAD3_WG_AC"]) 
                                             <= ac_max].\
                                                  copy()
    
    return scored_covg_filtered_df

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--scored-variant-file', '-vf', 
                        help="Path to TSV containing sSNV and Sequence Context information",
                        required=True)
    parser.add_argument('--tcga-files', '-tf', type=str, 
                        nargs="+", action="extend",
                        help="Paths to variant .freq files")
    parser.add_argument('--tcga-tags', '-tag', type=str, 
                        nargs="+", action="extend",
                        help="Tags to assign to variant files passed to --tcga-files. Use same ordering.")
    parser.add_argument('--ac-max', '-ac', type=int,
                        help="Maximum to set for gnomAD allele count when filtering")
    parser.add_argument('--merge-type','-m',choices=["all_syn","intx","all_tcga"],
                        help="Type of merged table to store: all in scoring table (all_syn), intersection (intx), all in TCGA table (all_tcga)")
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/4_tcga_analysis/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    #Set up cohort dict
    if len(args.tcga_tags) != len(args.tcga_files) :
        print("Input lists don't match in size.")
    else :
        cohort_dict = {tc:args.tcga_files[i]
                       for i, tc in enumerate(args.tcga_tags)}
    print("Using these cohort files:",
          cohort_dict)
    
    #Import the overall list of scored synonymous variants
    cc_score_df = pd.read_csv(args.scored_variant_file,
                          sep="\t",
                          dtype={"CHR":str})
    cc_score_df.head()
    
    #Read in the cohort files
    cohort_df_dict = {}
    for cohort in cohort_dict :
        cohort_filename = cohort_dict[cohort]
    
        cohort_df = read_in_merge_freq (cohort_filename, cohort)
        cohort_df_dict[cohort] = cohort_df
    #-merge dataframes
    cohort_merge_list = list(cohort_df_dict.keys())
    if len(cohort_merge_list) > 1 :
        merged_cohort_df = merge_freq_dfs (cohort_df_dict, 
                                           cohort_merge_list)
    else :
        merged_cohort_df = cohort_df
    
    #Match up scored ssnv and TCGA tables
    cc_score_tcga_df = cc_score_df.merge(merged_cohort_df,
                                         on=["CHR", "POS", "REF", "ALT"],
                                         how="outer", 
                                         indicator="_vcf")
    print("Matching sSNV table with TCGA table:")
    print(cc_score_tcga_df["_vcf"].value_counts())
    
    #-select all sSNVs, including those in TCGA set
    if args.merge_type == "all_syn" :
        cc_score_tcga_merged_df = cc_score_tcga_df[(cc_score_tcga_df["_vcf"] == "left_only")|
                                                   (cc_score_tcga_df["_vcf"] == "both")]
    elif args.merge_type == "intx" :
        cc_score_tcga_merged_df = cc_score_tcga_df[(cc_score_tcga_df["_vcf"] == "both")]
    else :
        cc_score_tcga_merged_df = cc_score_tcga_df[(cc_score_tcga_df["_vcf"] == "right_only")|
                                                   (cc_score_tcga_df["_vcf"] == "both")]
    
    cc_score_tcga_merged_gnomad_df = filter_scored_table(cc_score_tcga_merged_df,
                                                         args.ac_max)
    
    del cc_score_df, cc_score_tcga_df
    
    out_filename = (args.output_folder+
                    "context_scored_variant_frq_cohorts"+
                    "_".join(cohort_merge_list)+
                    "_AC"+str(args.ac_max)+
                    "_"+args.merge_type+
                    "_"+args.output_tag+
                    ".tsv")
    
    print("Writing to", out_filename)
    
    cc_score_tcga_merged_gnomad_df.to_csv(out_filename,
                                          sep="\t",
                                          index=False)
    
    