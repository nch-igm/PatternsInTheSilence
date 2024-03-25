#!/usr/bin/env python

""" 

Reads in thresholds file and index to select the set of thresholds, grabs variants
from specified cohort and context set and map to genes (background gene list)
Variants meeting specified threshold map to the input gene list. Run ORA against 
two databases. Optionally save version of scored cohort file with the variants used for
the analysis flagged.

Usage: python gene_ora_on_contexts.py \
    -cf <scored cohort file> -thf <thresholds file> -thi <thresholds index> \
        -tag <cohort tag> -cs <context set> -g <group column> \
            -s <score column> --save-cohort-table \
                -o <output folder> -t <output tag>
Reads: scored cohort file, thresholds file
Writes: <output folder>gseapy_ora_kegg_go_cohort<cohort tag>_<output tag>.tsv
<output folder>input_gene_list_cohort<cohort tag>_<output tag>.tsv
<output folder>bg_gene_list_cohort<cohort tag>_<output tag>.tsv
<output folder>scored_variants_flagged_filtered_cohort<cohort tag>_<output tag>.tsv

"""

import argparse
import pandas as pd
import gseapy as gp

def get_variant_set_by_filtering(th_df,
                                 df,
                                 th_index,
                                 comparison,
                                 score_column,
                                 group_column) :
    """
    from thresholds table, select thresholds and filter df according to threshold x 
    group on score_column
    """
    
    #Get thresholds and build map
    th_filtered_df = th_df.query("comparison == @comparison").\
        query("threshold_index == @th_index")
    thresholds_map = th_filtered_df.set_index("group")["threshold_value"].to_dict()
    print("Using thresholds on context scores:",
          thresholds_map)
    #-get relevant threshold for each row
    thresholds_column = df[group_column].map(thresholds_map)
        
    if comparison == "above" :
        flag_against_threshold = \
            (df[score_column] > thresholds_column)
    else :
        flag_against_threshold = \
            (df[score_column] <= thresholds_column)
                
    #Grab variants that pass
    #-flag variants that pass (modifies external df)
    df["gene_set_flag"] = 0
    df.loc[flag_against_threshold,"gene_set_flag"] = 1
    #-filter down to this set
    variant_set_df = df[flag_against_threshold].copy()
    
    return variant_set_df

def get_filtered_cohort_table (df, 
                               group_col,
                               context_set,
                               cohort_tags) :
    """
    filter df entries to context_set in group_col, and for records
    marked in the cohort set
    """
    
    df_on_context = df[df[group_col].isin(context_set)]
    select_against_any = (df_on_context[cohort_tags].sum(axis=1) > 0)
    df_filtered = df_on_context[select_against_any].copy()
    
    return df_filtered

def get_save_bg_list (df,
                      bg_filename) :
    """
    copy gene-related info from df and save data frame, and
    return gene symbol set
    """
    
    bg_gene_df = df[["entrezgene","name","symbol"]].\
            drop_duplicates(ignore_index=True)
            
    print("Saving background genes list to:", bg_filename)
    bg_gene_df.to_csv(bg_filename,
                      sep="\t",
                      index=False)
    #get gene list for now:
    bg_gene_set = list(df["symbol"].unique())
    
    return bg_gene_set

def save_input_gene_list (var_df, gene_filename) :
    """
    copy gene-related info from df and save data frame
    """
    
    gene_df = var_df[["entrezgene","name","symbol"]].\
            drop_duplicates(ignore_index=True)
    print("Saving input gene list to:", gene_filename)
    gene_df.to_csv(gene_filename,
                   sep="\t",
                   index=False)
    
    return 0

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--scored-cohort-file', '-cf', 
                        help="Path to TSV containing sSNV and Sequence Context information with specified cohorts tagged",
                        required=True)
    parser.add_argument('--thresholds-file', '-thf',
                        help="File containing score thresholds for various bins")
    parser.add_argument('--background-genes-file', '-bg',
                        help="Supply or write background gene list file to use with ORA")
    parser.add_argument('--context-set', '-cs', nargs="+",
                        help="Set of contexts to parse score enrichment on")
    parser.add_argument('--tcga-tags', '-tag', type=str, nargs="+",
                        help="Cohort tags to filter for. Note:All tagged variants are processed together.")
    parser.add_argument('--threshold-index', '-thi', type=int,
                        help="Entry in threshold_index column to set for each context in the set")
    parser.add_argument('--comparison', choices=['above', 'below'],
                        default="below",
                        help="Set whether to isolate scores 'above' or 'below' threshold")
    parser.add_argument('--group-column', '-g', type=str,
                        default="SNVContext",
                        help="Column in scored_cohort_file to group variants by")
    parser.add_argument('--score-column', '-s', 
                        help="Column name for score to be referenced")
    parser.add_argument('--save-cohort-table', action='store_true',
                        help="Set if you want to save cohort table with gene set variants marked")
    parser.add_argument('--output-folder', '-o', 
                        help="Path/ to label output files with",
                        default="../../data/4_tcga_analysis/")
    parser.add_argument('--output-tag', '-t', 
                        help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    #Read in scored variant table
    cc_score_cohort_df = pd.read_csv(args.scored_cohort_file,
                                     sep="\t",
                                     dtype={"CHR":str})
    
    #Set databases to query
    db_list = ["KEGG_2021_Human", "GO_Biological_Process_2023"]
    
    #Get threshold set
    if args.thresholds_file is not None :
        #Read in thresholds table and filter for desired values
        thresholds_df = pd.read_csv(args.thresholds_file,
                                sep="\t")
        
        #If no context set given, take whatever is in table
        if args.context_set is None :
            context_set = list(thresholds_df["group"].unique())
        else :
            context_set = args.context_set
            #filter to thresholds needed
            thresholds_df = thresholds_df.query("group.isin(@context_set)").copy()
    else :
        #Just use all contexts present
        context_set = list(cc_score_cohort_df[args.group_column].unique())
    
    #-
    print("Looking at variants with", args.group_column, 
          "in",context_set)
    
    #-Filter table by:
    #--context set
    #--cohort set
    cc_score_cohort_filtered_df = get_filtered_cohort_table(cc_score_cohort_df,
                                                            args.group_column,
                                                            context_set,
                                                            args.tcga_tags)
    
    if args.threshold_index is not None :
        variant_set_df = get_variant_set_by_filtering(thresholds_df,
                                                      cc_score_cohort_filtered_df,
                                                      args.threshold_index, 
                                                      args.comparison,
                                                      args.score_column,
                                                      args.group_column)
    else :
        #If threshold_index is not set, assume all variants are to be used
        variant_set_df = cc_score_cohort_filtered_df.copy()
    #Get list of genes
    input_gene_list = list(variant_set_df["symbol"].unique())
    
    #Save input gene list
    gene_filename = (args.output_folder+"input_gene_list"+
                     "_cohort"+"_".join(args.tcga_tags)+
                     "_"+args.output_tag+
                     ".tsv")
    save_input_gene_list(variant_set_df, gene_filename)
    
    #Set background gene list
    if args.background_genes_file is not None :
        bg_gene_df = pd.read_csv(args.background_genes_file,
                                 sep="\t")
        bg_gene_set = bg_gene_df["symbol"].to_list()
    else :
        bg_gene_filename = (args.output_folder+"bg_gene_list"+
                            "_cohort"+"_".join(args.tcga_tags)+
                            "_"+args.output_tag+
                            ".tsv")
        bg_gene_set = get_save_bg_list(cc_score_cohort_filtered_df,
                                       bg_gene_filename)

    #Get ORA results
    print("Running ORA with", len(input_gene_list), "genes",
          "against a background list of", len(bg_gene_set), "genes")
    enr_bg = gp.enrichr(gene_list = input_gene_list,
                        gene_sets = db_list,
                        background = bg_gene_set,
                        outdir=None)
    print(enr_bg.results.head())
    
    #Save output
    ora_filename = (args.output_folder+
                    "gseapy_ora_kegg_go"+
                    "_cohort"+"_".join(args.tcga_tags)+
                    "_"+args.output_tag+
                    ".tsv")
    print("Writing ORA results to:", ora_filename)
    enr_bg.results.to_csv(ora_filename,
                          sep="\t",
                          index=False)
    
    if (args.save_cohort_table) :
        cohort_table_filename = (args.output_folder+
                             "scored_variants_flagged_filtered"+
                             "_cohort"+"_".join(args.tcga_tags)+
                             "_"+args.output_tag+
                             ".tsv")
        print("Writing flagged cohort variant table to:",
              cohort_table_filename)
        cc_score_cohort_filtered_df.to_csv(cohort_table_filename,
                                           sep="\t",
                                           index=False)
    else :
        pass
        
    
    
        
    
        
        