#!/usr/bin/env python

""" 

Load scored variant table that has gene_set_flag column to be used as source for random variant 
sampling from unflagged part of table, and to tally flagged variants by context set.
Use variant sample to map to gene lists and run ORA on that gene list with provided background gene
list (repeat for specified number of random sets)

Usage: python gene_ora_on_random_gene_lists_byContext.py \
    -bvf <bg variant file> -bg <bg genes file> -g <group column> \
        -nr <# random sets> -tag <cohort tags> \
            -o <output folder> -t <output tag>
Reads: bg variant file, bg genes file, db filenames
Writes: <output folder>gseapy_sampled_ora_kegg_go_cohort<cohort tag>_n<# random sets>_<output tag>.tsv
<output folder>gene_sampled_matchedContext_cohort<cohort tag>_n<# random sets>_<output tag>.tsv 

"""

import argparse
import pandas as pd
import time
import gseapy as gp

random_states = range(1000)

def get_sampled_gene_set (base_df,
                          group_column,
                          target_counts,
                          random_state_i) :
    """
    For each group in target_counts, take that sup-group (based on group column) and
    sample to the target size set. Concat sub-samplings across sub-groups
    """
    
    sampled_dfs = []
    base_xGroup = base_df.groupby(group_column)
    
    for context, base_context_df in base_xGroup :
        if context in target_counts :
            sample_base_context_df = base_context_df.sample(n=target_counts[context],
                                                            replace=False,
                                                            random_state=random_state_i)
            sampled_dfs.append(sample_base_context_df)
        else :
            pass
    sampled_df = pd.concat(sampled_dfs,
                           axis=0)
    print(sampled_df[group_column].value_counts())
    return sampled_df

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--background-variant-file', '-bvf', 
                        help="Scored variant file with a gene_set_flag column",
                        required=True)
    parser.add_argument('--background-genes-file', '-bg',
                        help="Path to background gene list to be used with ORA")
    parser.add_argument('--group-column', '-g', type=str,
                        default="SNVContext",
                        help="Column in scored_cohort_file to group variants by")
    parser.add_argument('--num-random-sets', '-nr', type=int,
                        help="Number of random gene sets to generate")
    parser.add_argument('--tcga-tags', '-tag', type=str, nargs="+",
                        help="Cohort tags to filter for. Note:All tagged variants are processed together.")
    parser.add_argument('--db-filenames', '-db', nargs="+",
                        default=["../../data/4_tcga_analysis/enrichr_gene_sets/KEGG_2021_Human.tsv",
                                 "../../data/4_tcga_analysis/enrichr_gene_sets/GO_Biological_Process_2023.tsv"],
                        help="Enrichr database files that can be read into dictionaries")
    parser.add_argument('--save-gene-array', action='store_true',
                        help="Set to store version of background gene list marked by set inclusions.")
    parser.add_argument('--max-pv', type=float, default=0.05,
                        help="Max. adjusted p-value to filter ORA results")
    parser.add_argument('--output-folder', '-o', 
                        help="Path/ to label output files with",
                        default="../../data/4_tcga_analysis/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    #Read in background variant table
    #-assume column named "gene_set_flag"
    bg_var_df = pd.read_csv(args.background_variant_file,
                            sep="\t",
                            dtype={"CHR":str})
    num_genes = len(bg_var_df["symbol"].unique())
    
    #Set databases to query
    db_dict_list = []
    print("Reading databases from:", args.db_filenames)
    for db in args.db_filenames :
        db_df = pd.read_csv(db,
                            sep="\t").\
                                set_index("set_name")
        db_df["gene_list"] = db_df["gene_list"].apply(lambda x: x.split(";"))
        db_dict = db_df.squeeze().to_dict()
        db_dict_list.append(db_dict)
    #db_list = ["KEGG_2021_Human", "GO_Biological_Process_2023"]
    
    #Read in background genes file (should match the set of genes in bg_var_df)
    bg_gene_df = pd.read_csv(args.background_genes_file,
                             sep="\t")
    bg_gene_set = list(bg_gene_df["symbol"])
    if bg_gene_df.shape[0] != num_genes :
        print("Are these compatible lists?", 
              bg_gene_df.shape[0], "vs.", num_genes, "-",
              args.background_genes_file, args.background_variant_file)
    else :
        pass
    
    #-set target distribution of contexts
    target_context_set = bg_var_df.query("gene_set_flag == 1")[args.group_column].value_counts()
    print("Matching this context set:", print(target_context_set))
    print("Filtering ORA results <",args.max_pv)
    #-set table that we can draw from:
    base_variant_df = bg_var_df.query("gene_set_flag == 0")
    
    #Track which genes were sampled in which set
    sampled_gene_cols = []
    #-track significant ORA results
    ora_result_list = []
    for i in range(args.num_random_sets) :
        #Get sampled gene set
        sampled_df = get_sampled_gene_set (base_variant_df,
                                           args.group_column,
                                           target_context_set,
                                           random_states[i])
        sampled_gene_set = list(sampled_df["symbol"].unique())
        sampled_gene_set_i = 1*(bg_gene_df["symbol"].isin(sampled_gene_set))
        sampled_gene_cols.append(pd.Series(sampled_gene_set_i,
                                           name="set_"+str(i)).to_frame())
        
        #Run ORA on gene set
        #Get ORA results
        print("Set", i, ": running ORA with", len(sampled_gene_set), "genes",
              "against a background list of", len(bg_gene_set), "genes")
        enr_rand_bg = gp.enrichr(gene_list = sampled_gene_set,
                                 gene_sets = db_dict_list, #db_list,
                                 background = bg_gene_set,
                                 outdir=None)
        enr_rand_bg_df = enr_rand_bg.results
        enr_rand_bg_df["set_number"] = i
        enr_rand_bg_filtered = enr_rand_bg_df[enr_rand_bg_df["Adjusted P-value"] < args.max_pv]
        ora_result_list.append(enr_rand_bg_filtered)
        print(" smallest adjusted p-value:", enr_rand_bg_df["Adjusted P-value"].min())
        
    #Concat results
    sampled_ora_result_df = pd.concat(ora_result_list,
                                      axis=0,
                                      ignore_index=True)
    
    #-write results
    out_filename = (args.output_folder+
                    "gseapy_sampled_ora_kegg_go"+
                    "_cohort"+"_".join(args.tcga_tags)+
                    "_n"+str(args.num_random_sets)+
                    "_"+args.output_tag+
                    ".tsv")
    
    print("Writing", sampled_ora_result_df.shape[0], "results to",
          out_filename)
    sampled_ora_result_df.to_csv(out_filename,
                                 sep="\t",
                                 index=False)
    
    if args.save_gene_array :
        out_gene_filename = (args.output_folder+
                             "gene_sampled_matchedContext"+
                             "_cohort"+"_".join(args.tcga_tags)+
                             "_n"+str(args.num_random_sets)+
                             "_"+args.output_tag+
                             ".tsv")
        sampled_gene_df = pd.concat([bg_gene_df]+
                                    sampled_gene_cols,
                                    axis=1)
        sampled_gene_df.to_csv(out_gene_filename,
                               sep="\t",
                               index=False)
        print("Writing gene set members to:",out_gene_filename)
    else :
        pass
        
        
    
        
    
    