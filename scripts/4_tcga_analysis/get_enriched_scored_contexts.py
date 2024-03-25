#!/usr/bin/env python

""" 

For specified threshold (assessed per group in grouping column), measure enrichment
of scores meeting that threshold among variants in scored cohort file, for each individual cohort 
and for the combined cohort. Uses variants in table not marked in cohorts to model a null
distribution

Usage: python get_enriched_scored_contexts.py \
    -cf <scored cohort file> -thf <thresholds file> \
        -cs <context set> -thi <threshold index> -g <grouping column> \
            -tag <TCGA tag> [<TCGA tag>] -s <score column> \
                -o <output folder> -t <output tag>
Reads: scored cohort file, thresholds file
Writes: <output folder>scores_enriched_on_<threshold index>quantile_cohort_<tag>[_<tag>]_<output tag>.tsv

"""

import argparse
import pandas as pd
import numpy as np
import scipy.stats

def score_enrichment_at_threshold_bg_xContext (cc_df,
                                               cc_all_df,
                                               th_df,
                                               cohort_tag,
                                               score_column,
                                               group_column) :
    """
    Measure enrichment of scores in cc_df against threshold, using rate in cc_all_df
    to generated expected counts
    """
    
    results_list = []
    #Loop through contexts in thresholds table
    cc_xContext = cc_df.groupby(group_column)
    cc_all_xContext = cc_all_df.groupby(group_column)
    for i, row in th_df.iterrows() :
        context = row["group"]
        #-filter data frames to context
        cc_con_df = cc_df.loc[cc_xContext.groups[context]]
        cc_all_con_df = cc_all_df.loc[cc_all_xContext.groups[context]]
        #-get lengths of data frames
        total_cohort_context = cc_con_df.shape[0]
        total_context = cc_all_con_df.shape[0]
                
        if row["comparison"] == "above" :
            cohort_on_context = (cc_con_df[score_column] > 
                                 row["threshold_value"])
            all_on_context = (cc_all_con_df[score_column] > 
                              row["threshold_value"])
        else :
            cohort_on_context = (cc_con_df[score_column] <= 
                                 row["threshold_value"])
            all_on_context = (cc_all_con_df[score_column] <= 
                              row["threshold_value"])
        cohort_num_in = cohort_on_context.sum()
        cohort_num_out = total_cohort_context-cohort_num_in

        ratio_in = all_on_context.sum()/total_context
        exp_num_in = round(total_cohort_context*ratio_in)
        exp_num_out = total_cohort_context-exp_num_in
        
        res = scipy.stats.chisquare([cohort_num_in, cohort_num_out],
                                    [exp_num_in, exp_num_out])
        stat = res.statistic
        pv = res.pvalue
        
        result_dict = {"cohort":cohort_tag,
                       "group":group_column,
                       "context":context,
                       "threshold_quantile":row["threshold_quantile"],
                       "threshold_value":row["threshold_value"],
                       "comparison":row["comparison"],
                       "cohort_number_in":cohort_num_in,
                       "cohort_number_out":cohort_num_out,
                       "cohort_expected_in":exp_num_in,
                       "cohort_expected_out":exp_num_out,
                       "test_stat": stat,
                       "test_pv": pv}
        results_list.append(result_dict)
    return results_list
        
if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--scored-cohort-file', '-cf', 
                        help="Path to TSV containing sSNV and Sequence Context information with cohort variants tagged",
                        required=True)
    parser.add_argument('--thresholds-file', '-thf',
                        help="File containing score thresholds for various bins, groups corresponding to group column")
    parser.add_argument('--context-set', '-cs', nargs="+",
                        help="Set of contexts to parse score enrichment on")
    parser.add_argument('--threshold-index', '-thi', type=int,
                        help="Entry in threshold_index column to set for each context in the set")
    parser.add_argument('--comparison', choices=['above', 'below'],
                        default="below",
                        help="Set whether to isolate scores 'above' or 'below' threshold")
    parser.add_argument('--group-column', '-g', type=str,
                        default="SNVContext",
                        help="Column in scored_cohort_file to group variants by")
    parser.add_argument('--tcga-tags', '-tag', type=str, 
                        nargs="+",
                        help='Column names that tag TCGA cohort membership')
    parser.add_argument('--score-column', '-s', help="Column name for score to be referenced")
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/4_tcga_analysis/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    #Read in annotated variant table
    usecols = ["CHR", "POS", "REF", "ALT", "NM_ID",
               args.group_column, args.score_column,
               "y"] + list(args.tcga_tags)
    cc_score_cohort_df = pd.read_csv(args.scored_cohort_file,
                                     sep="\t",
                                     dtype={"CHR":str},
                                     usecols=usecols)
    #-if we didn't fill in the cohort columns, do so now
    cc_score_cohort_df.fillna(value={c:False for c in args.tcga_tags},
                              inplace=True)
    
    #Read in thresholds table and filter for desired values
    thresholds_df = pd.read_csv(args.thresholds_file,
                                sep="\t")
    #If no context set given, take whatever is in table
    if args.context_set is None :
        context_set = list(thresholds_df["group"].unique())
    else :
        context_set = args.context_set
    #-
    print("Looking at variants with", args.group_column, 
          "in",context_set)

    #Filter for selected contexts
    thresholds_filtered_df = thresholds_df.query("comparison == @args.comparison").\
        query("group.isin(@context_set)").\
            query("threshold_index == @args.threshold_index")
    #-
    cc_score_cohort_filtered_df = cc_score_cohort_df[cc_score_cohort_df[args.group_column].isin(context_set)]
            
    #Loop through cohorts
    results_list = []
    for cohort in args.tcga_tags :
        cc_at_cohort_df = cc_score_cohort_filtered_df[cc_score_cohort_filtered_df[cohort]]
        print("For cohort", cohort+": ", cc_at_cohort_df.shape[0], "variants")
                
        cohort_results_list = score_enrichment_at_threshold_bg_xContext(cc_at_cohort_df,
                                                                     cc_score_cohort_filtered_df,
                                                                     thresholds_filtered_df,
                                                                     cohort,
                                                                     args.score_column,
                                                                     args.group_column)
        results_list.extend(cohort_results_list)
    #Run with entire table if needed (combined cohort):
    if len(args.tcga_tags) > 1 :
        cohort="grouped_TCGA"
        select_against_any = (cc_score_cohort_filtered_df[args.tcga_tags].sum(axis=1) > 0)
        cc_at_cohort_df = cc_score_cohort_filtered_df[select_against_any]
        cohort_results_list = score_enrichment_at_threshold_bg_xContext(cc_at_cohort_df,
                                                                     cc_score_cohort_filtered_df,
                                                                     thresholds_filtered_df,
                                                                     cohort,
                                                                     args.score_column,
                                                                     args.group_column)
        results_list.extend(cohort_results_list)
    else :
        pass
        
    #Convert to dataframe
    results_df = pd.DataFrame(results_list)
    print("Summary:")
    print(results_df.head())
    
    #Save output
    out_filename = (args.output_folder+
                    "scores_enriched_on_"+str(args.threshold_index)+"quantile"+
                    "_cohort"+"_".join(args.tcga_tags)+
                    "_"+args.output_tag+
                    ".tsv")
    print("Writing to:", out_filename)
    results_df.to_csv(out_filename,
                      sep="\t",
                      index=False)
    