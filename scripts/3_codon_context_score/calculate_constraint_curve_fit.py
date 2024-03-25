#!/usr/bin/env python

"""

Fit linear model to constraint-score curves, generated per group in 
grouping column
 
Usage: python calculate_constraint_curve.py \
    -gvf <scored_grouped_variant_file> -s <score column> \
        -y <indicator column> -g <grouping column> \
            -bm <bin size min.> -o <output folder> -t <output tag>
Reads: scored_grouped_variant_file
Writes: <output folder>constraint_curve_fits_<bin size min.>min_<output tag>.tsv

"""

import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm

def get_linear_fit (x, y, w) :
    """
    Set up linear model of data, weighted by number of variants
    in score bin, fit model with WLS
    """
    X_1d = np.array(x)
    X_1d = sm.add_constant(X_1d)
    mod1d_wls = sm.WLS(y, X_1d, weights=np.array(w))
    res_1d_wls = mod1d_wls.fit()
    
    return res_1d_wls

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--scored-grouped-variant-file', '-gvf', help="Path to TSV containing grouped score and maf data",
                        required=True)
    parser.add_argument('--score-column', '-s', help="Column containing binned scores",
                        default="diff_sum_context_score_REF_Codon_zscore_binned_left")
    parser.add_argument('--indicator-column', '-y', help="Column label for indicator variable",
                        default="y")
    parser.add_argument('--group-column', '-g', help="Column used to group scores",
                        default="SNVContext")
    parser.add_argument('--bin-min', '-bm', help="Minimum number of data points to include bin in curve fit",
                        type=int,
                        default=1000)
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/3_codon_context_score/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()
    
    variant_scored_grouped = pd.read_csv(args.scored_grouped_variant_file,
                                         sep="\t")
    bin_minimum = args.bin_min
    variant_scored_grouped_filtered = variant_scored_grouped.query("N >= @bin_minimum")
    
    vsg_xGroup = variant_scored_grouped_filtered.groupby(args.group_column)
    y_column = args.indicator_column+"_mean"
    
    curve_fit_rows = []
    for grp, df_grp in vsg_xGroup:
        x = df_grp[args.score_column]
        y = df_grp[y_column]
        n = df_grp["N"]
        
        r1 = get_linear_fit(x,y,n)
        
        #For WLS, weights are proportional to the inverse of the error variance
        curve_fit_row = {"group_variable":args.group_column,
                         "group_label":grp,
                         "y_variable":y_column,
                         "r1_x0": r1.params.loc["const"],
                         "r1_x1": r1.params.loc["x1"],
                         "r1_x0_pv": r1.pvalues.loc["const"],
                         "r1_x1_pv": r1.pvalues.loc["x1"],
                         "r1_rsq_adj": r1.rsquared_adj}
        curve_fit_rows.append(curve_fit_row)
        
    curve_fit_df = pd.DataFrame(curve_fit_rows)
    output_filename = (args.output_folder + 
                       "constraint_curve_fits_"+
                       str(args.bin_min)+"min_"+
                       args.output_tag+".tsv")
    
    curve_fit_df.to_csv(output_filename,
                        sep="\t",
                        index=False)
