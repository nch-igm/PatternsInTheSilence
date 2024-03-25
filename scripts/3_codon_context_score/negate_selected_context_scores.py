#!/usr/bin/env python

""" 

Negate contents of specified column for specified subsets of group column

Usage: python negate_selected_context_scores.py \
    -vf <scored variant file> -s <score column> \
        -g <grouping column> -r <rescale groups> \
            -o <output folder> -f <output filename>
Reads: scored variant file
Writes: <output folder><output filename>

"""

import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--scored-variant-file', '-vf', 
                        help="Path to TSV containing sSNV and Sequence Context information",
                        required=True)
parser.add_argument('--score-column', '-s', 
                        help="Column name for score to be rescaled")
parser.add_argument('--group-column', '-g', type=str,
                        help="Column to group scores by when rescaling scores")
parser.add_argument('--rescale-groups', '-r', nargs="+",
                        help="Values of 'group column' to rescale context scores within.")
parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/3_codon_context_score/")
parser.add_argument('--output-filename', '-f', help="Name to write file to in folder.",
                        type=str,
                        default="")
args = parser.parse_args()
    
#Read in scored variant file
cc_score_df = pd.read_csv(args.scored_variant_file,
                          sep="\t",
                          dtype={"CHR":str})
cc_score_df.head()
    
#Negate scores of selected groups
print("Groups rescaled:", args.rescale_groups)
cc_score_df.loc[cc_score_df[args.group_column].isin(args.rescale_groups),
                args.score_column] *= -1

#Save the file
out_filename = args.output_folder+args.output_filename

cc_score_df.to_csv(out_filename,
                   sep="\t",
                   index=False)
