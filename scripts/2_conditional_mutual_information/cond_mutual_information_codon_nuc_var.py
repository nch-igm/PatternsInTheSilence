#!/usr/bin/env python

""" 

Calculate mutual information between codon distribution and a 
metric describing the codon context, then calculate the conditonal mutual 
information between codon and nucleotide distribution after accounting for metric

Usage: python cond_mutual_information_codon_nuc_var.py -f <codon_record_tsv> \
    -l <specify aa list> -a <specify aa column in file> \
        -p <# of positions in context> -b <# of equal distribution bins to use on metric> \
            -v <column to use/transform> -vf <function to apply to column> \
                -mf <metric file> -mw <codon window to calculate metric> \
                    -sl <substring list> -o <output folder> -t <output tag>
Reads: codon_record_tsv, metric file
Writes: <output folder>site_metric<output tag>.tsv
    <output folder>cond_mut_inf_codon_nuc_pos_var<output tag>.tsv
    <output folder>mut_inf_codon_var<output tag>.tsv

"""

import pandas as pd
import numpy as np
import argparse

import importlib.util
import sys

ccv_spec = importlib.util.spec_from_file_location("codon_context_variables", "../codon_context_variables.py")
ccv = importlib.util.module_from_spec(ccv_spec)
ccv_spec.loader.exec_module(ccv)

def get_codon_metric_avg(seq, metric_dict, window_size, third_base_pos=51) :
    """ 
    Get weights of codons in context, take average 
    """
    left_string = seq[:(third_base_pos-3)]
    rght_string = seq[third_base_pos:]

    left_string_select = left_string[-(3*window_size):]
    rght_string_select = rght_string[:(3*window_size)]
    string_select = left_string_select+rght_string_select

    try :
        metric_list = [metric_dict[string_select[x:(x+3)]]
                       for x in range(0,len(string_select),3)]
    except :
        print("Can't convert codon:",string_select)
        print("Seq:", seq)
        print("Available scores:", metric_dict)

    metric_avg = np.mean(metric_list)
    
    return metric_avg

def get_substring_count(seq, search_list, third_base_pos=51) :
    """ 
    count occurrences of strings in search_list in seq
    """
    left_string = seq[:(third_base_pos-3)]
    rght_string = seq[third_base_pos:]
    context_string = left_string+rght_string
    
    count = 0
    for substring in search_list :
        count += context_string.count(substring)
    
    return count

def get_equal_pop_bin_edges (s, num_bins) :
    """
    from distribution s, generate bins of equal population, 
    return edges of bins
    """
    s_binned, bins = pd.qcut(s, num_bins, retbins=True, duplicates="drop")
    #substract an amount from the first entry so that all
    # values are binned (function assumes right inclusion)
    bins[0] = bins[0]-1
    
    s_binned = pd.cut(s, bins=bins, labels=False, retbins=False)

    print("Bin boundaries ("+str(num_bins)+"):", bins)
    
    return s_binned, bins

def calc_joint_counts (codon_df,
                       amino_acid_list,
                       amino_acid_col,
                       num_pos,
                       num_metric_bins,
                       variable_col,
                       variable_func,
                       **func_kwargs) :
    """
    calculate counts of joint distribution of central codons, nucleotides,
    and variable
    """

    print("Parsing", variable_col, "column into", num_metric_bins, "bins")
    
    #Add column corresponding to variable calculation if necessary
    if variable_func is not None :
        codon_df[variable_col+"_on_seq"] = codon_df["REF_Sequence"].apply(lambda x:
                                                                          variable_func(x,**func_kwargs))
        metric_col = variable_col+"_on_seq"

        #Store this calculated metric
        codon_metric_df = codon_df[["CHR", "POS", "REF", "NM_ID", metric_col]]
    else :
        metric_col = variable_col
        codon_metric_df = None

    #Generate bin labels
    binned_metric_col = variable_col+"_binned"
    codon_df[binned_metric_col], metric_bins = get_equal_pop_bin_edges(codon_df[metric_col],
                                                                           num_metric_bins)
    
    #---Which AA and codons ---
    num_aa = len(amino_acid_list)
    codons_in_list = [x for x in ccv.codon_aa_dict if 
                      ((ccv.codon_aa_dict[x] in amino_acid_list) or
                       (ccv.codon_aa_sub_dict[x] in amino_acid_list))]
    codons_in_list_indx_dict = {x:i for i,x in enumerate(codons_in_list)}
    num_codons = len(codons_in_list)
    num_bases = ccv.num_bases

    #---Storage dict---
    count_aa_codon_nx_pi_var = {}

    #---Group records Amino Acid x Codon x binned Variable---
    print("Proceeding on this dataframe:")
    print(codon_df.head())
    codon_xAAxCodonxVar = codon_df.groupby([amino_acid_col, "REF_Codon", binned_metric_col])

    #---Calculate counts---
    # split reference sequence in array of letters
    codon_df["REF_Sequence_array"] = codon_df["REF_Sequence"].apply(lambda x: list(x))
    
    # loop through amino acids in list
    for a_i, amin in enumerate(amino_acid_list) :
        #how many codons = AA
        syn_codons = ccv.aa_comb_syn_codons[amin]
        num_syn_codons = len(syn_codons)
        #set up storage
        #- count:
        # Codon x Base observed x At position
        count_codon_nx_pi_var = np.zeros((num_syn_codons, 
                                          num_bases,
                                          num_pos,
                                          num_metric_bins))
        # loop through codons corresponding to this aa
        for c_i, codon in enumerate(syn_codons) :
            for var_j in range(num_metric_bins) :
                #get the index for this AA, Codon combination
                loc_aa_codon_var = (amin, codon, var_j)
                #get the indexes of all matching codon records
                #-check if this combination exists
                if loc_aa_codon_var in codon_xAAxCodonxVar.groups :
                    group_aa_codon_var = codon_xAAxCodonxVar.groups[loc_aa_codon_var]
                else :
                    #skip if no data
                    continue
                #stack the context arrays
                array_codon_pos = np.vstack(codon_df.loc[group_aa_codon_var, "REF_Sequence_array"])
                #convert array to data frame
                df_codon_pos = pd.DataFrame(array_codon_pos)
                # loop through positions in sequence context
                for p_i, pos in enumerate(df_codon_pos.columns) :
                    #for each (position) column, count the number of observations
                    # for each base
                    base_pos_counts = df_codon_pos[pos].value_counts()
                    #loop through each base
                    for base in base_pos_counts.index :
                        #store in overall array
                        count_codon_nx_pi_var[c_i, ccv.base_dict[base],
                                              p_i,
                                              var_j] = base_pos_counts[base]
        count_aa_codon_nx_pi_var[amin] = count_codon_nx_pi_var
    
    return count_aa_codon_nx_pi_var, metric_bins, codon_metric_df

def calc_cond_mutual_information_codon_nuc_var (counts_dict, third_base_pos) :
    """
    calculate conditional mutual information using joint distribution counts
    """
    #---Storage---
    cond_mut_inf_list = []
    
    #Assume order of <codon index>, <base index>, <position index>, <metric index>
    for amin in counts_dict.keys() :
        #Grab full joint table
        count_codon_nx_pi_var = counts_dict[amin]
        #-label dimensions
        count_array_shape = count_codon_nx_pi_var.shape
        num_syn_codons = count_array_shape[0]
        num_bases = count_array_shape[1]
        num_pos = count_array_shape[2]
        num_metric_bins = count_array_shape[3]

        #Calculate probabilities
        #-N(AA)
        count_aa = np.sum(count_codon_nx_pi_var[:,:,0,:])
        #-N(Var)
        count_var = np.sum(count_codon_nx_pi_var[:,:,0,:], axis=(0,1)).reshape(num_metric_bins,1)
        #-N(C, Var)
        #-shape codon x metric_bins
        count_c_var = np.sum(count_codon_nx_pi_var[:,:,0,:], axis=1)
        #-N(N_x, Var)
        #-shape base x pos x metric_bins
        count_nx_pi_var = np.sum(count_codon_nx_pi_var, axis=0)
        #-N(C, N_x, Var)
        #-shape codon x base x position x metric_bins
        #count_codon_nx_pi_var

        p_var = count_var/count_aa
        p_c_var = count_c_var/count_aa
        p_nx_pi_var = count_nx_pi_var/count_aa
        p_c_nx_pi_var = count_codon_nx_pi_var/count_aa
    
        mut_inf_lessVar = np.zeros(num_pos)
        for p_i in range(num_pos) :
            for c_i in range(num_syn_codons) :
                for nx_i in range(num_bases) :
                    for var_i in range(num_metric_bins) :
                        if p_c_nx_pi_var[c_i,nx_i, p_i, var_i] == 0 :
                            continue
                        else :
                            log_1 = np.log(p_c_nx_pi_var[c_i, nx_i, p_i, var_i])
                            log_2 = np.log(p_var[var_i])
                            log_3 = np.log(p_c_var[c_i, var_i])
                            log_4 = np.log(p_nx_pi_var[nx_i, p_i, var_i])
                            log_arg = log_1 + log_2 - log_3 - log_4
                            summand = p_c_nx_pi_var[c_i, nx_i, p_i, var_i]*log_arg
                            mut_inf_lessVar[p_i] += summand
        mut_inf_lessVar[(third_base_pos-3):third_base_pos] = np.nan
        cond_mut_inf_list.append(mut_inf_lessVar)

    #Combine to DataFrame
    cond_mut_inf_df = pd.DataFrame(cond_mut_inf_list,
                                   index=counts_dict.keys(),
                                   columns=["cmi_p"+str(n) for n in range(num_pos)])
    
    return cond_mut_inf_df

def calc_mutual_information_codon_var (counts_dict) :
    """
    calculate mutual information between codon distribution and 
    variable distribution
    """
    #---Storage---
    cond_mut_inf_list = []
    
    #Assume order of <codon index>, <base index>, <position index>, <metric index>
    for amin in counts_dict.keys() :
        #Grab full joint table
        count_codon_nx_pi_var = counts_dict[amin]
        #-label dimensions
        count_array_shape = count_codon_nx_pi_var.shape
        num_syn_codons = count_array_shape[0]
        num_bases = count_array_shape[1]
        num_pos = count_array_shape[2]
        num_metric_bins = count_array_shape[3]
        
        #Calculate probabilities
        #-N(AA)
        count_aa = np.sum(count_codon_nx_pi_var[:,:,0,:])
        #-N(C)
        count_c = np.sum(count_codon_nx_pi_var[:,:,0,:], axis=(1,2)).reshape(num_syn_codons,1)
        #-N(Var)
        count_var = np.sum(count_codon_nx_pi_var[:,:,0,:], axis=(0,1)).reshape(num_metric_bins,1)
        #-N(C, Var)
        #-shape codon x metric_bins
        count_c_var = np.sum(count_codon_nx_pi_var[:,:,0,:], axis=1)

        p_c = count_c/count_aa
        p_var = count_var/count_aa
        p_c_var = count_c_var/count_aa
    
        mut_inf_codon_var = 0.0
        for c_i in range(num_syn_codons) :
            for var_i in range(num_metric_bins) :
                if p_c_var[c_i, var_i] == 0:
                    continue
                else :
                    log1 = np.log(p_c_var[c_i, var_i])
                    log2 = np.log(p_var[var_i])
                    log3 = np.log(p_c[c_i])
                    log_arg = log1 - log2 - log3
                    summand = p_c_var[c_i, var_i]*log_arg
                    mut_inf_codon_var += summand
        cond_mut_inf_list.append(mut_inf_codon_var)

    #Combine to DataFrame
    cond_mut_inf_df = pd.DataFrame(cond_mut_inf_list,
                                   index=counts_dict.keys(),
                                   columns=["mi"])
    
    return cond_mut_inf_df

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--codon-file', '-f', help="Path to TSV containing Codon and Sequence Context information",
                        required=True)
    parser.add_argument('--aa-list', '-l', help="Flag which amino acid list to use (all, at_cp3), or enter list",
                        nargs="+",
                        default="all")
    parser.add_argument('--aa-column', '-a', help="Name of column in dataframe containing the selected amino acid labels",
                        type=str,
                        default="REF_AminoAcid")
    parser.add_argument('--num-pos', '-p', help="Number of positions in sequence context field",
                        type=int,
                        default=101)
    parser.add_argument('--num-bins', '-b', help="Number of bins to split the variable's distribution into",
                        type=int,
                        default=20)
    parser.add_argument('--var-column', '-v', help="Name of column in dataframe (existing or added) containing the sequence or metric to condition on",
                        type=str)
    parser.add_argument('--var-function', '-vf', help="Select option for calculating final metric on var-column. 1) substring_count, 2) codon_metric, 3) mapped_metric.",
                        default=None,
                        choices=["substring_count", "codon_metric", "mapped_metric", None])
    parser.add_argument('--metric-file', '-mf', help="File to path with key column(s) and then a metric column",
                        type=str)
    parser.add_argument('--metric-window', '-mw', help="Number of codons flanking central codon to average metric over.",
                        type=int,
                        default=12)
    parser.add_argument('--string-list', '-sl', help="List of substrings to search for",
                        nargs="+")
    parser.add_argument('--output-folder', '-o', help="Path/ to label output files with",
                        default="../../data/1_mutual_information/")
    parser.add_argument('--output-tag', '-t', help="Tag to attach to output files",
                        type=str,
                        default="")
    args = parser.parse_args()

    context_df = pd.read_csv(args.codon_file,
                             delimiter="\t",
                             dtype={"CHR":str})

    if args.aa_list == ["all"] :
        aa_list = ccv.aa_with_syn_noSTOP
    elif args.aa_list == ["at_cp3"] :
        aa_list = ccv.aa_sub_with_syn_noSTOP
    else :
        aa_list = args.aa_list

    print("Processing these amino acids:", aa_list)

    #Set function to calculate metric on context
    # add corresponding arguments
    third_base_pos = int((args.num_pos-1)/2 + 1)
    print("We're assuming a third_base_pos of:", third_base_pos)
        
    if args.var_function == "substring_count" :
        var_function = get_substring_count
        function_kwargs = {"search_list":args.string_list,
                           "third_base_pos":third_base_pos}
        print("-counting the number of occurrences of:", args.string_list)
    elif args.var_function == "codon_metric" :
        var_function = get_codon_metric_avg
        if args.metric_file is not None :
            metric_df = pd.read_csv(args.metric_file,
                                    sep="\t", index_col=0)
            metric_dict = metric_df.squeeze()
            print("Using this codon scoring dictionary:", metric_dict)
        else :
            print ("We'd expect a file to be provided to assess codons in the context.")
        function_kwargs = {"metric_dict":metric_dict,
                           "window_size":args.metric_window,
                           "third_base_pos":third_base_pos}
    elif args.var_function == "mapped_metric" :
        #read in file
        metric_df = pd.read_csv(args.metric_file,
                                sep="\t")
        print(metric_df.head())
        key_columns = metric_df.columns[:-1]
        if metric_df.columns[-1] != args.var_column :
            print("Did you mean to analyze this column?", args.var_column)
        else :
            pass
        if len(key_columns) == 1 :
            key_columns = key_columns[0]
        else :
            pass
        context_df = context_df.merge(metric_df,
                                      on=key_columns,
                                      how="left")
        var_function = None
        function_kwargs = {}
    else :
        var_function = None
        function_kwargs = {}
        
    #Generate counts dictionary
    counts_dict, metric_bins, codon_metric_df = calc_joint_counts (context_df,
                                                                   aa_list,
                                                                   args.aa_column,
                                                                   args.num_pos,
                                                                   args.num_bins,
                                                                   args.var_column,
                                                                   var_function,
                                                                   **function_kwargs)
    
    #Calculate conditional mutual information I(C,N_i|X)
    cmi_df = calc_cond_mutual_information_codon_nuc_var (counts_dict, third_base_pos)
    #Calculate mutual information I(C,X)
    mi_var_df = calc_mutual_information_codon_var (counts_dict)
    
    #Save dataframes
    print("Saving data frames.")
    metric_filename = args.output_folder+"site_metric"+args.output_tag+".tsv"
    cond_mut_inf_filename = args.output_folder+"cond_mut_inf_codon_nuc_pos_var"+args.output_tag+".tsv"
    mut_inf_var_filename = args.output_folder+"mut_inf_codon_var"+args.output_tag+".tsv"

    if codon_metric_df is not None :
        print ("Writing:", metric_filename)
        codon_metric_df.to_csv(metric_filename, sep="\t", index=False)
    else :
        pass

    cmi_df.to_csv(cond_mut_inf_filename,
                  sep="\t", index=True)
    mi_var_df.to_csv(mut_inf_var_filename,
                     sep="\t", index=True)
