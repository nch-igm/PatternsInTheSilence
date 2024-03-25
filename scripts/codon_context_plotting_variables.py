#!/usr/bin/env python3

#Updated: June 29, 2022

#import seaborn as sns
#import matplotlib.pyplot as plt
import matplotlib
print(matplotlib.get_cachedir())
print(matplotlib.font_manager.findfont("Arial"))
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

###Colors
#Nucleotides
nuc_colors = {"A":"#e66101",
              "T":"#fdb863",
              "C":"#813c99",
              "G":"#cfabd2"}

#Sequence context descriptors
context_variable_colors = {'ApTcount': '#dae6f1',
                           'TpAcount': '#b3cede',
                           'Ccount': '#78aac8',
                           'CpGcount': '#4884af',
                           'GCcount': '#225b91',
                           'end': '#f7ded2',
                           'cd': '#efb5a0',
                           'cfe': '#e88d75',
                           'meafe': '#dd6150',
                           'mfe': '#bf3838',
                           'efe': '#972328',
                           'tAIavg': '#c9cade',
                           'CSCavg': '#7f77aa'}

aa_polarity_colors = {"Nonpolar":"#d7d7d7",
                      "Polar":"#cdeafa",
                      "Basic":"#c8fa94",
                      "Acidic":"#f77cdc"}

aa_degeneracy_colors = {2:"#ecdbff",
                        3:"#a88fcc",
                        4:"#644a9c",
                        6:"#0a096d"}

aa_sub_colors = {"F":"#add36c",
                 "Y":"#c3bcc8",
                 "H":"#d8beb5",
                 "Q":"#ed9ab3",
                 "N":"#83d4b4",
                 "K":"#eaac5f",
                 "D":"#88c0ea",
                 "E":"#ccc958",
                 "C":"#d3ade8",
                 "I":"#49382a",
                 "V":"#ba506b",
                 "P":"#608388",
                 "T":"#c95030",
                 "A":"#7378be",
                 "G":"#9a793c",
                 "S2":"#f47fe5",
                 "S4":"#c84eb8",
                 "L2":"#82db69",
                 "L4":"#538c40",
                 "R2":"#8f97ff",
                 "R4":"#3452d0"}

aa_colors = {"F":"#add36c",
             "Y":"#c3bcc8",
             "H":"#d8beb5",
             "Q":"#ed9ab3",
             "N":"#83d4b4",
             "K":"#eaac5f",
             "D":"#88c0ea",
             "E":"#ccc958",
             "C":"#d3ade8",
             "I":"#49382a",
             "V":"#ba506b",
             "P":"#608388",
             "T":"#c95030",
             "A":"#7378be",
             "G":"#9a793c",
             "S":"#c84eb8",
             "L":"#538c40",
             "R":"#3452d0"}

codon_by_opt_colors = {}

snv_colors = {"A>C":"#e66101",
              "A>G":"#df905d",
              "A>T":"#FFF4BC",
              "T>A":"#fdb863",
              "T>C":"#ab683c",
              "T>G":"#983200",
              "C>A":"#5e3c99",
              "C>G":"#986ab9",
              "C>T":"#cd9cda", 
              "CpG>CpA":"#ffd2ff",
              "CpG>TpG":"#9890e1",
              "G>A":"#b2abd2",
              "G>C":"#747aae",
              "G>T":"#345eff"}

###Functions to formatting axis labels
def get_flat_sequence_labels_plain (seq_length, cp3_index) :
    left_length = (cp3_index-3)+1
    right_length = seq_length-(cp3_index+1)
    
    left_labels = ["N-"+str(x) for x in range(1,left_length+1)[::-1]]
    codon_labels = ["CP1", "CP2", "CP3"]
    rght_labels = ["N+"+str(x) for x in range(1,right_length+1)]
    
    labels = left_labels + codon_labels + rght_labels 
    
    return labels

def get_flat_sequence_labels_formatted (seq_length, cp3_index) :
    left_length = (cp3_index-3)+1
    right_length = seq_length-(cp3_index+1)
    
    left_labels = [f"$N_{{-{x}}}$" for x in range(1,left_length+1)[::-1]]
    codon_labels = ["CP1", "CP2", "CP3"]
    rght_labels = [f"$N_{{+{x}}}$" for x in range(1,right_length+1)]
    
    labels = left_labels + codon_labels + rght_labels 
    
    return labels

def get_flat_codon_sequence_labels_formatted (seq_length, central_index) :
    left_length = (central_index-1)+1
    right_length = seq_length-(central_index+1)
    
    left_labels = [f"$C_{{-{x}}}$" for x in range(1,left_length+1)[::-1]]
    codon_labels = ["Central codon"]
    rght_labels = [f"$C_{{+{x}}}$" for x in range(1,right_length+1)]
    
    labels = left_labels + codon_labels + rght_labels 
    
    return labels
