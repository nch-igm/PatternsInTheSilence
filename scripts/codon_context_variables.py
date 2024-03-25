#!/usr/bin/env python3

#Updated: June 22, 2021

###Codon variables
#number of nucleotides per codon
codon_size = 3
#number of available bases (A,T,C,G)
num_bases = 4

###Bases and base pairing
#relate base to index
base_dict = {'A':0, 'T':1, 'C':2, 'G':3}
#relate index to base
indx_dict = {'0':'A', '1':'T', '2':'C', '3':'G'}
#relate base to complement
complement_dict = {"A":"T", "C":"G", "T":"A", "G":"C"}

###Going between codons and amino acids
#Relate codon sequence to AA symbol
codon_aa_dict = {'TTT':'F', 'TTC':'F', 
                  'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 
                  'CTA':'L', 'CTG':'L', 
                  'ATT':'I', 'ATC':'I', 'ATA':'I',
                  'ATG':'M',
                  'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                  'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                  'AGT':'S', 'AGC':'S',
                  'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                  'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                  'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                  'TAT':'Y', 'TAC':'Y',
                  'TAA':'STOP', 'TAG':'STOP', 'TGA':'STOP',
                  'CAT':'H', 'CAC':'H', 
                  'CAA':'Q', 'CAG':'Q',
                  'AAT':'N', 'AAC':'N',
                  'AAA':'K', 'AAG':'K',
                  'GAT':'D', 'GAC':'D',
                  'GAA':'E', 'GAG':'E',
                  'TGT':'C', 'TGC':'C',
                  'TGG':'W', 
                  'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                  'AGA':'R', 'AGG':'R', 
                  'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

#Relate codon to 3rd-position-variant AA subclass symbolc
codon_aa_sub_dict = {'TTT':'F', 'TTC':'F', 
                     'TTA':'L2', 'TTG':'L2', 'CTT':'L4', 'CTC':'L4', 
                     'CTA':'L4', 'CTG':'L4', 
                     'ATT':'I', 'ATC':'I', 'ATA':'I',
                     'ATG':'M',
                     'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                     'TCT':'S4', 'TCC':'S4', 'TCA':'S4', 'TCG':'S4',
                     'AGT':'S2', 'AGC':'S2',
                     'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                     'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                     'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                     'TAT':'Y', 'TAC':'Y',
                     'TAA':'STOP2', 'TAG':'STOP2', 'TGA':'STOP1',
                     'CAT':'H', 'CAC':'H', 
                     'CAA':'Q', 'CAG':'Q',
                     'AAT':'N', 'AAC':'N',
                     'AAA':'K', 'AAG':'K',
                     'GAT':'D', 'GAC':'D',
                     'GAA':'E', 'GAG':'E',
                     'TGT':'C', 'TGC':'C',
                     'TGG':'W', 
                     'CGT':'R4', 'CGC':'R4', 'CGA':'R4', 'CGG':'R4',
                     'AGA':'R2', 'AGG':'R2', 
                     'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

###Amino Acid sets
#All AA symbols
aa_labels = ['F', 'L', 'I', 'M', \
                 'V', 'S', 'P', 'T', \
                 'A', 'Y', 'STOP', \
                 'H', 'Q', 'N', 'K', 'D', \
                 'E', 'C', 'W', 'R', 'G']
#All AA_sub symbols
aa_sub_labels = ['F', 'L2', 'L4', 'I', 'M', \
                 'V', 'S4', 'S2', 'P', 'T', \
                 'A', 'Y', 'STOP1', 'STOP2', \
                 'H', 'Q', 'N', 'K', 'D', \
                 'E', 'C', 'W', 'R4', 'R2', 'G']
#All AA classes with degeneracy and no Stop codons
aa_with_syn_noSTOP  = ['F', 'L', 'I', 'V',
                       'S', 'P', 'T', 'A', 'Y', 
                       'H', 'Q', 'N', 'K', 'D',
                       'E', 'C', 'R', 'G']

#All AA_sub classes with degeneracy
aa_sub_with_syn = ['F', 'L2', 'L4', 'I', 'V',
                   'S4', 'S2', 'P', 'T', 'A', 'Y', 
                   'STOP1', 'H', 'Q', 'N', 'K', 'D',
                   'E', 'C', 'R4', 'R2', 'G']
#All AA_sub classes with degeneracy and no Stop codons
aa_sub_with_syn_noSTOP = ['F', 'L2', 'L4', 'I', 'V',
                          'S4', 'S2', 'P', 'T', 'A', 'Y',
                          'H', 'Q', 'N', 'K', 'D',
                          'E', 'C', 'R4', 'R2', 'G']
#AA_sub classes where C isn't observed in 3rd position
aa_sub_no_C3 = ['L2', 'M', 'STOP1', 'STOP2', 'Q', 'K', 'E', 'W', 'R2']
#AA_sub classes that don't appear in normal set of AA labels
aa_sub_only = ['L2', 'L4', 'S2', 'S4', 'R2', 'R4', 'STOP1', 'STOP2']

#Relate each AA symbol to its set of corresponding codons
#AA -> [codons]
aa_syn_codons = {}
for aa in aa_labels :
    aa_syn_codons[aa] = [y for y in codon_aa_dict
                         if codon_aa_dict[y] == aa]
#AA_sub -> [codons]
aa_sub_syn_codons = {}
for aa_sub in aa_sub_labels :
    aa_sub_syn_codons[aa_sub] = [y for y in codon_aa_sub_dict
                                 if codon_aa_sub_dict[y] == aa_sub]
#Combine the above two dictionaries
aa_comb_syn_codons = aa_syn_codons.copy()
for aa_sub in aa_sub_only :
    aa_comb_syn_codons[aa_sub] = aa_sub_syn_codons[aa_sub]
    
#Amino Acids by multiplicity
aa_sub_mult2 = ["F", "L2", "S2", "Y", "H", "Q", "N", "K", "D", "E", "C", "R2"]
aa_mult3 = ["I"]
aa_sub_mult4 = ["L4", "V", "S4", "P", "T", "A", "R4", "G"]

aa_mult6 = ["L", "R", "S"]
aa_mult2 = ["F", "Y", "H", "Q", "N", "K", "D", "E", "C"]
aa_mult4 = ["V", "P", "T", "A", "G"]

###Index codons
#Get list of all codons
codons = codon_aa_dict.keys()
codons_with_sub_syn_noSTOP = [x for x in codon_aa_sub_dict
                              if codon_aa_sub_dict[x] in aa_sub_with_syn_noSTOP]
##Assign an index to each codon
codons_index_dict = {codon:c_i for c_i, codon
                     in enumerate(codons)}
codons_with_sub_syn_noSTOP_index_dict = {codon:c_i for c_i, codon
                                         in enumerate(codons_with_sub_syn_noSTOP)}
#Codons that code for common STOP codons
codons_STOP = ['TAA','TAG','TGA']

###Mapping codon substitutions
def codonchange_to_snv (codonchange) :
    #x in format of REF_Codon>ALT_Codon +/- CpG context
    ref_alt = codonchange.split(">")
    ref_codon = ref_alt[0][0:3]
    
    ref_allele = ref_alt[0][2:]
    alt_allele = ref_alt[1][2:]
    
    if (ref_allele == "G" and 
        alt_allele == "A" and
        ref_codon[1] == "C") :
        ref_allele = "CpG"
        alt_allele = "CpA"
    else :
        pass
    
    snv_context = ref_allele+">"+alt_allele
    return snv_context

def codonchange_to_aa_sub (codonchange) :
    #x in format of REF_Codon>ALT_Codon +/- CpG context
    ref_codon = codonchange[0:3]
    ref_aa_sub = codon_aa_sub_dict[ref_codon]
    
    return ref_aa_sub
    
    
    

