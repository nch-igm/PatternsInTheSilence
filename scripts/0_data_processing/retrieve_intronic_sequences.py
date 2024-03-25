#!/usr/bin/env python

"""
From section of intronic range, write sliding window of intronic sub-sequences 
 with some basic annotation of the central "codon"/trinucleotide
 
Usage: python retrieve_intronic_sequences.py --transcript-file <transcript_file> \
    --intron-coord-file <intron_coord_file> --reference <ref_fasta> \
        --output <output_path> --num-pos <window_size> \
            --boundary <bp_bnd> --displacement <window_slide>
Reads: transcipt_file, intron_coord_file, ref_fasta
Writes: output_path

"""

import argparse
import importlib.util

import pandas as pd
from pyfaidx import Fasta

ccv_spec = importlib.util.spec_from_file_location("codon_context_variables", "../codon_context_variables.py")
ccv = importlib.util.module_from_spec(ccv_spec)
ccv_spec.loader.exec_module(ccv)

base_list = ["A", "C", "G", "T"]

CHR_list = [str(x) for x in range(1,22+1)] + ["X", "Y"]

def get_CHR_from_UCSC (x) :
    try :
        y = int(x.split(".")[0].split("_")[1])
        z = CHR_list[y-1]
    except:
        z = x
        
    return z

comp_dict = {"A":"T",
                 "T":"A",
                 "C":"G",
                 "G":"C",
                 "N":"N"}

def reverse_complement (seq) :
    """
    get reverse complement of seq
    """
    seq_comp = [comp_dict[a] for a in seq]
    seq_rev_comp = "".join(seq_comp[::-1])
    
    return seq_rev_comp
    
def get_nonover_subsequences (seq_len, seq_bnd, disp, df, ref_dict,
                              write_out) :
    """
    from single intron sequence, generate sub-sequencess of length seq_len,
    omitting seq_bnd base pairs from intron ends, sliding the sequence window
    disp bases and writing out codon records
    """
#starts attributes are 1-based
#end attributes are 0-base
#reference_dict['1'][982939:982942] will print entries 982940, 982941, 982942
    ref_pos = int((seq_len-1)/2)

    num_introns = df.shape[0]
    for i in range(num_introns) :
        intron_chr = df.iloc[i]["CHR"]
        intron_left = int(df.iloc[i]["Intron_Start"])+seq_bnd
        intron_right = int(df.iloc[i]["Intron_Stop"])-seq_bnd
        intron_sense = df.iloc[i]["Sense"]
    
        intron_seq = ref_dict[intron_chr][(intron_left-1):(intron_right)].seq.upper()    
        #Check for too-short sequence or non-standard bases
        if (len(intron_seq) < seq_len or 
            len(set(intron_seq).difference(base_list)) > 0):
            continue
        else : 
            pass
        
        if intron_sense == "-" :
            intron_seq = reverse_complement(intron_seq)
        else :
            pass
        
        #Remove last element if not perfectly divisible
        subseq_list = [intron_seq[x:(x+seq_len)] for x 
                       in range(0, len(intron_seq)-(seq_len-disp), disp)]
        if len(subseq_list[-1]) != seq_len :
            subseq_list = subseq_list[:-1]
        else :
            pass
        
        for j in subseq_list :
            REF_Base = j[ref_pos]
            REF_Codon = j[(ref_pos-2):(ref_pos+1)]
            np1 = j[ref_pos+1]
            amin = ccv.codon_aa_dict[REF_Codon]
            amin_sub = ccv.codon_aa_sub_dict[REF_Codon]
            if intron_sense == "+" :
                REF = REF_Base
            else :
                REF = comp_dict[REF_Base]
                
            line_out = "\t".join([intron_chr,
                                  REF_Base,
                                  REF_Codon,
                                  amin,
                                  amin_sub,
                                  np1,
                                  intron_sense,
                                  j])
            write_out.write(line_out+"\n")
    
    return None

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcript-file", "-t", help="Path to data frame with selected transcript info")
    parser.add_argument("--intron-coord-file", "-ic", help="Path to data frame with the genomic coordinates for the exons")
    parser.add_argument("--reference", "-ref", help="Path to reference genome Fasta file")
    parser.add_argument("--output", "-o", help="Path to write output data frame to")
    parser.add_argument('--num-pos', '-p', help="Number of positions in sequence context",
                        type=int,
                        default=101)
    parser.add_argument('--boundary', '-bnd', help="Number of nucleotides near intron boundary to exclude",
                        type=int,
                        default=50)
    parser.add_argument('--displacement', '-d', help="Number of bases to slide window over for neighboring sub-sequences",
                        type=int,
                        default=50)
    
    args = parser.parse_args()
    
    sequence_length = args.num_pos
    boundary = args.boundary
    print("Writing subsequences of length", sequence_length,
          "nt and a minimum of", boundary,
          "nt from the splice junction")

    #Read in files
    transcript_select_df = pd.read_csv(args.transcript_file,
                                       sep="\t",
                                       dtype={"CHR":str})
    intron_coord_df = pd.read_csv(args.intron_coord_file,
                                  sep="\t")
    #-convert chromosome format
    intron_coord_df["CHR"] = intron_coord_df["UCSC_CHR"].apply(get_CHR_from_UCSC)

    print(intron_coord_df.head())

    #Match transcript range annotations with selected transcripts
    #-not a unique key on the right, but check if selection requires more specification
    intron_select_coord_df = transcript_select_df.merge(intron_coord_df,
                                                        on=["NM_ID", "CHR", "Sense"],
                                                        how="left").reset_index()
    
    #Check that we didn't multi-map any chromosomes
    # Ex: NM_004536.2 mapped twice to Chr5
    print("Size of selected transcript data frame before merge:", transcript_select_df.shape)
    print("-> size after merge:", intron_select_coord_df.shape)
    print("Top transcription mapping (=1?):", intron_select_coord_df.groupby(["NM_ID", "Intron_Number"]).count().max()["CHR"])

    #Check for NA values
    print("Transcripts without (matched) introns?:", sum(intron_select_coord_df["Intron_Number"].isna()))
    
    intron_select_coord_df = intron_select_coord_df[~intron_select_coord_df["Intron_Number"].isna()]
    
    #Generate dict for genome
    reference_dict = Fasta(args.reference)
    
    with open(args.output,"w") as intron_seq_out :
        #Write header
        header_line = "\t".join(["CHR",
                                 "REF_Base",
                                 "REF_Codon",
                                 "REF_AminoAcid",
                                 "REF_AminoAcid_sub",
                                 "N+1",
                                 "Sense",
                                 "REF_Sequence"])
        intron_seq_out.write(header_line+"\n")
        get_nonover_subsequences(sequence_length, boundary, args.displacement,
                                 intron_select_coord_df,
                                 reference_dict, intron_seq_out)