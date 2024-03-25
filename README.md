# Setup

This walks through scripts and data published with Moreland _et al._ (2023).

## Conda environment

```bash
conda create --name codon_context python=3.8 scipy=1.9.3 statsmodels=0.13.2 pandas=1.5.1 numpy=1.23 matplotlib=3.6 mygene seaborn pip 
pip install pyfaidx
conda install -c bioconda gseapy
pip install jupyterlab
conda install ipykernel
conda install -c conda-forge scikit-learn=1.2.0
```

# 0. Data processing

## 0.1 From SNV parquet generate transcript list

Parquet file contains entries for all possible SNVs variants generated on loci in RefSeq transcripts.

```bash
n1_get_transcript_list.ipynb

# Input: 
#   <Gaither et al. parquet>
# Output:
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_summary.tsv
```

## 0.2 Get Entrez gene annotations

Get Entrez Gene IDs that correspond with the RefSeq transcript names.

```bash
#From scripts/0_data_processing/
python annotate_transcript_list_with_entrez.py \
    ../../data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_summary.tsv \
    ../../data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_summary.tsv

# Input:
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_summary.tsv
# Output:
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_summary.tsv
```

## 0.3 Filter on transcript attributes

Filter transcripts to select one representative transcript per gene.

```bash
n2_select_transcripts.ipynb

# Input:
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_summary.tsv
#   data/0_data_processing/MANE.GRCh38.v0.9.summary.txt.gz
#   <Gaither et al. parquet>
# Output:
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_flagged_selected_summary.tsv
```

## 0.4 Filter on variant attributes

Filter variants so that records don't overlap, annotate records with additional columns, and write codon and variant tables.

```bash
n3_filter_variants_and_annotate.ipynb

# Input:
#   <Gaither et al. parquet>
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_flagged_selected_summary.tsv
# Output:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_synonymous.tsv
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_CP3.tsv
```

## 0.5 Filter synonymous variants

Filter codon records, apply to synonymous variant records, write out final codon table, write per-Amino acid variant tables.

```bash
0a_filter_syn_codons_for_analysis.ipynb

# Input:
#   data/0_data_processing/transcripts/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_CP3.tsv
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_summary.tsv
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_synonymous.tsv
# Output:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_CP3.tsv
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_selected_0a.tsv
#   data/0_data_processing/rna_stability_exports/byAminoAcid_sub/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_synonymous_noEdge_noSTOP_CP3_seqCol_AminoAcid<AA>.tsv
#   data/0_data_processing/rna_stability_exports/byAminoAcid_sub/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_synonymous_noEdge_noSTOP_CP3_nonseqCol_AminoAcid<AA>.tsv
```

## 0.6 Sub-sample records by amino acid sub-class

Generate version of codon record table where each amino acid sub-class is sampled to the same size.

```bash
python generate_subsampled_context_distribution.py \
    -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
    -l aa_sub_with_syn_noSTOP \
    -a REF_AminoAcid_sub \
    -o ../../data/0_data_processing/rna_stability_exports/ \
    -t RNAStability_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3_minSubsample

# Input:
#   data/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/1_mutual_information/RNAStability_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3_minSubsample.tsv
```

## 0.7 Generate records from intron sequences

Go through transcripts and write out genomic coordinates of the exons.

```bash
#From scripts/0_data_processing/
python annotate_transcripts_with_exons_genomic_range.py \
    -g ../../data/0_data_processing/GCF_000001405.33_knownrefseq_alignments.gff3 \
    -f ../../data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments

# Input:
#   data/external/GCF_000001405.33_knownrefseq_alignments.gff3
# Output:
#   data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments_Exon_pos.tsv
#   data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments_Trans_genomic_coords.tsv
```

From exon coordinates, infer bounds of intronic ranges, and write table of intronic ranges.

```bash
python write_intronic_ranges.py \
    ../../data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments_Exon_pos.tsv \
    ../../data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments_Intron_pos.tsv

# Input:
#   data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments_Exon_pos.tsv
# Output:
#   data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments_Intron_pos.tsv
```

From intronic ranges, generate set of intronic sub-ranges by sliding a window by a certain displacement, and write out table of intronic sequences from sub-ranges with some annotations.

```bash
python retrieve_intronic_sequences.py \
    -t ../../data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_selected_0a.tsv \
    -ic ../../data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments_Intron_pos.tsv \
    -ref ../../data/0_data_processing/hs38DH.nochr.fa \
    -o ../../data/0_data_processing/transcripts/intron_subsequences_101nt_50bnd_20w_transcript_selected.tsv \
    -bnd 50 \
    -p 101 \
    -d 20

# Input:
#   data/0_data_processing/transcripts/RNAStability_v10.5.1_transcript_entrez_selected_0a.tsv
#   data/0_data_processing/transcripts/GCF_000001405.33_knownrefseq_alignments_Intron_pos.tsv
#   data/0_data_processing/hs38DH.nochr.fa
# Output:
#   data/0_data_processing/transcripts/intron_subsequences_101nt_50bnd_20w_transcript_selected.tsv
```

## 0.8 Sub-sample intronic records

From set of intron-sourced records, generate a data set with the same proportions of codons as the original table.

```bash
python generate_matched_context_distribution.py \
      -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
      -l codons_with_sub_syn_noSTOP \
      -a REF_Codon \
      -f2 ../../data/0_data_processing/transcripts/intron_subsequences_101nt_50bnd_20w_transcript_selected.tsv \
      -o ../../data/0_data_processing/transcripts/ \
      -t intron_subsequences_101nt_50bnd_20w_CP3_matched_sample
# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
#   data/0_data_processing/transcripts/intron_subsequences_101nt_50bnd_20w_transcript_selected.tsv
# Output:
#   data/1_mutual_information/intron_subsequences_101nt_50bnd_20w_CP3_matched_sample.tsv
```

# 1. Mutual information

## 1.1 Calculate codon-nucleotide mutual information

From table of codon contexts, calculate codon and nucleotide distributions and mutual information.

```bash
python mutual_information_codon_nuc.py \
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -o ../../data/1_mutual_information/ \
       -t _AAsub
# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/1_mutual_information/shannon_entropy_codon_AAsub.tsv
#   data/1_mutual_information/shannon_entropy_nuc_pos_101bp_AAsub.tsv
#   data/1_mutual_information/shannon_entropy_codon_nuc_pos_101bp_AAsub.tsv
#   data/1_mutual_information/mut_info_codon_nuc_pos_101bp_AAsub.tsv
```

Calculate mutual information on the sub-sampled table.

```bash
python mutual_information_codon_nuc.py \
    -f ../../data/0_data_processing/rna_stability_exports/RNAStability_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3_minSubsample.tsv \
    -l at_cp3 \
    -a REF_AminoAcid_sub \
    -p 101 \
    -o ../../data/1_mutual_information/ \
    -t _AAsub_subsampled

# Input
#   data/0_data_processing/rna_stability_exports/RNAStability_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3_minSubsample.tsv
# Output:
#   data/1_mutual_information/mut_info_codon_nuc_pos_101bp_AAsub_subsampled.tsv
#   data/1_mutual_information/shannon_entropy_codon_AAsub_subsampled.tsv
#   data/1_mutual_information/shannon_entropy_codon_nuc_pos_101bp_AAsub_subsampled.tsv
#   data/1_mutual_information/shannon_entropy_nuc_pos_101bp_AAsub_subsampled.tsv
```

## 1.2 Calculate codon-context codon mutual information

From table of codon contexts, calculate codon distribution and distribution of codons in the sequence contexts and mutual information.

```bash
python mutual_information_codon_cxtCodon.py \
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 33 \
       -o ../../data/1_mutual_information/ \
       -t _AAsub

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/1_mutual_information/shannon_entropy_cxtCodon_33cod_AAsub.tsv
#   data/1_mutual_information/shannon_entropy_codon_cxtCodon_33cod_AAsub.tsv
#   data/1_mutual_information/mut_info_codon_cxtCodon_33cod_AAsub.tsv
```

## 1.3 Calculate mutual information on shuffled contexts

Shuffle the codon table a specified number of times and calculate the codon-nucleotide mutual information on each shuffled data set.

```bash
python mutual_information_on_shuffled_contexts.py \
    -f ../../data/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
    -l at_cp3 \
    -a REF_AminoAcid_sub \
    -p 101 \
    -ns 100 \
    --seed-start 0 \
    -o ../../data/1_mutual_information/mut_info_context_shuffled/ \
    -t _AAsub

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/1_mutual_information/mut_info_context_shuffled/mut_info_codon_nuc_pos_101bp_AAsub_shuffle<num>.tsv
```

## 1.4 Calculate mutual information on intronic sequences

With this intron-sourced table, repeat codon-nucleotide mutual information calculations.

```bash
python mutual_information_codon_nuc.py \
       -f ../../data/1_mutual_information/intron_subsequences_101nt_50bnd_20w_CP3_matched_sample.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -o ../../data/1_mutual_information/ \
       -t _AAsub_intronic_50bnd_20w

# Input:
#   data/1_mutual_information/intron_subsequences_101nt_50bnd_20w_CP3_matched_sample.tsv
# Output:
#   data/1_mutual_information/shannon_entropy_codon_AAsub_intronic_50bnd_20w.tsv
#   data/1_mutual_information/shannon_entropy_nuc_pos_101bp_AAsub_intronic_50bnd_20w.tsv
#   data/1_mutual_information/shannon_entropy_codon_nuc_pos_101bp_AAsub_intronic_50bnd_20w.tsv
#   data/1_mutual_information/mut_info_codon_nuc_pos_101bp_AAsub_intronic_50bnd_20w.tsv
```

## 1.5 Formatting tables

* Calculate average $H_A(C,N_i)$ across ranges of $i$ for each amino acid, $A$.
* Calculate summations of $H_A(C,N_i)$ per positions in codons, $\sum_{i \in c_j} H_A(C,N_i)$, compare to $H_A(C,C_i)$.
* Aggregate $H_A^{\text{shuffled}}(C,N_i)$ calculations over shuffles of codon data set

```bash
1_mutual_information_processing.ipynb
```

# 2. Conditional mutual information

## 2.1 Format external tables

Source tables for codon metrics and format them for analysis

```bash
2a_process_sequence_metric_files.ipynb

# Inputs:
#   Download tables at given source
# Outputs:
#   data/2_conditional_mutual_information/csc_wu_2019.tsv
#   data/2_conditional_mutual_information/tai_tuller_2010.tsv
```

## 2.2 Calculate conditional mutual information

For each sequence context metric, repeat calculation of conditional mutual information between codon and nucleotide distributions given the context metric.

### Nucleotide counts

G/C content

```bash
python cond_mutual_information_codon_nuc_var.py \
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_Sequence \
       -vf substring_count \
       -sl G C \
       -o ../../data/2_conditional_mutual_information/ \
       -t _GCcount_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_GCcount_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_GCcount_20bins.tsv
#   data/2_conditional_mutual_information/site_metric_GCcount_20bins.tsv
```

CpG counts

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_Sequence \
       -vf substring_count \
       -sl CG \
       -o ../../data/2_conditional_mutual_information/ \
       -t _CpGcount_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_CpGcount_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_CpGcount_20bins.tsv
#   data/2_conditional_mutual_information/site_metric_CpGcount_20bins.tsv
```

TpA counts

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_Sequence \
       -vf substring_count \
       -sl TA \
       -o ../../data/2_conditional_mutual_information/ \
       -t _TpAcount_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_TpAcount_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_TpAcount_20bins.tsv
#   data/2_conditional_mutual_information/site_metric_TpAcount_20bins.tsv
```

ApT counts

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_Sequence \
       -vf substring_count \
       -sl AT \
       -o ../../data/2_conditional_mutual_information/ \
       -t _ApTcount_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_ApTcount_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_ApTcount_20bins.tsv
#   data/2_conditional_mutual_information/site_metric_ApTcount_20bins.tsv
```

C nucleotide counts

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_Sequence \
       -vf substring_count \
       -sl C \
       -o ../../data/2_conditional_mutual_information/ \
       -t _Ccount_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_Ccount_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_Ccount_20bins.tsv
#   data/2_conditional_mutual_information/site_metric_Ccount_20bins.tsv
```

### RNA stability metrics

local predicted MFE

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_mfeValue \
       -o ../../data/2_conditional_mutual_information/ \
       -t _mfe_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_mfe_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_mfe_20bins.tsv
```

local predicted CFE

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_cfeValue \
       -o ../../data/2_conditional_mutual_information/ \
       -t _cfe_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_cfe_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_cfe_20bins.tsv
```

local predicted MEAFE

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_meafeValue \
       -o ../../data/2_conditional_mutual_information/ \
       -t _meafe_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_meafe_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_meafe_20bins.tsv
```

local predicted EFE

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_efeValue \
       -o ../../data/2_conditional_mutual_information/ \
       -t _efe_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_efe_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_efe_20bins.tsv
```

local predicted CD

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_cdValue \
       -o ../../data/2_conditional_mutual_information/ \
       -t _cd_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_cd_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_cd_20bins.tsv
```

local predicted END

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_endValue \
       -o ../../data/2_conditional_mutual_information/ \
       -t _end_20bins

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_end_20bins.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_end_20bins.tsv
```

### Codon metrics

average tAI of surrounding sequence

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_Sequence \
       -vf codon_metric \
       -mf ../../data/2_conditional_mutual_information/tai_tuller_2010.tsv \
       -mw 12 \
       -o ../../data/2_conditional_mutual_information/ \
       -t _tAIavg_12cod_20bins_221204

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_tAIavg_12cod_20bins_221204.tsv 
#   data/2_conditional_mutual_information/mut_inf_codon_var_tAIavg_12cod_20bins_221204.tsv  
#   data/2_conditional_mutual_information/site_metric_tAIavg_12cod_20bins_221204.tsv
```

average CSC of surrounding sequence

```bash
python cond_mutual_information_codon_nuc_var.py\
       -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
       -l at_cp3 \
       -a REF_AminoAcid_sub \
       -p 101 \
       -b 20 \
       -v REF_Sequence \
       -vf codon_metric \
       -mf ../../data/2_conditional_mutual_information/csc_wu_2019.tsv \
       -mw 12 \
       -o ../../data/2_conditional_mutual_information/ \
       -t _CSCavg_12cod_20bins_221204

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_CSCavg_12cod_20bins_221204.tsv
#   data/2_conditional_mutual_information/mut_inf_codon_var_CSCavg_12cod_20bins_221204.tsv 
#   data/2_conditional_mutual_information/site_metric_CSCavg_12cod_20bins_221204.tsv
```

## 2.3 Set codon context range

Use conditional mutual information between codon and nucleotide distributions, conditioned on GC content, to measure a drop-off range around central codon.

```bash
2b_measure_GC_CMI_range.ipynb

# Input:
#   data/2_conditional_mutual_information/cond_mut_inf_codon_nuc_pos_var_GCcount_20bins.tsv
# Output:
#   data/2_conditional_mutual_information/cmi_codon_nuc_pos_GCcount_20bins_range.tsv
```

## 2.4 Formatting tables

* Convert CMI tables from wide to long format, add codon position information
* Combine CMI data for three categories of sequence context descriptions
* Combine MI data into one table for all sequence context descriptions

# 3. Codon context score

## 3.1 Calculate codon bias factors

Re-generate codon and nucleotide probabilities, and calculate codon-nucleotide bias factors.

```bash
python codon_nucleotide_bias_factors.py \
    -f ../../data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv \
    -l at_cp3 \
    -a REF_AminoAcid_sub \
    -p 101 \
    -pr 12 \
    -o ../../data/3_codon_context_score/ \
    -t ""

# Input:
#   data/0_data_processing/rna_stability_exports/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_noEdge_wSyn_noSTOP_noStruct_CP3.tsv
# Output:
#   data/3_codon_context_score/codon_nuc_mi_bias_factors_AminoAcid_sub_12nt.tsv
```

## 3.2 Generate codon context scores

Go through synonymous variant records and assign weights for reference codon, alternate codon, and nucleotides in the sequence context.

```bash
mkdir ../../data/3_codon_context_score/syn_variant_codon_context_score_byAminoAcid_sub

for X in F L2 L4 I V S4 S2 P T A Y H Q N K D E C R4 R2 G
do
    python annotate_variant_table_codon_score.py \
       -f ../../data/0_data_processing/rna_stability_exports/byAminoAcid_sub/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_synonymous_noEdge_noSTOP_CP3_seqCol_AminoAcid"$X".tsv \
       -ff ../../data/3_codon_context_score/codon_nuc_mi_bias_factors_AminoAcid_sub_12nt.tsv \
       -p 101 \
       -pr 12 \
       -o ../../data/3_codon_context_score/syn_variant_codon_context_score_byAminoAcid_sub/ \
       -t _AminoAcid"$X"_seqCol
done

# Input:
#   data/0_data_processing/rna_stability_exports/byAminoAcid_sub/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_synonymous_noEdge_noSTOP_CP3_seqCol_AminoAcid"$X".tsv
#   data/3_codon_context_score/codon_nuc_mi_bias_factors_AminoAcid_sub_12nt.tsv
# Output:
#   data/3_codon_context_score/syn_variant_codon_context_score_byAminoAcid_sub/ssnv_codon_context_score_mi_12nt_AminoAcid$"X"_seqCol.tsv
#   data/3_codon_context_score/syn_variant_codon_context_score_byAminoAcid_sub/ssnv_codon_context_score_mi_annotated_12nt_AminoAcid"$X"_seqCol.tsv
```

Combine the scored variant files with the variant tables with additional information columns and combine across amino acids.

```bash
python combine_scored_variant_files.py \
   -f1 "../../data/0_data_processing/rna_stability_exports/byAminoAcid_sub/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_synonymous_noEdge_noSTOP_CP3_nonseqCol_AminoAcid*.tsv" \
   -f2 "../../data/3_codon_context_score/syn_variant_codon_context_score_byAminoAcid_sub/ssnv_codon_context_score_mi_12nt_AminoAcid*_seqCol.tsv" \
   -l at_cp3 \
   -o ../../data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3.tsv

# Input:
#   data/0_data_processing/rna_stability_exports/byAminoAcid_sub/RNAStability_v10.5.1_hg38_filterDups_noOverlaps_synonymous_noEdge_noSTOP_CP3_nonseqCol_AminoAcid*.tsv
#   data/3_codon_context_score/syn_variant_codon_context_score_byAminoAcid_sub/ssnv_codon_context_score_mi_12nt_AminoAcid*_seqCol.tsv
# Output:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3.tsv
```

## 3.3 Normalize context score and assign gnomAD status

Apply standard normalization to context score. Flag variants with adequate coverage in gnomAD and by non-zero MAF. 

```bash
python context_score_normalization_maf.py \
   -vf ../../data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3.tsv \
   -s diff_sum_context_score \
   -g REF_Codon \
   -o ../../data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon.tsv
# Input:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3.tsv
# Output:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon.tsv
```

## 3.4 Measure constraint-score correlation

Generate constraint-codon context score curves per SNVContext.

```bash
python calculate_constraint_curve.py \
    -vf ../../data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon.tsv \
    -s diff_sum_context_score_REF_Codon_zscore \
    -y y y_rand \
    -g SNVContext \
    -b 100 \
    -t 12nt_CP3_xSNVContext_y_yrand

# Input:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon.tsv
# Output:
#   data/3_codon_context_score/constraint_curve_diff_sum_context_score_REF_Codon_zscore_100bins_12nt_CP3_xSNVContext_y_yrand.tsv
```

Fit linear model to constraint-codon context score curves. One run for the gnomAD data:

```bash
python calculate_constraint_curve_fit.py \
   -gvf ../../data/3_codon_context_score/constraint_curve_diff_sum_context_score_REF_Codon_zscore_100bins_12nt_CP3_xSNVContext_y_yrand.tsv \
   -s diff_sum_context_score_REF_Codon_zscore_binned_left \
   -y y \
   -g SNVContext \
   -bm 1000 \
   -t dsCS_RCodonZ_100b_12nt_CP3_xSNVContext_y

# Input:
#   data/3_codon_context_score/constraint_curve_diff_sum_context_score_REF_Codon_zscore_100bins_12nt_CP3_xSNVContext_y_yrand.tsv
# Output:
#   data/3_codon_context_score/constraint_curve_fits_1000min_dsCS_RCodonZ_100b_12nt_CP3_xSNVContext_y.tsv
```

And another from the shuffled data:

```bash
python calculate_constraint_curve_fit.py \
   -gvf ../../data/3_codon_context_score/constraint_curve_diff_sum_context_score_REF_Codon_zscore_100bins_12nt_CP3_xSNVContext_y_yrand.tsv \
   -s diff_sum_context_score_REF_Codon_zscore_binned_left \
   -y y_rand \
   -g SNVContext \
   -bm 1000 \
   -t dsCS_RCodonZ_100b_12nt_CP3_xSNVContext_yrand

# Input:
#   data/3_codon_context_score/constraint_curve_diff_sum_context_score_REF_Codon_zscore_100bins_12nt_CP3_xSNVContext_y_yrand.tsv
# Output:
#   data/3_codon_context_score/constraint_curve_fits_1000min_dsCS_RCodonZ_100b_12nt_CP3_xSNVContext_yrand.tsv
```

## 3.5 Assess score thresholds

Negate codon context scores for specified SNVContexts.

```bash
python negate_selected_context_scores.py \
    -vf ../../data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon.tsv \
    -s diff_sum_context_score_REF_Codon_zscore \
    -g SNVContext \
    -r "C>T" "G>A" \
    -o ../../data/3_codon_context_score/ \
    -f syn_variant_12nt_codon_context_score_CP3_zCodon_SNVConNeg.tsv

# Input:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon.tsv
# Output:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon_SNVConNeg.tsv
```

Range through percentiles of the constraint score and assess how well the codon
context score classifies variants that are observed or unobserved in gnomAD. At each threshold, measure enrichment of flagged variants.

```bash
python measure_constraint_thresholds_on_variants.py \
    -vf ../../data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon_SNVConNeg.tsv \
    -s diff_sum_context_score_REF_Codon_zscore \
    -g SNVContext \
    -o ../../data/3_codon_context_score/ \
    -t diff_sum_context_score_posScaled

# Input:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon_SNVConNeg.tsv
# Output:
#   data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv
```

## 3.6 Summarize context score attributes and thresholds

* Filter variant table for variants with gnomAD coverage, add annotations
* Describe distribution of context scores by position in sequence context
* Determine mutational signatures with significant scaling between constraint and context score
* Determine thresholds for context score, per mutational signature that fit specified minima for specificity and enrichment

```bash
3_codon_context_score_processing.ipynb

# Input:
#   data/3_codon_context_score/syn_variant_codon_context_score_byAminoAcid_sub/ssnv_codon_context_score_mi_annotated_12nt_AminoAcid*_seqCol.tsv
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon.tsv
#   data/3_codon_context_score/constraint_curve_fits_1000min_dsCS_RCodonZ_100b_12nt_CP3_xSNVContext_y.tsv
#   data/3_codon_context_score/constraint_curve_fits_1000min_dsCS_RCodonZ_100b_12nt_CP3_xSNVContext_yrand.tsv
#   data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv
# Output:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_scoreAvg_xPosition_xSNVContext.tsv
#   data/3_codon_context_score/constraint_curve_fits_1000min_dsCS_RCodonZ_100b_12nt_CP3_xSNVContext_y_yrand_summary.tsv
```

# 4. TCGA analysis

List of participant and sample IDs used can be found: 
* `data/4_tcga_analysis/project_id_brca_PASS_merge.fam`
* `data/4_tcga_analysis/project_id_ucec_PASS_merge.fam`

## 4.1 Score and filter somatic sSNVss

Take file of all scored sSNVs and merge with somatic variant tables to annotate if variant appears in somatic cohort.

```bash
python score_filter_tcga_variants.py \
    -vf ../../data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon_SNVConNeg.tsv \
    -tf [BRCA variant list] \
    -tf [UCEC variant list] \
    -tag BRCA \
    -tag UCEC \
    -ac 0 -m all_syn \
    -o ../../data/4_tcga_analysis/ \
    -t TCGA_ConNegCP3

# Input:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon_SNVConNeg.tsv
#   [BRCA variant list]
#   [UCEC variant list]
# Output:
#   data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_all_syn_TCGA_ConNegCP3.tsv
```

Run again to store just the intersection of the TCGA table with the scored variant table:

```bash
python score_filter_tcga_variants.py \
    -vf ../../data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon_SNVConNeg.tsv \
    -tf [BRCA variant list] \
    -tf [UCEC variant list] \
    -tag BRCA \
    -tag UCEC \
    -ac 0 -m intx \
    -o ../../data/4_tcga_analysis/ \
    -t TCGA_ConNegCP3

# Input:
#   data/3_codon_context_score/syn_variant_12nt_codon_context_score_CP3_zCodon_SNVConNeg.tsv
#   [BRCA variant list]
#   [UCEC variant list]
# Output:
#   data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_intx_TCGA_ConNegCP3.tsv
```

## 4.2 Measure score enrichment

For each cohort and the combined cohort, measure enrichment of high context effect scores
among somatic variants compared to complement of variants.

```bash
python get_enriched_scored_contexts.py \
    -cf ../../data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_all_syn_TCGA_ConNegCP3.tsv \
    -thf ../../data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv \
    -cs "A>C" "A>T" "C>A" "C>G" "C>T" "CpG>CpA" "CpG>TpG" "G>A" "G>C" \
    -thi 5 \
    -g SNVContext \
    -tag BRCA UCEC \
    -s diff_sum_context_score_REF_Codon_zscore \
    -o ../../data/4_tcga_analysis/ \
    -t TCGA_AC0bg_ConNegCP3

# Input:
#   data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_all_syn_TCGA_ConNegCP3.tsv
#   data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv
# Output:
#   data/4_tcga_analysis/scores_enriched_on_5quantile_cohortBRCA_UCEC_TCGA_AC0bg_ConNegCP3.tsv
```

## 4.3 Pathway analysis on enriched contexts

For each cohort and set of enriched contexts, map high context effect variants to genes and run gene set through Over-Representation Analysis.

```bash
python gene_ora_on_contexts.py \
    -cf ../../data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_intx_TCGA_ConNegCP3.tsv \
    -thf ../../data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv \
    -tag BRCA \
    -cs "C>G" "G>A" \
    -thi 5 \
    -g SNVContext \
    -s diff_sum_context_score_REF_Codon_zscore \
    --save-cohort-table \
    -o ../../data/4_tcga_analysis/ \
    -t AC0_5q_CG_GA

# Input:
#   data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_intx_TCGA_ConNegCP3.tsv
#   data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv
# Output:
#   data/4_tcga_analysis/bg_gene_list_cohortBRCA_AC0_5q_CG_GA.tsv
#   data/4_tcga_analysis/gseapy_ora_kegg_go_cohortBRCA_AC0_5q_CG_GA.tsv
#   data/4_tcga_analysis/input_gene_list_cohortBRCA_AC0_5q_CG_GA.tsv
#   data/4_tcga_analysis/scored_variants_flagged_filtered_cohortBRCA_AC0_5q_CG_GA.tsv

python gene_ora_on_contexts.py \
    -cf ../../data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_intx_TCGA_ConNegCP3.tsv \
    -thf ../../data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv \
    -tag UCEC \
    -cs "C>T" "G>A" \
    -thi 5 \
    -g SNVContext \
    -s diff_sum_context_score_REF_Codon_zscore \
    --save-cohort-table \
    -o ../../data/4_tcga_analysis/ \
    -t AC0_5q_CT_GA

# Input
#   data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_intx_TCGA_ConNegCP3.tsv
#   data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv
# Output
#   4_tcga_analysis/bg_gene_list_cohortUCEC_AC0_5q_CT_GA.tsv
#   4_tcga_analysis/gseapy_ora_kegg_go_cohortUCEC_AC0_5q_CT_GA.tsv
#   4_tcga_analysis/input_gene_list_cohortUCEC_AC0_5q_CT_GA.tsv
#   4_tcga_analysis/scored_variants_flagged_filtered_cohortUCEC_AC0_5q_CT_GA.tsv

python gene_ora_on_contexts.py \
    -cf ../../data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_intx_TCGA_ConNegCP3.tsv \
    -thf ../../data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv \
    -tag BRCA UCEC \
    -cs "C>G" "G>A" "C>T" \
    -thi 5 \
    -g SNVContext \
    -s diff_sum_context_score_REF_Codon_zscore \
    --save-cohort-table \
    -o ../../data/4_tcga_analysis/ \
    -t AC0_5q_CG_GA_CT

# Input
#   data/4_tcga_analysis/context_scored_variant_frq_cohortsBRCA_UCEC_AC0_intx_TCGA_ConNegCP3.tsv
#   data/3_codon_context_score/context_score_threshold_quantile_xSNVContext_diff_sum_context_score_posScaled.tsv
# Output
#   4_tcga_analysis/bg_gene_list_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv
#   4_tcga_analysis/gseapy_ora_kegg_go_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv
#   4_tcga_analysis/input_gene_list_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv
#   4_tcga_analysis/scored_variants_flagged_filtered_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv
```

## 4.4 Pathway analysis on random gene sets

Take scored variant table flagged with variants used in the original ORA and subsample
variants from the non-flagged remainder (which should represent a background data set).
Map those variants to the set of genes and use that gene set as input (with same background gene list) for ORA (repeat specified number of times).

```bash
python gene_ora_on_random_gene_lists_byContext.py \
    -bvf ../../data/4_tcga_analysis/scored_variants_flagged_filtered_cohortBRCA_AC0_5q_CG_GA.tsv \
    -bg ../../data/4_tcga_analysis/bg_gene_list_cohortBRCA_AC0_5q_CG_GA.tsv \
    -g SNVContext \
    -nr 100 \
    -tag BRCA \
    -o ../../data/4_tcga_analysis/ \
    -t AC0_5q_CG_GA

# Input:
#   data/4_tcga_analysis/scored_variants_flagged_filtered_cohortBRCA_AC0_5q_CG_GA.tsv
#   data/4_tcga_analysis/bg_gene_list_cohortBRCA_AC0_5q_CG_GA.tsv
# Output:
#   data/4_tcga_analysis/gseapy_sampled_ora_kegg_go_cohortBRCA_n100_AC0_5q_CG_GA.tsv
 
python gene_ora_on_random_gene_lists_byContext.py \
    -bvf ../../data/4_tcga_analysis/scored_variants_flagged_filtered_cohortUCEC_AC0_5q_CT_GA.tsv \
    -bg ../../data/4_tcga_analysis/bg_gene_list_cohortUCEC_AC0_5q_CT_GA.tsv \
    -g SNVContext \
    -nr 100 \
    -tag UCEC \
    -o ../../data/4_tcga_analysis/ \
    -t AC0_5q_CT_GA

#Input:
#   data/4_tcga_analysis/scored_variants_flagged_filtered_cohortUCEC_AC0_5q_CT_GA.tsv
#   data/4_tcga_analysis/bg_gene_list_cohortUCEC_AC0_5q_CT_GA.tsv
#Output:
#   data/4_tcga_analysis/gseapy_sampled_ora_kegg_go_cohortUCEC_n100_AC0_5q_CT_GA.tsv

python gene_ora_on_random_gene_lists_byContext.py \
    -bvf ../../data/4_tcga_analysis/scored_variants_flagged_filtered_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv \
    -bg ../../data/4_tcga_analysis/bg_gene_list_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv \
    -g SNVContext \
    -nr 100 \
    -tag BRCA UCEC \
    -o ../../data/4_tcga_analysis/ \
    -t AC0_5q_CG_GA_CT

# Input:
#   data/4_tcga_analysis/scored_variants_flagged_filtered_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv
#   data/4_tcga_analysis/bg_gene_list_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv
# Output:
#   data/4_tcga_analysis/gseapy_sampled_ora_kegg_go_cohortBRCA_UCEC_n100_AC0_5q_CG_GA_CT.tsv
```

## 4.5 Process and summarize analyses

* Determine mutational signatures and cohorts that are enriched for the specified set of context scores

* Determine gene sets that are significantly enriched among genes that are represented by high context effect variant sets from each cohort

```bash
4_tcga_analysis_processing.ipynb

# Input:
#   data/4_tcga_analysis/scores_enriched_on_5quantile_cohortBRCA_UCEC_TCGA_AC0bg_ConNegCP3.tsv
#   data/4_tcga_analysis/gseapy_ora_kegg_go_cohortBRCA_AC0_5q_CG_GA.tsv
#   data/4_tcga_analysis/gseapy_ora_kegg_go_cohortUCEC_AC0_5q_CT_GA.tsv
#   data/4_tcga_analysis/gseapy_ora_kegg_go_cohortBRCA_UCEC_AC0_5q_CG_GA_CT.tsv
#   data/4_tcga_analysis/gseapy_sampled_ora_kegg_go_cohortBRCA_n100_AC0_5q_CG_GA.tsv
#   data/4_tcga_analysis/gseapy_sampled_ora_kegg_go_cohortUCEC_n100_AC0_5q_CT_GA.tsv
#   data/4_tcga_analysis/gseapy_sampled_ora_kegg_go_cohortBRCA_UCEC_n100_AC0_5q_CG_GA_CT.tsv
# Output:
#   data/4_tcga_analysis/gseapy_ora_summary_kegg_go_cohortBRCA_UCEC_n100_AC0_5q_CT_CG_GA.tsv   
```

