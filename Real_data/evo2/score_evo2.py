#!/usr/bin/env python3

# Required imports
import os
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from evo2.models import Evo2

# --- Configuration ---
WINDOW_SIZE = 8192
# Keep FASTA_PATH absolute so it can be found from any trait directory
FASTA_PATH = '/scratch/gilbreth/chen4422/hapnest/to_depot_locus_to_evo2/hg38.ml.fa'

print("Loading hg38 reference genome index...")
genome_index = SeqIO.index(FASTA_PATH, "fasta")

def parse_sequences(chrom, pos, ref, alt):
    """
    Parse reference and variant sequences from the reference genome sequence.
    Returns: ref_seq, var_seq, and a multiplier (1 or -1) to fix allele swaps.
    """
    chrom_str = str(chrom)
    if chrom_str not in genome_index and f"chr{chrom_str}" in genome_index:
        chrom_str = f"chr{chrom_str}"
        
    full_seq = str(genome_index[chrom_str].seq).upper()
    p = int(pos) - 1 # Convert to 0-indexed position

    ref_seq_start = max(0, p - WINDOW_SIZE//2)
    ref_seq_end = min(len(full_seq), p + WINDOW_SIZE//2)
    ref_seq = full_seq[ref_seq_start:ref_seq_end]
    
    snv_pos_in_ref = min(WINDOW_SIZE//2, p)
    genome_allele = ref_seq[snv_pos_in_ref]
    
    # ---------------------------------------------------------
    # DYNAMIC MULTIPLIER LOGIC
    # ---------------------------------------------------------
    if genome_allele == ref:
        # Standard case: Genome matches REF. No score flipping needed.
        var_seq = ref_seq[:snv_pos_in_ref] + alt + ref_seq[snv_pos_in_ref+1:]
        multiplier = 1
        
    elif genome_allele == alt:
        # Swapped case: Genome matches ALT. Inject REF and mark for flipping.
        var_seq = ref_seq[:snv_pos_in_ref] + ref + ref_seq[snv_pos_in_ref+1:]
        multiplier = -1
        
    else:
        # True mismatch
        raise AssertionError(
            f"Severe mismatch at chr{chrom}:{pos}. Expected {ref} or {alt}, but genome has {genome_allele}."
        )

    # Sanity check
    assert len(var_seq) == len(ref_seq), "Sequence lengths do not match."

    return ref_seq, var_seq, multiplier

def main():
    # Set up argument parsing for SLURM compatibility
    parser = argparse.ArgumentParser(description="Score GWAS summary stats with Evo2")
    parser.add_argument('--input', required=True, help="Input GWAS summary stats file (e.g., eur_loci_1_summ)")
    parser.add_argument('--output', required=True, help="Output CSV filename (e.g., loci_1_summ_scored.csv)")
    args = parser.parse_args()

    print(f"\n--- Starting processing for {args.input} ---")
    
    # Read the data
    gwas_df = pd.read_csv(args.input)
    
    # Data structures for tracking
    ref_seqs = []
    ref_seq_to_index = {}
    ref_seq_indexes = []
    var_seqs = []
    
    valid_indices = [] 
    failed_rows = []   
    multipliers = [] # Keep track of the multiplier for every valid row

    # Parse sequences
    print("Parsing sequences and checking alleles against hg38...")
    for idx, row in gwas_df.iterrows():
        try:
            # Unpacking 3 variables from the function
            ref_seq, var_seq, mult = parse_sequences(row['CHR'], row['POS'], row['REF'], row['ALT'])

            if ref_seq not in ref_seq_to_index:
                ref_seq_to_index[ref_seq] = len(ref_seqs)
                ref_seqs.append(ref_seq)
            
            ref_seq_indexes.append(ref_seq_to_index[ref_seq])
            var_seqs.append(var_seq)
            valid_indices.append(idx)
            multipliers.append(mult) # Store the multiplier (1 or -1)
            
        except AssertionError as e:
            failed_rows.append({'index': idx, 'error': str(e)})

    ref_seq_indexes = np.array(ref_seq_indexes)
    multipliers_array = np.array(multipliers)

    # Log parsing results
    print(f"Successfully parsed: {len(valid_indices)} variants.")
    swap_count = np.sum(multipliers_array == -1)
    print(f" -> Standard matches: {len(valid_indices) - swap_count}")
    print(f" -> Allele swaps detected: {swap_count}")
    
    if failed_rows:
        print(f"Failed to parse: {len(failed_rows)} variants (saved to error log).")

    # Load Model
    print("\nLoading Evo2 7B model...")
    model = Evo2('evo2_7b')

    # Score Sequences
    print(f'Scoring likelihoods of {len(ref_seqs)} reference sequences...')
    ref_scores = model.score_sequences(ref_seqs)

    print(f'Scoring likelihoods of {len(var_seqs)} variant sequences...')
    var_scores = model.score_sequences(var_seqs)

    # Calculate Deltas
    print("\nCalculating final corrected delta scores...")
    raw_delta_scores = np.array(var_scores) - np.array(ref_scores)[ref_seq_indexes]

    # DYNAMIC CORRECTION: Multiply the raw scores by the tracked array of 1s and -1s
    corrected_delta_scores = raw_delta_scores * multipliers_array

    # Clean the dataframe (drop rows that failed to parse so lengths match)
    gwas_df_clean = gwas_df.loc[valid_indices].copy()
    gwas_df_clean['evo2_delta_score'] = corrected_delta_scores

    # Save to disk
    gwas_df_clean.to_csv(args.output, index=False)
    print(f"Success! Data saved to: {os.path.abspath(args.output)}\n")

if __name__ == "__main__":
    main()
