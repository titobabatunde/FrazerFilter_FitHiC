#!/usr/bin/env python3
"""
Filter individual Fit-Hi-C loop replicate files by FDR threshold and neighboring bin criteria.

This script:
1. Filters each replicate file by FDR threshold (default 0.0001)
2. Filters each replicate by requiring significant interactions with neighbors around opposing anchors (Frazer filtering)
3. Outputs filtered files for each replicate independently
"""

import pandas as pd
import numpy as np
import argparse
import os
from pathlib import Path
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Filter individual fithic loop replicate files by FDR and neighbor significance')
    parser.add_argument('--input_file', 
                       required=True,
                       help='Input fithic significance file (.txt)')
    parser.add_argument('--replicate_name', 
                       required=True,
                       help='Replicate name (e.g., Th2-1)')
    parser.add_argument('--chromosome', 
                       type=str,
                       default='chr1',
                       help='Chromosome (default: chr1)')
    parser.add_argument('--fdr_threshold', 
                       type=float,
                       default=0.0001,
                       help='FDR threshold for significance (default: 0.0001)')
    parser.add_argument('--min_neighbors', 
                       type=int,
                       default=3,
                       help='Minimum number of significant neighbors required (default: 3)')
    parser.add_argument('--total_neighbors', 
                       type=int,
                       default=5,
                       help='Total number of neighbors to check (default: 5, excluding anchor)')
    parser.add_argument('--resolution', 
                       type=int,
                       required=True,
                       help='Resolution in bp (e.g., 10000)')

    parser.add_argument('--output_dir', 
                       required=True,
                       type=str,
                       help='Output directory')
    parser.add_argument('--verbose', 
                       action='store_true',
                       help='Print verbose output')
    return parser.parse_args()

def load_and_sort_fithic_data(file_path, fdr_threshold, chromosome, verbose=False):
    """
    Load and filter fithic data by FDR threshold and distance.
    
    Args:
        file_path: Path to fithic significance file (.txt or .txt.gz)
        fdr_threshold: FDR threshold for significance filtering
        chromosome: Chromosome to filter for (e.g., 'chr1')
        verbose: Whether to print verbose output
    
    Returns:
        df: Full dataframe sorted by coordinates
        significant_df: Filtered dataframe with significant interactions (FDR < threshold, distance >= 32000 bp)
    """
    if verbose:
        print(f"Loading: {file_path}")
    
    # Handle gzipped files
    compression = 'gzip' if file_path.endswith('.gz') else None
    
    # Load the data
    df = pd.read_csv(file_path, sep='\t', compression=compression)
    
    if verbose:
        print(f"  Initial: {len(df)} interactions")
        print(f"  Columns: {list(df.columns)}")

    # Filter for specified chromosome (both anchors must be on the same chromosome)
    df = df[df['chr1'] == chromosome]
    df = df[df['chr2'] == chromosome]

    # Sort by coordinates for efficient neighbor searching
    df = df.sort_values(['chr1', 'fragmentMid1', 'chr2', 'fragmentMid2']).reset_index(drop=True).copy()
    
    # Filter by FDR threshold
    significant_df = df[df['q-value'] < fdr_threshold].copy()
    if verbose:
        print(f"  After FDR filtering: {len(significant_df)} interactions")
    
    # Filter out interactions with absolute distance between anchors less than 32000 bp
    # This removes very short-range interactions that are likely technical artifacts
    significant_df = significant_df[abs(significant_df['fragmentMid1'] - significant_df['fragmentMid2']) >= 32000]
    if verbose:
        print(f"  After distance filtering (>= 32000 bp): {len(significant_df)} interactions")
    
    return df.copy(), significant_df.copy()
# end def

def filter_by_neighbor_significance(significant_df, df_sorted, resolution, qvalue_threshold, min_neighbors, max_neighbors, verbose=False):
    """
    Filter significant interactions by checking if neighbors around each opposing anchor are significant.
    
    This function implements the Frazer filtering approach: for an interaction A ↔ B to be kept,
    both anchors must have significant interactions with neighbors around the opposing anchor.
    
    For interaction A ↔ B:
    - Anchor A: Check if anchor A has ≥min_neighbors significant interactions with upstream/downstream 
      bins around anchor B (i.e., A ↔ B-1, A ↔ B-2, ..., A ↔ B+1, A ↔ B+2, ...)
    - Anchor B: Check if anchor B has ≥min_neighbors significant interactions with upstream/downstream 
      bins around anchor A (i.e., A-1 ↔ B, A-2 ↔ B, ..., A+1 ↔ B, A+2 ↔ B, ...)
    
    An interaction passes if BOTH anchors pass the neighbor significance test.
    
    Args:
        significant_df: DataFrame with significant interactions (already filtered by FDR)
        df_sorted: Full sorted DataFrame (for neighbor lookup)
        resolution: Resolution in base pairs
        qvalue_threshold: FDR threshold for significance
        min_neighbors: Minimum number of significant neighbors required
        max_neighbors: Maximum number of neighbors to check in each direction
        verbose: Whether to print verbose output
    
    Returns:
        filtered_df: DataFrame with interactions that pass the neighbor significance test
    """
    if verbose:
        print(f"Filtering significant interactions by checking if the anchors have significant interactions with neighbors around the other anchor")
        print(f"Only keep interactions where both anchors pass the neighbor significance test")

        print(f" Head of significant_df:\n {significant_df.head()}")
        print(f" Head of df_sorted:\n {df_sorted.head()}")
    
    # For each significant interaction, check if the anchor has significant interactions with neighbors around the other anchor
    significant_neighbors = []
    if len(significant_df) == 0:
        if verbose:
            print("No significant interactions to check neighbors for")
        return significant_df
        
    if verbose:
        print(f"Checking neighbors for {len(significant_df)} significant interactions")
        
    for _, row in significant_df.iterrows():
        chrA, midA = row['chr1'], row['fragmentMid1']
        chrB, midB = row['chr2'], row['fragmentMid2']

        # Check anchor A: how many upstream/downstream bins from anchor B are significant?
        upstream = [midB - i * resolution for i in range(1, max_neighbors + 1)]
        downstream = [midB + i * resolution for i in range(1, max_neighbors + 1)]
        
        # Count significant upstream neighbors for anchor A (A ↔ B-1, A ↔ B-2, ...)
        upstream_df_A = df_sorted[(df_sorted['chr1'] == chrA) & (df_sorted['fragmentMid2'].isin(upstream))]
        upstream_count_A = (upstream_df_A['q-value'].astype(float) < qvalue_threshold).sum()
        
        # Count significant downstream neighbors for anchor A (A ↔ B+1, A ↔ B+2, ...)
        downstream_df_A = df_sorted[(df_sorted['chr1'] == chrA) & (df_sorted['fragmentMid2'].isin(downstream))]
        downstream_count_A = (downstream_df_A['q-value'].astype(float) < qvalue_threshold).sum()

        # Anchor A passes if it has enough significant neighbors in either direction
        anchorA_pass = upstream_count_A >= min_neighbors or downstream_count_A >= min_neighbors

        # Check anchor B: how many upstream/downstream bins from anchor A are significant?
        # We check if anchor B has significant interactions with neighbors around anchor A
        upstream_A = [midA - i * resolution for i in range(1, max_neighbors + 1)]
        downstream_A = [midA + i * resolution for i in range(1, max_neighbors + 1)]

        # Count significant upstream neighbors for anchor B (A-1 ↔ B, A-2 ↔ B, ...)
        upstream_df_B = df_sorted[(df_sorted['chr2'] == chrB) & (df_sorted['fragmentMid1'].isin(upstream_A))]
        upstream_count_B = (upstream_df_B['q-value'].astype(float) < qvalue_threshold).sum()
        
        # Count significant downstream neighbors for anchor B (A+1 ↔ B, A+2 ↔ B, ...)
        downstream_df_B = df_sorted[(df_sorted['chr2'] == chrB) & (df_sorted['fragmentMid1'].isin(downstream_A))]
        downstream_count_B = (downstream_df_B['q-value'].astype(float) < qvalue_threshold).sum()

        # Anchor B passes if it has enough significant neighbors in either direction
        anchorB_pass = upstream_count_B >= min_neighbors or downstream_count_B >= min_neighbors

        # Keep interaction only if BOTH anchors pass the neighbor significance test
        passes_filter = anchorA_pass and anchorB_pass
        significant_neighbors.append(passes_filter)
        
        if verbose:
            total_neighbors = max_neighbors # upstream or downstream
            print(f"  {chrA}:{midA} ↔ {chrB}:{midB}: A={upstream_count_A}/{total_neighbors} B={upstream_count_B}/{total_neighbors} → {'PASS' if passes_filter else 'FAIL'}")
            print(f"  {chrA}:{midA} ↔ {chrB}:{midB}: A={downstream_count_A}/{total_neighbors} B={downstream_count_B}/{total_neighbors} → {'PASS' if passes_filter else 'FAIL'}")
    # end for
    filtered_df = significant_df[significant_neighbors]
    
    if verbose:
        if len(significant_df) > 0:
            print(f"Filtered {len(significant_df)} → {len(filtered_df)} interactions ({len(filtered_df)/len(significant_df)*100:.1f}% retained)")
        else:
            print(f"No significant interactions found to filter")
    
    return filtered_df
# end def


def fithic_df_to_bedpe(fithic_df, resolution):
    """
    Convert fithic dataframe to BEDPE format.
    
    BEDPE format columns: chr1, start1, end1, chr2, start2, end2, name, score
    - start/end coordinates are calculated from fragmentMid ± resolution/2
    - name: loop_{row_number}
    - score: -log10(q-value) for visualization
    
    Args:
        fithic_df: DataFrame with fithic results (must contain fragmentMid1, fragmentMid2, q-value)
        resolution: Resolution in base pairs
    
    Returns:
        bedpe_df: DataFrame in BEDPE format
    """
    # Select required columns
    fithic_df = fithic_df[['chr1', 'fragmentMid1', 'chr2', 'fragmentMid2', 'q-value']].copy()
    
    # Calculate start and end coordinates from fragment midpoints
    fithic_df['start1'] = fithic_df['fragmentMid1'] - resolution//2
    fithic_df['end1'] = fithic_df['fragmentMid1'] + resolution//2
    fithic_df['start2'] = fithic_df['fragmentMid2'] - resolution//2
    fithic_df['end2'] = fithic_df['fragmentMid2'] + resolution//2
    
    # Add loop names and scores
    fithic_df['name'] = [f'loop_{i}' for i in range(len(fithic_df))]
    fithic_df['score'] = -np.log10(fithic_df['q-value'])

    # Return in BEDPE format
    fithic_df = fithic_df[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score']]
    return fithic_df
# end def

def main():
    args = parse_arguments()
    
    if args.verbose:
        print("=== Fit-Hi-C Replicate Filtering ===")
        print(f"Input file: {args.input_file}")
        print(f"Replicate name: {args.replicate_name}")
        print(f"Chromosome: {args.chromosome}")
        print(f"FDR threshold: {args.fdr_threshold}")
        print(f"Resolution: {args.resolution}")
        print(f"Min neighbors: {args.min_neighbors}")
        print(f"Total neighbors: {args.total_neighbors}")
        print()
    # end if

    # Validate input file (must be .txt or .txt.gz)
    if not (args.input_file.endswith('.txt') or args.input_file.endswith('.txt.gz')):
        print(f"Error: Input file {args.input_file} must have a .txt or .txt.gz extension")
        exit(1)
    # end if
    
    # Step 1: Load fithic data
    if args.verbose:
        print("Step 1: Loading fithic data...")
    
    # Create output directory
    output_path = f"{args.output_dir}/frazer-chrs"
    os.makedirs(output_path, exist_ok=True)

    # Step 1: Load data, sort it, and identify significant interactions
    fithic_df = {}
    fithic_df_significant = {}
    extension = f'{args.chromosome}.frazer.fdr{args.fdr_threshold}.txt'
    fithic_df, fithic_df_significant = load_and_sort_fithic_data(
        args.input_file, args.fdr_threshold, args.chromosome, args.verbose)

    # Step 2: Filter by neighbor significance (Frazer filtering)
    if args.verbose:
        print("Step 2: Filtering by neighbor significance (Frazer filtering)...")
    

    if args.verbose:
        print(f"  Processing {args.replicate_name}...")
    fithic_df_significant = filter_by_neighbor_significance(fithic_df_significant, 
                                                            fithic_df, 
                                                            args.resolution, 
                                                            args.fdr_threshold, 
                                                            args.min_neighbors, 
                                                            args.total_neighbors, 
                                                            args.verbose)

    # Step 3: Output filtered files for each replicate independently
    if args.verbose:
        print("Step 3: Creating per-replicate output files...")
    
    total_loops = 0

    output_fithic_file = os.path.join(output_path, f'{args.replicate_name}.{extension}')
    fithic_df_significant.to_csv(output_fithic_file, sep='\t', index=False)

    output_bedpe_file = os.path.join(output_path, f'{args.replicate_name}.{extension.replace("txt", "bedpe")}')
    fithic_df_significant_bedpe = fithic_df_to_bedpe(fithic_df_significant, args.resolution)
    fithic_df_significant_bedpe.to_csv(output_bedpe_file, sep='\t', index=False)
    
    total_loops = len(fithic_df_significant)

    # Print summary
    print(f"\n=== Filtering Summary ===")
    print(f"Input file: {args.input_file}")
    print(f"FDR threshold: {args.fdr_threshold}")
    print(f"Total filtered loops: {total_loops}")
    print()
    


if __name__ == "__main__":
    main()

