#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os
import numpy as np

def parse_arguments():
    """Set up command line arguments"""
    parser = argparse.ArgumentParser(description="Filter sgRNA guides based on feature-free regions")
    parser.add_argument("-g", "--guides", required=True, help="Path to input guides TSV file")
    parser.add_argument("-f", "--feature-free", required=True, help="Path to feature-free regions TSV file")
    parser.add_argument("-o", "--output", help="Optional: Path for output file. Default: input_name_safe_filtered.tsv")
    return parser.parse_args()

def filter_guides_by_safe_regions(guides_df, feature_free_df):
    """Filter guides based on whether their insertion sites fall in feature-free regions"""
    print("\nFiltering guides based on feature-free regions...")
    
    # Add predicted insertion site column
    guides_df['predicted_insertion_site'] = np.where(
        guides_df['sp_dir'] == 'F',
        guides_df['tar_start'] + 82,
        guides_df['tar_end'] - 82
    )
    
    # Sort both dataframes by chromosome and position
    guides_df = guides_df.sort_values(['chr', 'predicted_insertion_site'])
    feature_free_df = feature_free_df.sort_values(['chr', 'start'])
    
    safe_guides = []
    for chrom in guides_df['chr'].unique():
        print(f"Processing chromosome {chrom}...")
        chrom_guides = guides_df[guides_df['chr'] == chrom]
        chrom_regions = feature_free_df[feature_free_df['chr'] == chrom]
        
        if chrom_regions.empty:
            print(f"No feature-free regions found for chromosome {chrom}")
            continue
        
        # Convert to numpy arrays for faster comparison
        insertion_sites = chrom_guides['predicted_insertion_site'].values
        region_starts = chrom_regions['start'].values
        region_ends = chrom_regions['end'].values
        
        # Process guides in sorted order
        current_region_idx = 0
        safe_mask = np.zeros(len(chrom_guides), dtype=bool)
        
        for i, insertion_site in enumerate(insertion_sites):
            # Move to next region if we're past current one
            while (current_region_idx < len(region_starts) and 
                   region_ends[current_region_idx] < insertion_site):
                current_region_idx += 1
            
            # Check if current region contains insertion site
            if current_region_idx < len(region_starts):
                if region_starts[current_region_idx] <= insertion_site <= region_ends[current_region_idx]:
                    safe_mask[i] = True
        
        safe_guides.append(chrom_guides[safe_mask])
        print(f"Found {safe_mask.sum():,} safe guides out of {len(chrom_guides):,}")
    
    result_df = pd.concat(safe_guides) if safe_guides else pd.DataFrame()
    print(f"\nTotal guides with safe insertion sites: {len(result_df):,}")
    return result_df

def main():
    try:
        args = parse_arguments()
        print("Arguments parsed successfully:")
        print(f"  Guides file: {args.guides}")
        print(f"  Feature-free regions file: {args.feature_free}")
        
        # Load input files
        print("\nLoading input files...")
        guides_df = pd.read_csv(args.guides, sep='\t')
        feature_free_df = pd.read_csv(args.feature_free, sep='\t')
        
        print(f"Loaded {len(guides_df):,} guides")
        print(f"Loaded {len(feature_free_df):,} feature-free regions")
        
        # Filter guides
        filtered_guides = filter_guides_by_safe_regions(guides_df, feature_free_df)
        
        # Set output path
        if not args.output:
            base = os.path.splitext(args.guides)[0]
            output_path = f"{base}_safe_filtered.tsv"
        else:
            output_path = args.output
        
        # Write output
        filtered_guides.to_csv(output_path, sep='\t', index=False)
        print(f"\nWrote {len(filtered_guides):,} filtered guides to {output_path}")
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        import traceback
        print("\nFull traceback:")
        print(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()