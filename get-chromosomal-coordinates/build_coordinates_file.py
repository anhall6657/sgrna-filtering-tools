#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import sys
import os
from Bio import SeqIO
from Bio import SeqFeature

def parse_arguments():
    parser = argparse.ArgumentParser(description='Build a coordinates file from guide TSV and GenBank files')
    parser.add_argument('guide_tsv', help='Path to the guide TSV file')
    parser.add_argument('genbank_file', help='Path to the GenBank file')
    parser.add_argument('--output', '-o', help='Output file name (default: GENBANKFILENAME_coordinates.tsv)',
                       default=None)
    return parser.parse_args()

def get_output_filename(genbank_file):
    """Generate output filename based on input GenBank filename."""
    base = os.path.splitext(os.path.basename(genbank_file))[0]
    return f"{base}_coordinates.tsv"

def process_guide_tsv(tsv_file):
    """Extract and process required columns from the guide TSV file."""
    try:
        # First try with default engine
        df = pd.read_csv(tsv_file, sep='\t')
    except:
        # If that fails, try with python engine
        df = pd.read_csv(tsv_file, sep='\t', engine='python')
    
    # Extract required columns
    required_cols = ['chr', 'locus_tag', 'gene', 'tar_dir']
    if not all(col in df.columns for col in required_cols):
        missing = [col for col in required_cols if col not in df.columns]
        sys.exit(f"Error: Missing required columns in guide TSV: {missing}")
    
    coords_df = df[required_cols].copy()
    
    # Remove rows with empty/NA values in locus_tag or gene
    coords_df = coords_df.dropna(subset=['locus_tag', 'gene'])
    coords_df = coords_df[coords_df['locus_tag'].str.strip() != '']
    coords_df = coords_df[coords_df['gene'].str.strip() != '']
    coords_df = coords_df[coords_df['locus_tag'].str.lower() != 'none']
    coords_df = coords_df[coords_df['gene'].str.lower() != 'none']
    
    # De-duplicate
    coords_df = coords_df.drop_duplicates()
    
    # Rename tar_dir column
    coords_df = coords_df.rename(columns={'tar_dir': 'gene_dir'})
    
    return coords_df

def parse_location(location_str):
    """Parse location string to get start and end coordinates."""
    if 'join' in location_str:
        # For features crossing the origin
        numbers = [int(x) for x in re.findall(r'\d+', location_str)]
        start = numbers[0]  # First number is the start
        end = numbers[-1]   # Last number is the end
    else:
        # For simple locations
        numbers = [int(x) for x in re.findall(r'\d+', location_str)]
        start, end = numbers[0], numbers[1]
    return start, end

def extract_coordinates(genbank_file, coords_df):
    """Extract feature coordinates from GenBank file using BioPython."""
    # Initialize new columns
    coords_df['feature_start'] = None
    coords_df['feature_end'] = None
    coords_df['feature_type'] = None
    
    # Create dictionary to store locus tag information
    locus_info = {}
    
    # Track chromosomes processed
    chr_count = 0
    
    # Parse GenBank file
    for record in SeqIO.parse(genbank_file, "genbank"):
        chr_count += 1
        print(f"Processing chromosome/record: {record.id}")
        
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers['locus_tag'][0]
                location_str = str(feature.location)
                start, end = parse_location(location_str)
                locus_info[locus_tag] = (start, end, feature.type)
    
    print(f"Total chromosomes/records processed: {chr_count}")
    
    # Update dataframe with coordinates and feature type
    for idx, row in coords_df.iterrows():
        locus_tag = row['locus_tag']
        if locus_tag in locus_info:
            coords_df.at[idx, 'feature_start'] = locus_info[locus_tag][0]
            coords_df.at[idx, 'feature_end'] = locus_info[locus_tag][1]
            coords_df.at[idx, 'feature_type'] = locus_info[locus_tag][2]
    
    # Print stats about matches
    total_loci = len(coords_df)
    matched_loci = coords_df['feature_start'].notna().sum()
    print(f"\nMatched {matched_loci} out of {total_loci} locus tags from the guide TSV")
    
    # Remove rows where coordinates weren't found
    coords_df = coords_df.dropna(subset=['feature_start', 'feature_end'])
    
    return coords_df

def main():
    args = parse_arguments()
    
    # Process guide TSV
    print("Processing guide TSV file...")
    coords_df = process_guide_tsv(args.guide_tsv)
    
    # Extract coordinates from GenBank
    print("\nExtracting coordinates from GenBank file...")
    coords_df = extract_coordinates(args.genbank_file, coords_df)
    
    # Determine output filename
    if args.output:
        output_file = args.output
    else:
        output_file = get_output_filename(args.genbank_file)
    
    # Add 1 to start coordinates to match GenBank format
    coords_df['feature_start'] = coords_df['feature_start'] + 1
    
    # Save to file
    coords_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nCoordinates file saved as: {output_file}")
    print(f"Number of features processed: {len(coords_df)}")
    print("\nFeature type breakdown:")
    print(coords_df['feature_type'].value_counts())
    
    # Print chromosome breakdown
    print("\nChromosome breakdown:")
    print(coords_df['chr'].value_counts())

if __name__ == "__main__":
    main()