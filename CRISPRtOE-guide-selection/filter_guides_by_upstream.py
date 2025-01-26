#!/usr/bin/env python3

import argparse
import pandas as pd
import os
from Bio import SeqIO
import sys
import numpy as np
from multiprocessing import Pool

def parse_arguments():
    """Set up command line arguments"""
    parser = argparse.ArgumentParser(description="Filter sgRNA guides based on their position relative to genes")
    parser.add_argument("-g", "--guides", required=True, help="Path to input guides TSV file")
    parser.add_argument("-c", "--coordinates", required=True, help="Path to coordinates TSV file")
    parser.add_argument("-b", "--genbank", required=True, help="Path to GenBank file (needed for chromosome sizes)")
    parser.add_argument("-l", "--loci", help="Optional: Path to file containing genes/loci of interest")
    parser.add_argument("-o", "--output", help="Optional: Path for output file. Default: input_name_upstream_filtered.tsv")
    parser.add_argument("-p", "--proximal", type=int, default=102, 
                       help="Distance from gene start for proximal boundary of search region (default: 102)")
    parser.add_argument("-d", "--distal", type=int, default=302, 
                       help="Distance from gene start for distal boundary of search region (default: 302)")
    parser.add_argument("-t", "--threads", type=int, default=None,
                       help="Number of processor threads to use (default: auto)")
    return parser.parse_args()

def get_chromosome_sizes(genbank_path):
    """Extract chromosome sizes from GenBank file"""
    chromosome_sizes = {}
    for record in SeqIO.parse(genbank_path, "genbank"):
        chromosome_sizes[record.id] = len(record.seq)
    return chromosome_sizes

def validate_chromosome_ids(guides_df, coordinates_df, chromosome_sizes):
    """Validate that chromosome IDs in guides and coordinates files exist in GenBank"""
    guide_chroms = set(guides_df['chr'].unique())
    coord_chroms = set(coordinates_df['chr'].unique())
    genbank_chroms = set(chromosome_sizes.keys())
    
    missing_chroms = (guide_chroms | coord_chroms) - genbank_chroms
    
    if missing_chroms:
        print("WARNING: Found chromosome IDs in input files that don't exist in GenBank:")
        for chrom in missing_chroms:
            print(f"  {chrom}")
        print("\nChromosomes in GenBank:")
        for chrom, size in chromosome_sizes.items():
            print(f"  {chrom}: {size:,} bp")
            
    return missing_chroms

def load_input_files(guides_path, coordinates_path, loci_path=None):
    """Load and validate input files"""
    # Check if files exist
    if not os.path.exists(guides_path):
        raise FileNotFoundError(f"Guides file not found: {guides_path}")
    if not os.path.exists(coordinates_path):
        raise FileNotFoundError(f"Coordinates file not found: {coordinates_path}")
    if loci_path and not os.path.exists(loci_path):
        raise FileNotFoundError(f"Loci file not found: {loci_path}")
        
    print("Found all input files")
    
    # Load main input files
    try:
        guides_df = pd.read_csv(guides_path, sep='\t')
        print(f"Loaded guides file with columns: {guides_df.columns.tolist()}")
    except Exception as e:
        raise Exception(f"Error loading guides file: {str(e)}")
        
    try:
        coordinates_df = pd.read_csv(coordinates_path, sep='\t')
        print(f"Loaded coordinates file with columns: {coordinates_df.columns.tolist()}")
    except Exception as e:
        raise Exception(f"Error loading coordinates file: {str(e)}")
    
    # Process optional loci file if provided
    missing_loci = []
    if loci_path:
        loci_df = pd.read_csv(loci_path, sep='\t')
        if 'locus_tag' in loci_df.columns:
            match_col = 'locus_tag'
        elif 'gene' in loci_df.columns:
            match_col = 'gene'
        else:
            raise ValueError("Loci file must have a header 'locus_tag' or 'gene'")
            
        # Get list of loci to find
        loci_to_find = set(loci_df[match_col])
        loci_found = set(coordinates_df[match_col])
        missing_loci = list(loci_to_find - loci_found)
        
        # Filter coordinates to only include requested loci
        coordinates_df = coordinates_df[coordinates_df[match_col].isin(loci_to_find)]
        
        # Report findings
        print(f"\nLoci file summary:")
        print(f"Total loci in input file: {len(loci_to_find)}")
        print(f"Loci found in coordinates: {len(loci_to_find) - len(missing_loci)}")
        print(f"Loci not found: {len(missing_loci)}")
    
    # Make sure we're returning all three values
    return guides_df, coordinates_df, missing_loci
        
def get_true_start_and_upstream_region(start, end, direction, proximal, distal, chromosome_size):
    """Determine true start position and upstream region based on direction"""
    crosses_origin = False
    if end < start:
        crosses_origin = True
        
    if direction == 'F':
        true_start = start
        if crosses_origin and start - distal < 0:
            # Handle upstream region crossing origin for forward features
            upstream_start = chromosome_size - (distal - (chromosome_size - start))
            upstream_end = chromosome_size - (proximal - (chromosome_size - start))
        else:
            upstream_start = start - distal
            upstream_end = start - proximal
    else:  # direction == 'R'
        true_start = end
        upstream_start = true_start + proximal
        upstream_end = true_start + distal
        
    return true_start, upstream_start, upstream_end

def is_position_in_range(pos, start, end, chromosome_size=None):
    """Check if a position falls within a range, handling origin crossing"""
    if start <= end:
        return start <= pos <= end
    else:  # Range crosses origin
        return pos >= start or pos <= end

def process_chromosome(args):
    """Process a single chromosome's worth of guides and genes"""
    chrom_guides, chrom_coords, proximal, distal, chrom_size = args
    filtered_guides = []
    
    # Pre-filter guides by direction
    grouped_guides = dict(tuple(chrom_guides.groupby('sp_dir')))
    
    # Process each gene
    for _, gene in chrom_coords.iterrows():
        # Only look at guides with matching direction
        if gene['gene_dir'] not in grouped_guides:
            continue
            
        matching_guides = grouped_guides[gene['gene_dir']]
        
        # Calculate upstream region for gene
        _, upstream_start, upstream_end = get_true_start_and_upstream_region(
            gene['feature_start'],
            gene['feature_end'],
            gene['gene_dir'],
            proximal,
            distal,
            chrom_size
        )
        
        # Get guide start positions based on direction
        guide_starts = np.where(gene['gene_dir'] == 'R', 
                              matching_guides['tar_end'], 
                              matching_guides['tar_start'])
        
        # Vectorized position check
        if upstream_start <= upstream_end:
            mask = (guide_starts >= upstream_start) & (guide_starts <= upstream_end)
        else:  # Region crosses origin
            mask = (guide_starts >= upstream_start) | (guide_starts <= upstream_end)
            
        # Get matching guides
        matches = matching_guides[mask].copy()
        if not matches.empty:
            matches['upstream_of_locus'] = gene['locus_tag']
            matches['upstream_of_gene'] = gene['gene']
            matches['gene_dir'] = gene['gene_dir']
            matches['gene_start'] = gene['feature_start']
            matches['gene_end'] = gene['feature_end']
            
            # Calculate bp_upstream for each guide individually
            if gene['gene_dir'] == 'F':
                matches['bp_upstream'] = matches.apply(
                    lambda row: (chrom_size - row['tar_start'] + gene['feature_start']) 
                    if row['tar_start'] > gene['feature_start']
                    else (gene['feature_start'] - row['tar_start']),
                    axis=1
                )
            else:  # gene_dir == 'R'
                matches['bp_upstream'] = matches.apply(
                    lambda row: (chrom_size - gene['feature_end'] + row['tar_end'])
                    if row['tar_end'] < gene['feature_end']
                    else (row['tar_end'] - gene['feature_end']),
                    axis=1
                )
            
            filtered_guides.append(matches)
    
    return pd.concat(filtered_guides) if filtered_guides else pd.DataFrame()

def filter_guides_parallel(guides_df, coordinates_df, chromosome_sizes, proximal, distal, n_processes=None):
    """Modified parallel implementation of guide filtering"""
    # Add predicted insertion site
    guides_df['predicted_insertion_site'] = np.where(
        guides_df['sp_dir'] == 'F',
        guides_df['tar_start'] + 82,
        guides_df['tar_end'] - 82
    )

    # Remove unnecessary columns
    columns_to_remove = ['locus_tag', 'gene', 'feature_types', 'offset', 'overlap', 'tar_dir']
    guides_df = guides_df.drop(columns=[col for col in columns_to_remove if col in guides_df.columns])
    
    # Process each chromosome
    chr_args = []
    for chrom in coordinates_df['chr'].unique():
        chrom_coords = coordinates_df[coordinates_df['chr'] == chrom]
        chrom_guides = guides_df[guides_df['chr'] == chrom]
        
        if not chrom_guides.empty and not chrom_coords.empty:
            chr_args.append((
                chrom_guides,
                chrom_coords,
                proximal,
                distal,
                chromosome_sizes[chrom]
            ))
    
    print(f"Starting parallel processing with {len(chr_args)} chromosome arguments")
    # Process chromosomes in parallel
    with Pool(processes=n_processes) as pool:
        results = pool.map(process_chromosome, chr_args)
    
    print("Parallel processing complete")
    # Combine results
    final_df = pd.concat(results) if results else pd.DataFrame()
    print(f"Combined results into dataframe with {len(final_df)} rows")
    return final_df

def report_guides_per_gene(filtered_guides, output_path, coordinates_df, loci_df=None):
    """Generate a summary of how many guides were found per gene"""
    # Start with all genes from coordinates (or subset if loci provided)
    if loci_df is not None:
        # Get the matching column (locus_tag or gene)
        if 'locus_tag' in loci_df.columns:
            match_col = 'locus_tag'
        else:
            match_col = 'gene'
        base_genes = coordinates_df[coordinates_df[match_col].isin(loci_df[match_col])]
    else:
        base_genes = coordinates_df
    
    # Create base dataframe with all genes
    summary_df = pd.DataFrame({
        'upstream_of_locus': base_genes['locus_tag'],
        'upstream_of_gene': base_genes['gene']
    })
    
    # Count guides per gene
    if not filtered_guides.empty:
        guide_counts = filtered_guides.groupby(['upstream_of_locus', 'upstream_of_gene']).size().reset_index()
        guide_counts.columns = ['upstream_of_locus', 'upstream_of_gene', 'guide_count']
        
        # Merge with base genes
        summary_df = summary_df.merge(guide_counts, 
                                    on=['upstream_of_locus', 'upstream_of_gene'], 
                                    how='left')
    else:
        summary_df['guide_count'] = 0
    
    # Fill NaN guide counts with 0
    summary_df['guide_count'] = summary_df['guide_count'].fillna(0).astype(int)
    
    # Sort by guide count (descending) and locus tag
    summary_df = summary_df.sort_values(['guide_count', 'upstream_of_locus'], 
                                      ascending=[False, True])
    
    # Print summary statistics
    total_genes = len(summary_df)
    genes_with_guides = (summary_df['guide_count'] > 0).sum()
    total_guides = summary_df['guide_count'].sum()
    
    print("\nGuides found per gene:")
    print(f"Total genes analyzed: {total_genes}")
    print(f"Genes with at least one guide: {genes_with_guides}")
    print(f"Genes with no guides: {total_genes - genes_with_guides}")
    print(f"Total guides found: {total_guides}")
    if genes_with_guides > 0:
        print(f"Average guides per gene (among genes with guides): {total_guides/genes_with_guides:.1f}")
    
    print("\nTop 10 genes by guide count:")
    print(summary_df.head(10).to_string(index=False))
    
    # Save the summary to a file
    summary_path = os.path.splitext(output_path)[0] + "_summary.tsv"
    summary_df.to_csv(summary_path, sep='\t', index=False)
    print(f"\nFull summary saved to: {summary_path}")
    
def write_outputs(filtered_guides, output_path, coordinates_df, missing_loci=None):
    """Write filtered guides and missing loci report"""
    # Write main output
    filtered_guides.to_csv(output_path, sep='\t', index=False)
    print(f"\nWrote {len(filtered_guides)} filtered guides to {output_path}")
    
    # Generate and write guide count summary
    report_guides_per_gene(filtered_guides, output_path, coordinates_df)
    
    # Write missing loci report if applicable
    if missing_loci:
        missing_loci_path = os.path.splitext(output_path)[0] + "_missing_loci.txt"
        with open(missing_loci_path, 'w') as f:
            f.write("# Loci not found in coordinates file:\n")
            for locus in missing_loci:
                f.write(f"{locus}\n")
        print(f"Wrote list of {len(missing_loci)} missing loci to {missing_loci_path}")
        

def main():
    try:
        args = parse_arguments()
        print("Arguments parsed successfully:")
        print(f"  Guides file: {args.guides}")
        print(f"  Coordinates file: {args.coordinates}")
        print(f"  GenBank file: {args.genbank}")
        print(f"  Proximal distance: {args.proximal}")
        print(f"  Distal distance: {args.distal}")
        if args.loci:
            print(f"  Loci file: {args.loci}")
        
        # Get chromosome sizes from GenBank
        print("Reading chromosome sizes from GenBank...")
        chromosome_sizes = get_chromosome_sizes(args.genbank)
        print(f"Found {len(chromosome_sizes)} chromosomes:")
        for chrom, size in chromosome_sizes.items():
            print(f"  {chrom}: {size:,} bp")
        
        # Load input files
        print("Loading input files...")
        guides_df, coordinates_df, missing_loci = load_input_files(
            args.guides, args.coordinates, args.loci
        )
        
        # Validate chromosome IDs
        missing_chroms = validate_chromosome_ids(guides_df, coordinates_df, chromosome_sizes)
        if missing_chroms:
            print("WARNING: Proceeding with analysis, but some chromosomes may be skipped.")
        
        # Filter guides based on upstream regions
        print("Starting guide filtering...")
        filtered_guides = filter_guides_parallel(
            guides_df, coordinates_df, chromosome_sizes,
            args.proximal, args.distal,
            n_processes=args.threads if hasattr(args, 'threads') else None
        )
        
        # Set output path
        if not args.output:
            base = os.path.splitext(args.guides)[0]
            output_path = f"{base}_upstream_filtered.tsv"
        else:
            output_path = args.output
            
        write_outputs(filtered_guides, output_path, coordinates_df, missing_loci)
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        import traceback
        print("\nFull traceback:")
        print(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()