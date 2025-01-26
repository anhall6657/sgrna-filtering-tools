#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
import sys
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Identify feature-free regions from GenBank files")
    parser.add_argument("-i", "--input", required=True, help="Path to input GenBank file")
    parser.add_argument("-o", "--output", help="Path for output files. Default: input_name_safe_regions")
    parser.add_argument("-u", "--upstream_buffer", type=int, default=0,
                       help="Buffer distance (bp) to add upstream of features (default: 0)")
    parser.add_argument("-d", "--downstream_buffer", type=int, default=0,
                       help="Buffer distance (bp) to add downstream of features (default: 0)")
    parser.add_argument("--include-types", help="Comma-separated list of feature types to include. Default: all except 'source'")
    parser.add_argument("--exclude-types", help="Comma-separated list of feature types to exclude. Default: source")
    return parser.parse_args()
    
def get_features_from_genbank(genbank_path, include_types=None, exclude_types=None):
    """Extract features from GenBank file, splitting origin-crossing features"""
    features_list = []
    
    for record in SeqIO.parse(genbank_path, "genbank"):
        chrom = record.id
        chrom_size = len(record.seq)
        
        # Process each feature
        for feature in record.features:
            print(f"\nFeature type: {feature.type}")
            print(f"Location: {feature.location}")
            
            # Skip excluded features
            if feature.type in exclude_types:
                continue
            if include_types and feature.type not in include_types:
                continue
                
            # Handle features that cross origin
            if isinstance(feature.location, CompoundLocation):
                print("Found origin-crossing feature:")
                print(f"Original location: {feature.location}")
                # Get the parts of the compound location
                for part in feature.location.parts:
                    start = part.start + 1  # Convert to 1-based
                    end = part.end
                    print(f"Processing part: {start}..{end}")
                    strand = '+' if part.strand == 1 else '-'
                    
                    features_list.append({
                        'chr': chrom,
                        'feature_type': feature.type,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'locus_tag': feature.qualifiers.get('locus_tag', [''])[0] if 'locus_tag' in feature.qualifiers else '',
                        'gene': feature.qualifiers.get('gene', [''])[0] if 'gene' in feature.qualifiers else '',
                        'chrom_size': chrom_size
                    })
            else:
                # Handle normal features
                start = feature.location.start + 1  # Convert to 1-based
                end = feature.location.end
                strand = '+' if feature.location.strand == 1 else '-'
                
                features_list.append({
                    'chr': chrom,
                    'feature_type': feature.type,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'locus_tag': feature.qualifiers.get('locus_tag', [''])[0] if 'locus_tag' in feature.qualifiers else '',
                    'gene': feature.qualifiers.get('gene', [''])[0] if 'gene' in feature.qualifiers else '',
                    'chrom_size': chrom_size
                })
    
    return pd.DataFrame(features_list)
    
def find_safe_regions(features_df, upstream_buffer=0, downstream_buffer=0):
    """Identify regions without any features, using directional buffers"""
    safe_regions = []
    feature_stats = {}
    
    for chrom in features_df['chr'].unique():
        print(f"\nProcessing chromosome {chrom}")
        chrom_features = features_df[features_df['chr'] == chrom].copy()
        chrom_size = chrom_features['chrom_size'].iloc[0]
        print(f"Chromosome size: {chrom_size}")
        
        # Count features by type for this chromosome
        type_counts = chrom_features['feature_type'].value_counts()
        feature_stats[chrom] = type_counts.to_dict()
        
        # Get all feature regions including directional buffers
        regions = []
        for _, feature in chrom_features.iterrows():
            if feature['strand'] == '+':
                start = max(1, feature['start'] - upstream_buffer)
                end = min(chrom_size, feature['end'] + downstream_buffer)
            else:
                start = max(1, feature['start'] - downstream_buffer)
                end = min(chrom_size, feature['end'] + upstream_buffer)
            regions.append((start, end))
            
        print(f"Created {len(regions)} regions from features")
        
        if not regions:
            safe_regions.append({
                'chr': chrom,
                'start': 1,
                'end': chrom_size,
                'size': chrom_size
            })
            continue
            
        # Sort regions
        regions.sort()
        print(f"First few regions: {regions[:5]}")
        
        # Merge overlapping regions
        merged = []
        current_start, current_end = regions[0]
        
        for start, end in regions[1:]:
            if start <= current_end + 1:
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        merged.append((current_start, current_end))
        
        print(f"Merged into {len(merged)} non-overlapping regions")
        print(f"First few merged regions: {merged[:5]}")
        
        # Find gaps (safe regions)
        if merged[0][0] > 1:
            safe_regions.append({
                'chr': chrom,
                'start': 1,
                'end': merged[0][0] - 1,
                'size': merged[0][0] - 1
            })
        
        for i in range(len(merged)-1):
            safe_regions.append({
                'chr': chrom,
                'start': merged[i][1] + 1,
                'end': merged[i+1][0] - 1,
                'size': merged[i+1][0] - merged[i][1] - 1
            })
            
        if merged[-1][1] < chrom_size:
            safe_regions.append({
                'chr': chrom,
                'start': merged[-1][1] + 1,
                'end': chrom_size,
                'size': chrom_size - merged[-1][1]
            })
        
        print(f"Found {len(safe_regions)} safe regions for this chromosome")
        if safe_regions:
            print(f"First few safe regions: {safe_regions[:5]}")
    
    # Create DataFrame with all required columns even if empty
    result_df = pd.DataFrame(safe_regions, columns=['chr', 'start', 'end', 'size']) if safe_regions else pd.DataFrame(columns=['chr', 'start', 'end', 'size'])
    return result_df, feature_stats

def write_outputs(safe_regions_df, feature_stats, output_base, buffers, include_types, exclude_types):
    """Write output files"""
    # Add metadata columns to safe regions dataframe
    upstream_buffer, downstream_buffer = buffers
    
    # Create clear description of what features were considered
    safe_regions_df['features_considered'] = ','.join(sorted(include_types)) if include_types else 'all'
    safe_regions_df['upstream_buffer'] = upstream_buffer
    safe_regions_df['downstream_buffer'] = downstream_buffer
    
    # Write safe regions
    output_path = f"{output_base}.tsv"
    safe_regions_df.to_csv(output_path, sep='\t', index=False)
    print(f"\nWrote {len(safe_regions_df)} safe regions to {output_path}")
    
    # Write metadata
    metadata_path = f"{output_base}_metadata.txt"
    with open(metadata_path, 'w') as f:
        f.write("Safe Regions Analysis Metadata\n")
        f.write("=============================\n\n")
        f.write(f"Upstream buffer: {upstream_buffer}bp\n")
        f.write(f"Downstream buffer: {downstream_buffer}bp\n")
        f.write(f"Features considered: {','.join(sorted(include_types)) if include_types else 'all'}\n\n")
        
        f.write("Feature counts by chromosome:\n")
        f.write("-----------------------------\n")
        for chrom, counts in feature_stats.items():
            f.write(f"\n{chrom}:\n")
            for feature_type, count in counts.items():
                f.write(f"  {feature_type}: {count:,}\n")
        
        f.write("\nSafe region statistics:\n")
        f.write("----------------------\n")
        for chrom in safe_regions_df['chr'].unique():
            chrom_regions = safe_regions_df[safe_regions_df['chr'] == chrom]
            total_safe = chrom_regions['size'].sum()
            f.write(f"\n{chrom}:\n")
            f.write(f"  Number of safe regions: {len(chrom_regions):,}\n")
            f.write(f"  Total safe bases: {total_safe:,}\n")
            if len(chrom_regions) > 0:
                f.write(f"  Average region size: {total_safe/len(chrom_regions):,.0f}\n")
                f.write(f"  Smallest region: {chrom_regions['size'].min():,}\n")
                f.write(f"  Largest region: {chrom_regions['size'].max():,}\n")

def main():
    try:
        args = parse_arguments()
        
        # Process feature type arguments
        exclude_types = {'source'}
        if args.exclude_types:
            exclude_types.update(args.exclude_types.split(','))
        
        include_types = None
        if args.include_types:
            include_types = set(args.include_types.split(','))
        
        # Set output base name
        if not args.output:
            base = os.path.splitext(args.input)[0]
            output_base = f"{base}_safe_regions"
        else:
            # Remove any extension from provided output name
            output_base = os.path.splitext(args.output)[0]
        
        # Extract features
        print(f"Reading features from {args.input}...")
        features_df = get_features_from_genbank(args.input, include_types, exclude_types)
        print(f"Found {len(features_df)} features across {len(features_df['chr'].unique())} chromosomes")
        
        # Find safe regions
        print(f"\nIdentifying safe regions (upstream buffer: {args.upstream_buffer}bp, downstream buffer: {args.downstream_buffer}bp)...")
        safe_regions_df, feature_stats = find_safe_regions(
            features_df, 
            upstream_buffer=args.upstream_buffer,
            downstream_buffer=args.downstream_buffer
        )
        
        # Write outputs
        write_outputs(safe_regions_df, feature_stats, output_base, 
                     (args.upstream_buffer, args.downstream_buffer),
                     args.include_types, args.exclude_types)
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        import traceback
        print("\nFull traceback:")
        print(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()