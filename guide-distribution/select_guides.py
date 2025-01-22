#!/usr/bin/env python3

"""
# Guide RNA Selection Script

## Overview
This script selects guide RNAs (gRNAs) for genes using a two-step approach:
1. Gap-based initial selection
2. Distribution optimization using simulated annealing

## Approach

### Step 1: Gap-based Selection
The script first identifies significant gaps in guide coverage using two criteria:
- Gap size relative to gene length
- Gap size relative to ideal guide spacing

For example, if we want 10 guides across a 1000bp gene:
- Ideal spacing would be 100bp
- A gap of 200bp would represent 2 missing guide positions
- If this gap is also >15% of gene length, it's considered significant

### Step 2: Simulated Annealing Optimization
After initial selection, the script optimizes guide distribution using simulated annealing:
- Tries random guide swaps
- Initially accepts some suboptimal swaps
- Gradually becomes more selective
- Aims to minimize variance in guide spacing

## Parameters

### Gap Analysis Parameters
1. `--min_gap_fraction` (default: 0.15)
   - What fraction of total gene length a gap must exceed to be considered significant
   - Default means gaps must be >15% of gene length
   - Lower values will identify more gaps
   - Higher values will only catch very large gaps

2. `--min_guides_missing` (default: 1.5)
   - How many ideal guide positions must be missing to consider it a gap
   - Based on ideal spacing (gene_length / n_guides)
   - Default of 1.5 means gap must be 1.5x ideal spacing
   - Lower values will identify more gaps
   - Higher values will only catch gaps missing many potential guides

### Simulated Annealing Parameters
1. `--max_iterations` (default: 1000)
   - Maximum number of optimization attempts
   - Higher values give more chances to find optimal distribution
   - But increase runtime
   - Default balances optimization and runtime

2. `--initial_temp` (default: 1.0)
   - Initial "temperature" for simulated annealing
   - Higher values increase initial randomness
   - Lower values make algorithm more selective from start
   - Default allows good initial exploration

3. `--cooling_rate` (default: 0.95)
   - How quickly temperature decreases
   - Higher values (closer to 1) cool slower, allow more exploration
   - Lower values cool faster, settle on solution more quickly
   - Default allows thorough exploration while ensuring convergence

### General Parameters
1. `--guides_per_pair` (default: 10)
   - Number of guides to select per gene/orientation pair
   - Affects both gap analysis and final selection

## Usage
Basic usage:
```bash
python select_guides.py input.tsv
```

With custom parameters:
```bash
python select_guides.py input.tsv \
    --guides_per_pair 20 \
    --min_gap_fraction 0.2 \
    --min_guides_missing 2.0 \
    --max_iterations 2000 \
    --initial_temp 2.0 \
    --cooling_rate 0.98
```

## Output
- Creates new TSV file with '_filtered' suffix
- Adds columns for original and selected guide counts
- Prints summary statistics including:
  - Original number of guides
  - Number of guides after filtering
  - Number of unique genes processed
  - Average guides per gene/strand before and after selection

Additional Parameters:
--processes : Number of CPU processes to use (default: total CPUs - 1)
--max_no_improvement : Number of iterations without improvement before early stopping (default: 100)
"""

import pandas as pd
import argparse
import numpy as np
import os
import random
from typing import List, Tuple
from multiprocessing import Pool, cpu_count
from functools import partial

def parse_args():
    parser = argparse.ArgumentParser(description='Filter guide RNAs based on gene and orientation')
    parser.add_argument('input_file', help='Input TSV file')
    
    # General parameters
    parser.add_argument('--guides_per_pair', type=int, default=10,
                       help='Number of guides to select per gene/orientation pair')
    
    # Gap analysis parameters
    parser.add_argument('--min_gap_fraction', type=float, default=0.15,
                       help='Minimum gap size as fraction of gene length')
    parser.add_argument('--min_guides_missing', type=float, default=1.5,
                       help='Number of ideal guide positions that need to be missing to consider it a gap')
    
    # Simulated annealing parameters
    parser.add_argument('--max_iterations', type=int, default=1000,
                       help='Maximum number of iterations for optimization')
    parser.add_argument('--initial_temp', type=float, default=1.0,
                       help='Initial temperature for simulated annealing')
    parser.add_argument('--cooling_rate', type=float, default=0.95,
                       help='Cooling rate for simulated annealing')
    parser.add_argument('--max_no_improvement', type=int, default=100,
                       help='Stop if no improvement after this many iterations')
    
    # Parallel processing parameters
    parser.add_argument('--processes', type=int, default=max(1, cpu_count() - 1),
                       help='Number of processes to use')
    
    return parser.parse_args()

def generate_output_filename(input_file):
    """Generate output filename by adding '_filtered' before the extension."""
    base, ext = os.path.splitext(input_file)
    return f"{base}_filtered{ext}"

def find_significant_gaps(positions: List[float], gene_length: float, n_guides: int, 
                         min_gap_fraction: float, min_guides_missing: float) -> List[Tuple[float, float]]:
    """
    Find significant gaps in guide positions.
    Returns list of (start, end) positions for significant gaps.
    """
    ideal_spacing = gene_length / n_guides
    gaps = []
    for i in range(len(positions) - 1):
        gap_size = positions[i+1] - positions[i]
        if (gap_size > (ideal_spacing * min_guides_missing) and 
            gap_size > (gene_length * min_gap_fraction)):
            gaps.append((positions[i], positions[i+1]))
    
    return gaps

def calculate_evenness_score(positions: List[float], gene_length: float) -> float:
    """Calculate how evenly distributed the guides are."""
    if len(positions) < 2:
        return float('inf')
    
    distances = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
    return np.std(distances)

def initial_gap_based_selection(group: pd.DataFrame, n_guides: int, 
                              min_gap_fraction: float, min_guides_missing: float) -> List[pd.Series]:
    """Make initial guide selection based on gaps."""
    positions = sorted(group['offset'].values)
    
    # Handle cases with too few guides
    if len(positions) == 0:
        return []
    if len(positions) == 1:
        return [group.iloc[0]]
    
    gene_length = positions[-1] - positions[0]
    
    # If gene_length is 0 (all guides at same position), return first guide
    if gene_length == 0:
        return [group.iloc[0]]
    
    
    # Find significant gaps
    gaps = find_significant_gaps(positions, gene_length, n_guides, 
                               min_gap_fraction, min_guides_missing)
    
    # Select guides at gap edges
    selected_indices = set()
    selected_guides = []
    
    # Always take first and last guides
    first_idx = group['offset'].idxmin()
    last_idx = group['offset'].idxmax()
    selected_indices.add(first_idx)
    selected_indices.add(last_idx)
    selected_guides.extend([group.loc[first_idx], group.loc[last_idx]])
    
    # Handle gaps
    for gap_start, gap_end in gaps:
        if len(selected_guides) >= n_guides:
            break
            
        # Find guides closest to gap edges
        start_guide_idx = group['offset'].sub(gap_start).abs().idxmin()
        end_guide_idx = group['offset'].sub(gap_end).abs().idxmin()
        
        for idx in [start_guide_idx, end_guide_idx]:
            if idx not in selected_indices and len(selected_guides) < n_guides:
                selected_indices.add(idx)
                selected_guides.append(group.loc[idx])
    
    # Fill remaining positions if needed
    remaining_guides = group[~group.index.isin(selected_indices)]
    while len(selected_guides) < n_guides and len(remaining_guides) > 0:
        next_guide_idx = remaining_guides['offset'].sample(1).index[0]
        selected_guides.append(remaining_guides.loc[next_guide_idx])
        remaining_guides = remaining_guides.drop(next_guide_idx)
    
    return selected_guides

def optimize_distribution(initial_guides: List[pd.Series], 
                        available_guides: pd.DataFrame,
                        n_guides: int,
                        max_iterations: int,
                        initial_temp: float,
                        cooling_rate: float,
                        max_no_improvement: int) -> List[pd.Series]:
    """Use simulated annealing to optimize guide distribution."""
    current_solution = initial_guides
    current_score = calculate_evenness_score([g['offset'] for g in current_solution], 
                                           available_guides['offset'].max())
    best_solution = current_solution
    best_score = current_score
    temperature = initial_temp
    no_improvement_count = 0
    
    for iteration in range(max_iterations):
        if no_improvement_count >= max_no_improvement:
            break
            
        # Try swapping a random guide
        if len(current_solution) < n_guides:
            continue
            
        # Don't swap first or last guides
        swap_idx = random.randint(1, len(current_solution)-2)
        available_positions = available_guides[
            ~available_guides.index.isin([g.name for g in current_solution])
        ]
        
        if len(available_positions) == 0:
            continue
            
        new_guide = available_positions.sample(1).iloc[0]
        new_solution = current_solution.copy()
        new_solution[swap_idx] = new_guide
        
        # Calculate new score
        new_score = calculate_evenness_score([g['offset'] for g in new_solution],
                                           available_guides['offset'].max())
        
        # Decide whether to accept new solution
        if new_score < current_score or random.random() < np.exp((current_score - new_score) / temperature):
            current_solution = new_solution
            current_score = new_score
            if current_score < best_score:
                best_solution = current_solution
                best_score = current_score
                no_improvement_count = 0
            else:
                no_improvement_count += 1
        else:
            no_improvement_count += 1
        
        temperature *= cooling_rate
    
    return best_solution
    
def process_group(group_data: Tuple[Tuple, pd.DataFrame], args) -> pd.DataFrame:
    """Process a single group of guides. Used for parallel processing."""
    name, group = group_data
    
    # Store original count
    original_count = len(group)
    
    # Handle empty groups or groups with no valid guides
    if len(group) == 0:
        # Create empty DataFrame with correct columns
        columns = list(group.columns) + ['original_guide_count', 'selected_guide_count']
        return pd.DataFrame(columns=columns)
    
    # Only consider non-negative offsets
    group = group[group['offset'] >= 0].copy()
    
    # If no guides remain after filtering negatives
    if len(group) == 0:
        # Create empty DataFrame with correct columns
        columns = list(group.columns) + ['original_guide_count', 'selected_guide_count']
        return pd.DataFrame(columns=columns)
    
    # If we have fewer guides than requested, return all of them
    if len(group) <= args.guides_per_pair:
        group['original_guide_count'] = original_count
        group['selected_guide_count'] = len(group)
        return group
    
    # Convert offset to integer type
    group['offset'] = group['offset'].astype(np.int32)
    
    # Get initial selection using gap-based approach
    initial_selection = initial_gap_based_selection(
        group, 
        args.guides_per_pair,
        args.min_gap_fraction,
        args.min_guides_missing
    )
    
    # If initial selection failed or returned too few guides,
    # just take the first n guides
    if len(initial_selection) < min(args.guides_per_pair, len(group)):
        initial_selection = [group.iloc[i] for i in range(min(args.guides_per_pair, len(group)))]
    
    # Only proceed with optimization if we have enough guides
    if len(initial_selection) > 2:
        final_selection = optimize_distribution(
            initial_selection, 
            group, 
            args.guides_per_pair,
            args.max_iterations,
            args.initial_temp,
            args.cooling_rate,
            args.max_no_improvement
        )
    else:
        final_selection = initial_selection
    
    # Create result DataFrame
    result = pd.concat([pd.DataFrame([g]) for g in final_selection], ignore_index=True)
    result['original_guide_count'] = original_count
    result['selected_guide_count'] = len(final_selection)
    
    return result

def main():
    args = parse_args()
    output_file = generate_output_filename(args.input_file)
    
    print(f"Reading input file: {args.input_file}")
    # Read TSV file
    df = pd.read_csv(args.input_file, sep='\t')
    
    # Remove rows where locus_tag is empty/NA/None
    df = df.dropna(subset=['locus_tag'])
    df = df[df['locus_tag'].str.strip() != '']
    
    # Convert offset to integer type
    df['offset'] = df['offset'].astype(np.int32)
    
    # Sort by locus_tag and sp_dir
    df = df.sort_values(['locus_tag', 'sp_dir', 'offset'])
    
    # Prepare groups for parallel processing
    groups = list(df.groupby(['locus_tag', 'sp_dir']))
    
    print(f"Processing {len(groups)} gene/orientation pairs using {args.processes} processes")
    
    # Process groups in parallel with error handling
    with Pool(processes=args.processes) as pool:
        try:
            results = pool.map(partial(process_group, args=args), groups)
        except Exception as e:
            print(f"Error during processing: {str(e)}")
            raise
    
    # Filter out empty results and combine
    results = [r for r in results if not r.empty]
    if not results:
        print("No guides were selected. Check input data.")
        return
        
    result = pd.concat(results, ignore_index=True)
    
    # Save to output file
    result.to_csv(output_file, sep='\t', index=False)
    
    # Print summary
    print(f"\nOriginal number of guides: {len(df)}")
    print(f"Number of guides after filtering: {len(result)}")
    print(f"Number of unique genes processed: {len(df['locus_tag'].unique())}")
    print("\nDistribution statistics:")
    print(f"Average original guides per gene/strand: {result['original_guide_count'].mean():.1f}")
    print(f"Average selected guides per gene/strand: {result['selected_guide_count'].mean():.1f}")
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    main()