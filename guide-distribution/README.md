# Guide RNA Selection Script

## Overview
This script selects and evenly distributes guide RNAs (gRNAs) across genes while accounting for constraints like gaps in guide availability. It employs a two-step approach combining gap analysis with simulated annealing optimization to achieve the most even possible distribution of guides.

## Usage
Basic usage:
```bash
python select_guides.py input.tsv --guides_per_pair 10
```

With custom parameters:
```bash
python select_guides.py input.tsv \
    --guides_per_pair 20 \
    --min_gap_fraction 0.2 \
    --min_guides_missing 2.0 \
    --max_iterations 2000 \
    --initial_temp 2.0 \
    --cooling_rate 0.98 \
    --processes 8
```

## Input
- TSV file containing guide RNA information
- Columns directly called:
   - locus_tag: Gene identifier
   - sp_dir: Strand of Spacer (F or R)
   - offset: Position relative to gene start (must be ≥ 0)

## Output
- Creates new TSV file with '_filtered' suffix
- Maintains all original columns
- Adds new columns:
   - original_guide_count: Number of guides available for that gene/orientation
   - selected_guide_count: Number of guides selected after filtering

## Distribution Approach

### Step 1: Gap-Based Initial Selection
The script will perform all of the described actions for guide candidates for each pair of locus_tag and sp_dir, to achieve evenly distributed guides across each gene on each strand.

The script first identifies significant gaps in guide coverage using two criteria:
- Gap size relative to gene length
- Gap size relative to ideal guide spacing

For example, if we want 10 guides across a 1000bp gene:
- Ideal spacing would be 100bp
- A gap of 200bp would represent 2 missing guide positions
- If this gap is also >15% of gene length, it's considered significant

### Step 2: Simulated Annealing Optimization
After initial selection, the script optimizes guide distribution using simulated annealing:
- Tries random guide swaps (keeping first and last guides fixed)
- Initially accepts some suboptimal swaps to avoid local minima
- Gradually becomes more selective
- Aims to minimize variance in guide spacing

### Example
Consider a gene with guides at these positions:
```
Position: 0, 10, 15, 20 |gap 200bp| 220, 230 |gap 400bp| 630, 640 |gap 150bp| 790, 800, 900, 1000
```

The algorithm would:
1. Identify significant gaps (200bp and 400bp gaps qualify)
2. Select guides at gene boundaries (0 and 1000)
3. Select guides flanking significant gaps (20, 220, 230, 630)
4. Use simulated annealing to optimize remaining selections
   - Try different combinations of guides
   - Evaluate evenness using standard deviation of distances
   - Accept better distributions immediately
   - Sometimes accept worse distributions early in the process
5. Final selection might be: 0, 20, 220, 630, 800, 900, 1000

## Parameters

### Gap Analysis Parameters
A gap must meet both of the following criteria to be considered significant:

1. `--min_gap_fraction` (default: 0.15)
   - Gap must be larger than this fraction of total gene length
   - Smaller values identify more gaps
   - Default means gaps must be >15% of gene length

2. `--min_guides_missing` (default: 1.5)
   - Gap must be larger than this many times the ideal guide spacing
   - Based on ideal spacing (gene_length / n_guides)
   - Default of 1.5 means gap must be 1.5x ideal spacing
   - For example, if ideal spacing is 100bp, gap must be >150bp

### Simulated Annealing Parameters
1. `--max_iterations` (default: 1000)
   - Maximum optimization attempts
   - Higher values may find better distributions but increase runtime

2. `--initial_temp` (default: 1.0)
   - Initial "temperature" for simulated annealing
   - Higher values increase initial randomness
   - Lower values make algorithm more selective from start

3. `--cooling_rate` (default: 0.95)
   - How quickly temperature decreases
   - Higher values (closer to 1) cool slower, allow more exploration
   - Lower values cool faster, settle on solution more quickly

4. `--max_no_improvement` (default: 100)
   - Stop if no improvement after this many iterations
   - Helps prevent unnecessary computation

### Processing Parameters
1. `--processes` (default: CPU count - 1)
   - Number of parallel processes to use
   - Defaults to using all but one CPU core

## Rationale for Approach
This hybrid approach was chosen because:
1. Simple even spacing doesn't account for variation in guide distribution across a given gene.
2. A gap-based selection can circumvent the problems of a simple, even-spacing approach, but may miss the optimal distribution.
3. Simulated annealing uses the gap-based guide selection and iteratively explores alternate distributions, allowing for an efficient search for the most even spacing.

## Evenness Score
The script quantifies how evenly distributed guides are using a standard deviation-based metric:

1. For each set of selected guides, calculate the distances between consecutive guides
2. Take the standard deviation of these distances
3. Lower scores indicate more even distribution

For example:
Perfect distribution (score = 0):
Positions: 0, 200, 400, 600, 800, 1000
Distances: 200, 200, 200, 200, 200

Uneven distribution (high score):
Positions: 0, 100, 150, 800, 900, 1000
Distances: 100, 50, 650, 100, 100

## Edge Cases
The script handles various edge cases:
- Genes with no guides
- Genes with fewer guides than requested
- Guides all at same position
- Large gaps in guide availability

## Limitations and Constraints

1. **Fixed End Points**
   - First and last guides are always selected and cannot be moved during optimization
   - This may occasionally prevent finding a more even distribution that excludes extreme positions

2. **Optimization vs Runtime**
   - Simulated annealing parameters balance finding optimal distributions with computational efficiency
   - More thorough exploration of possible distributions requires longer runtime
   - Early stopping when no improvement is found may occasionally miss better distributions

3. **Gap Analysis Trade-offs**
   - Gap identification requires meeting both size criteria (relative to gene length and ideal spacing)
   - This stringent approach might miss biologically relevant gaps that only meet one criterion
   - However, since gap-based selections can be modified during optimization, the impact is minimized

4. **Even Spacing vs Biological Context**
   - The evenness score only considers distances between consecutive guides
   - Does not account for other biological features like protein domains or regulatory regions
   - More sophisticated scoring could incorporate additional biological constraints

5. **Computational Resources**
   - Memory usage scales with number of guides per gene
   - Runtime increases with number of genes and requested guides
   - Parallel processing can strain system resources on large datasets
  
6. **Input Dataset**
   - The script was developed for the outputs of guide design tools from the Peters lab, which can be found here: https://github.com/ryandward/crispr_experiment/tree/main
   - The TSV with the output guide design does not include gene length or chromosomal coordinates, so gene length is approximated by the offset positions of the first and last available guide.