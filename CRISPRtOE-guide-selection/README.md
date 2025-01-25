# get_CRISPRtOE_guides.py

This script filters sgRNA guides based on their position relative to genes, specifically identifying guides that lie in upstream regions. It is particularly useful for finding guides that could be used for CRISPRi/CRISPRa applications where targeting the upstream region of a gene is desired.

## Prerequisites

- Python 3.x
- Required Python packages:
  - pandas
  - biopython
  - numpy

## Installation

1. Clone this repository or download the script
2. Install required packages:
```bash
pip install pandas biopython numpy
```

## Input Files

The script requires three input files:

1. **Guides TSV file** (required): A tab-separated file containing sgRNA guide information with columns including:
   - spacer: Guide sequence
   - chr: Chromosome ID
   - tar_start, tar_end: Target coordinates
   - sp_dir: Guide direction (F/R)
   - intergenic: Boolean (1/0) indicating if the guide is in an intergenic region
   - genes: Number of genes the guide targets (should be 0 for valid guides)
   - (other optional columns that will be preserved in the output)

2. **Coordinates TSV file** (required): A tab-separated file containing gene coordinates with columns:
   - chr: Chromosome ID
   - locus_tag: Gene locus tag
   - gene: Gene name
   - gene_dir: Gene direction (F/R)
   - feature_start: Start coordinate
   - feature_end: End coordinate
   - feature_type: Type of feature

3. **GenBank file** (required): The GenBank file that was used to generate the guides TSV. This is needed to determine chromosome sizes for proper coordinate handling.

4. **Loci file** (optional): A tab-separated file with a header of either 'locus_tag' or 'gene' containing the identifiers of genes of interest.

## Usage

Basic usage:
```bash
python get_CRISPRtOE_guides.py -g guides.tsv -c coordinates.tsv -b genome.gb
```

With all optional parameters:
```bash
python get_CRISPRtOE_guides.py \
  -g guides.tsv \
  -c coordinates.tsv \
  -b genome.gb \
  -l genes_of_interest.tsv \
  -o custom_output.tsv \
  -p 150 \
  -d 350 \
  -t 8
```

### Command Line Arguments

- `-g, --guides`: Path to input guides TSV file (required)
- `-c, --coordinates`: Path to coordinates TSV file (required)
- `-b, --genbank`: Path to GenBank file (required)
- `-l, --loci`: Path to file containing genes/loci of interest (optional)
- `-o, --output`: Path for output file (default: input_name_upstream_filtered.tsv)
- `-p, --proximal`: Distance from gene start for proximal boundary (default: 102)
- `-d, --distal`: Distance from gene start for distal boundary (default: 302)
- `-t, --threads`: Number of processor threads to use (default: auto)

## Output Files

The script generates several output files:

1. **Main output** (`*_upstream_filtered.tsv`): Contains filtered guides meeting these criteria:
   - Guide is in an intergenic region (intergenic=1)
   - Guide doesn't target any genes (genes=0)
   - Guide direction matches gene direction
   - Guide position is within specified upstream region
   
   Added columns in output:
   - upstream_of_locus: Locus tag of the gene the guide is upstream of
   - upstream_of_gene: Gene name the guide is upstream of
   - gene_dir: Direction of the gene
   - gene_start: Start coordinate of the gene
   - gene_end: End coordinate of the gene
   - bp_upstream: Distance in base pairs from the start of the guide to the start of the gene (taking into account both gene and guide orientation)

2. **Summary file** (`*_summary.tsv`): Contains statistics about guides found per gene:
   - upstream_of_locus: Locus tag
   - upstream_of_gene: Gene name
   - guide_count: Number of guides found upstream of this gene
   
   Note: This file includes ALL genes from the coordinates file (or subset if loci file provided), even those with zero guides found.

3. **Missing loci file** (`*_missing_loci.txt`, optional): Created only if a loci file was provided, lists any loci that weren't found in the coordinates file.

## Notes

- The script handles features that cross the origin of circular chromosomes
- The upstream region is calculated based on gene orientation:
  - For forward (F) genes: upstream is before the start coordinate
  - For reverse (R) genes: upstream is after the end coordinate
- Guide direction (sp_dir) must match gene direction (gene_dir) for a guide to be considered
- Only guides marked as intergenic (intergenic=1) and not targeting any genes (genes=0) are considered
- Parallel processing is used to improve performance with large datasets

## Error Handling

The script includes comprehensive error checking for:
- Missing input files
- Mismatched chromosome IDs between files
- Invalid file formats
- Missing required columns
- Origin-crossing coordinates

Error messages will indicate the specific issue and where it occurred in the process.