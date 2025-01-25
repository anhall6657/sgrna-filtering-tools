# Build Coordinates File

This script (`build_coordinates_file.py`) extracts feature coordinates from a GenBank file and matches them with locus tags from a guide TSV file. It was designed to support CRISPR guide RNA design by providing accurate genomic coordinates for features of interest.

## Purpose

The script serves to:
1. Extract coordinate information for genomic features from a GenBank file
2. Match these coordinates with locus tags from a guide TSV file
3. Handle special cases such as:
   - Features that cross the origin of replication
   - Features on both forward and reverse strands
   - Multiple chromosomes/records in a single GenBank file

## Input Files

### Guide TSV File
- Tab-separated file containing guide RNA information
- Required columns:
  - `chr`: Chromosome identifier
  - `locus_tag`: Unique identifier for genomic features
  - `gene`: Gene name
  - `tar_dir`: Target direction (will be renamed to `gene_dir` in output)

### GenBank File
- Standard GenBank format file containing feature annotations
- Must contain feature locations and locus tags that correspond to those in the guide TSV

## Usage

Basic usage:
```bash
python build_coordinates_file.py path/to/guide.tsv path/to/genome.gb
```

With custom output filename:
```bash
python build_coordinates_file.py path/to/guide.tsv path/to/genome.gb -o custom_output.tsv
```

## Output

The script generates a TSV file containing:
- All matched features from the input guide TSV
- Additional columns:
  - `feature_start`: Start coordinate from GenBank
  - `feature_end`: End coordinate from GenBank
  - `feature_type`: Type of feature (e.g., CDS, gene, etc.)
- Output filename format: `GENBANKID_coordinates.tsv` (unless specified with -o)

## Installation

1. Ensure Python 3.6 or higher is installed
2. Install required packages:
```bash
pip install -r requirements.txt
```

## Requirements

See `requirements.txt` for specific version information. Key dependencies:
- pandas: For data frame operations
- BioPython: For GenBank file parsing

## Error Handling

The script will:
- Check for required columns in the guide TSV
- Remove rows with empty/None locus tags or genes
- Skip features where coordinates cannot be found
- Provide progress updates and statistics about processed features