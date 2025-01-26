# CRISPR Guide Selection Pipeline

A set of scripts for selecting CRISPR guides based on their upstream location relative to genes and ensuring their insertion sites fall in feature-free regions. The pipeline consists of three steps:

## 1. Find Feature-Free Regions
```bash
find_feature_free_regions.py -i genome.gb -u UPSTREAM -d DOWNSTREAM
```
Identifies regions in the genome free of any genomic features, considering strand-specific buffers. Outputs a TSV file of safe regions and detailed metadata.

## 2. Filter Guides by Safe Regions
```bash
filter_guides_by_safe_regions.py -g guides.tsv -f genome_safe_regions.tsv
```
Takes a guide file and filters for guides whose predicted insertion sites (±82bp depending on direction) fall within feature-free regions.

## 3. Filter by Upstream Location
```bash
filter_guides_by_upstream.py -g guides_safe_filtered.tsv -c coordinates.tsv -b genome.gb -p PROXIMAL -d DISTAL
```
Filters guides to keep only those that fall within a specified window upstream of genes (default: 102-302bp).

## Pipeline Features
- Handles features that cross chromosome origin
- Considers strand-specific buffers around features
- Accounts for guide direction in insertion site prediction
- Provides detailed statistics and metadata for each step
- Supports parallel processing where applicable

## Notes
- Feature-free filtering should be run once per genome/buffer combination
- Upstream filtering parameters can be adjusted without rerunning feature-free analysis
- All coordinates use 1-based GenBank-style notation