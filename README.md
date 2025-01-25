# sgRNA Filtering Tools

A collection of tools for filtering and optimizing sgRNA/gRNA selections for CRISPR applications.

## Overview

This repository contains various tools to help researchers optimize their guide RNA selections for CRISPR experiments. Each tool is contained in its own directory with detailed documentation on usage and implementation.

## Available Tools

### guide-distribution

Located in: `/guide-distribution`

A Python-based tool for selecting and optimally distributing guide RNAs across genes. Built for CRISPRt library design, but could be applied to CRISPRi library design. The tool uses a sophisticated two-step approach:
1. Gap analysis to identify regions without guide coverage
2. Simulated annealing optimization to achieve the most even possible distribution of guides

Key features:
- Even distribution of sgRNAs across target genes
- Accounts for gaps in guide availability
- Optimization through simulated annealing
- Input format: TSV file containing sgRNA data

For detailed usage instructions, see the [guide-distribution README](guide-distribution/README.md).

### get-chromosomal-coordinates

Located in: `/get-chromosomal-coordinates`

A Python-based tool for extracting chromosomal coordinates of genetic features from a GenBank file. This is a required precursor for CRISPRtOE guide selection, given that the original sgRNA TSV file does not contain this information.

Key features:
- Takes sgRNA TSV file as input
- Finds chromosomal coordinates for all features in the input TSV file by finding features by locus_tag in a corresponding GenBank file
- Outputs a coordinates.tsv file with the columns chr, locus_tag, gene, gene_dir, feature_start, feature_end, and feature_type
- This output is used as an input in the CRISPRtOE guide selection code

For detailed usage instructions, see the [get-chromosomal-coordinates README](get-chromosomal-coordinates/README.md).

### CRISPRtOE-guide-selection

Located in: `/CRISPRtOE-guide-selection`

A Python-based tool for finding sgRNAs that could be used in CRISPRtOE experiments. This tool finds guides in regions upstream of genes that are on the same strand as the gene of interest. The default region of interest is defined as 102-302bp upstream of the start of the gene.

Key features:
- Removes sgRNAs that target within a gene. CHANGE: Need to exclude guides that would result in transposition into a gene, _not_ those that match to a gene themselves. Site of transposition is what matters.
- Ensures guide direction and gene direction match
- Outputs all guides that meet the criteria of interest, with added columns that describe the gene the guide is upstream of.


For detailed usage instructions, see the [CRISPRtOE-guide-selection README](CRISPRtOE-guide-selection/README.md).

## Tool Directory Structure

Each tool directory contains:
- A Python script implementing the tool's functionality
- A README with detailed usage instructions

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
