# sgRNA Filtering Tools

A collection of tools for filtering and optimizing sgRNA/gRNA selections for CRISPR applications.

## Overview

This repository contains various tools to help researchers optimize their guide RNA selections for CRISPR experiments. Each tool is contained in its own directory with detailed documentation on usage and implementation.

## Available Tools

### guide-distribution

Located in: `/guide-distribution`

A Python-based tool for selecting and optimally distributing guide RNAs across genes. The tool uses a sophisticated two-step approach:
1. Gap analysis to identify regions without guide coverage
2. Simulated annealing optimization to achieve the most even possible distribution of guides

Key features:
- Even distribution of gRNAs across target genes
- Accounts for gaps in guide availability
- Optimization through simulated annealing
- Input format: TSV file containing sgRNA data

For detailed usage instructions, see the [guide-distribution README](guide-distribution/README.md).

## Tool Directory Structure

Each tool directory contains:
- A Python script implementing the tool's functionality
- A README with detailed usage instructions

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
