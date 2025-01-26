# Changelog

## [2.0.0] - 2025-01-25

### Major Changes
- Separated functionality into three distinct scripts for better modularity and reusability
- Added prediction and validation of transposon insertion sites
- Removed requirement for guides to be intergenic
- Added proper handling of features that cross chromosome origin

### New Scripts
- `find_feature_free_regions.py`: New script to identify genomic regions free of features
- `filter_guides_by_safe_regions.py`: New script to filter guides based on insertion sites
- Renamed original script to `filter_guides_by_upstream.py` for clarity

### Feature Additions
- Added strand-aware buffer zones around features
- Added predicted insertion site calculations (±82bp based on guide direction)
- Added comprehensive metadata output for feature-free regions
- Added detailed progress reporting for each step

### Performance Improvements
- Optimized guide filtering using numpy arrays for position comparisons
- Added chromosome-by-chromosome processing for memory efficiency
- Separated feature-free analysis to avoid redundant computations

### Changes from Original Workflow
#### Removed
- Intergenic-only guide filtering
- Direct insertion site checks during upstream filtering

#### Modified
- Split processing into separate steps
- Changed output file naming conventions
- Updated progress reporting format

### Technical Details
- Forward guide insertion sites: tar_start + 82bp
- Reverse guide insertion sites: tar_end - 82bp
- Origin-crossing features are now properly split and processed
- All coordinates use 1-based GenBank-style notation

### Notes
- This is a breaking change from the previous single-script workflow
- Output files maintain same format but include additional insertion site information
- Existing bp_upstream calculations remain unchanged for compatibility