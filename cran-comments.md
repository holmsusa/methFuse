## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

Package: fuseR
Version: 1.0.0

fuseR implements the FUSE pipeline for genomic methylation segmentation, providing functions to cluster, segment, and summarize methylation data.

This is the first stable release (1.0.0). Major new features include:
- plot.fuse_summary() S3 method for visualization of segment-level methylation estimates
- Optional plotting of per-CpG beta values
- Improved input validation and S3 method registration
- Updated vignette with full example workflow

The package has been tested with R CMD check on R 4.5 and all examples and vignettes run successfully. Dependencies are limited to CRAN packages: stats, graphics, and rmarkdown.

We kindly submit this package for inclusion on CRAN.

