## Stable age-associated gene identification

This repository provides the R implementation of a stability-aware framework
for identifying robust age-associated genes from large-scale human transcriptomic data,
as used in the HPOD database.

### Method overview
The pipeline integrates:
- edgeR TMM normalization and voom transformation
- natural spline modeling of age effects (ns(Age, df))
- joint F-tests on spline terms using limma
- stratified resampling-based stability selection
- sign-consistency filtering based on Spearman correlations
- full-data refitting for final significance control

Model complexity (spline df) is selected based on maximal mean pairwise Jaccard
stability across resampling iterations.

### Reproducibility
All random seeds are fixed. A one-click script is provided to reproduce
the full analysis and figures.

