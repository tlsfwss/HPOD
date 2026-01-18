## Stable age-associated gene identification

This directory provides the R implementation of a stability-aware framework
used in the HPOD database to identify robust age-associated genes from
large-scale human transcriptomic data.
The pipeline is designed to be reusable for independent human transcriptomic
datasets with continuous age annotations.

### Method overview

The pipeline integrates:

* edgeR TMM normalization and voom transformation
* natural spline modeling of age effects (ns(Age, df))
* joint F-tests on spline terms using limma
* stratified resampling-based stability selection
* sign-consistency filtering based on Spearman correlations
* full-data refitting for final significance control

Model complexity (spline degrees of freedom) is selected by maximizing the
mean pairwise Jaccard similarity of selected gene sets across resampling
iterations, favoring stable and parsimonious models.

### Reproducibility

All random seeds are fixed to ensure reproducibility.
A one-click script (`run_stable_age_pipeline_oneclick.R`) is provided to
reproduce the full analysis workflow and associated figures.
