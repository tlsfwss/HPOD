#!/usr/bin/env Rscript
# ============================================================
# One-click runner: stable age genes + plots
# - Put these files in the same folder:
#     1) stable_age_genes_final.R
#     2) stable_age_genes_plots.R
#     3) counts_mat.txt
#     4) pheno.txt
#
# Run (terminal):
#   Rscript run_stable_age_pipeline_oneclick.R
#
# Or in R:
#   source("run_stable_age_pipeline_oneclick.R")
#
# ============================================================

# --------- User settings ----------
counts_file <- "counts_mat.txt"
pheno_file  <- "pheno.txt"
outdir      <- "stable_age_outputs_final"

B <- 150
df_grid <- c(2,3,4,5)

do_prefilter <- TRUE
prefilter_top_prop <- 0.10
prefilter_min_abs_rho <- 0.20

alpha <- 0.05
use_sex <- FALSE

freq_thresh <- 0.80
sign_consistency_thresh <- 0.90
require_full_refit_fdr <- 0.05

seed <- 123

top_n_trends <- 12
# -----------------------------------------------

# Load functions
source("stable_age_genes_final.R")
source("stable_age_genes_plots.R")

# 1) Run stability selection
res <- run_stable_age_genes_final(
  counts_file = counts_file,
  pheno_file  = pheno_file,
  outdir      = outdir,
  B = B,
  df_grid = df_grid,
  do_prefilter = do_prefilter,
  prefilter_top_prop = prefilter_top_prop,
  prefilter_min_abs_rho = prefilter_min_abs_rho,
  alpha = alpha,
  use_sex = use_sex,
  freq_thresh = freq_thresh,
  sign_consistency_thresh = sign_consistency_thresh,
  require_full_refit_fdr = require_full_refit_fdr,
  seed = seed
)

# 2) Generate plots
plots_res <- make_stability_plots(
  outdir       = outdir,
  counts_file  = counts_file,
  pheno_file   = pheno_file,
  top_n_trends = top_n_trends
)

message("\nAll done.")
message("Results: ", normalizePath(outdir))
message("Plots:   ", normalizePath(file.path(outdir, "plots")))
