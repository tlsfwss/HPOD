#!/usr/bin/env Rscript
# ============================================================
# Final stable age-associated gene identification (limma-voom)
# - Inputs:
#   1) counts_mat.txt : gene x sample raw counts (first col = gene; colnames = GSM IDs)
#   2) pheno.txt      : sample phenotype table with columns:
#        GSM.ID, GSE.ID, Age, Gender, Group (others allowed)
# - Core:
#   edgeR TMM -> voom -> limma linear model
#   resampling-based stability selection across df grid
# - Outputs (outdir):
#   df_stability_summary.csv
#   df{d}_selection_frequency.csv
#   df{d}_stable_genes.txt
#   final_stable_genes.csv / .txt
#   final_full_refit_table.csv
#   final_stable_genes_direction.csv
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(splines)
  library(dplyr)
  library(readr)
  library(tibble)
})

# -----------------------------
# Utility helpers
# -----------------------------

align_counts_pheno <- function(counts_mat, pheno) {
  stopifnot("SampleID" %in% names(pheno))
  common <- intersect(colnames(counts_mat), pheno$SampleID)
  if (length(common) < 10) stop("Too few overlapping samples between counts and pheno.")
  counts_mat <- counts_mat[, common, drop = FALSE]
  pheno <- pheno[match(common, pheno$SampleID), , drop = FALSE]
  stopifnot(all(colnames(counts_mat) == pheno$SampleID))
  list(counts_mat = counts_mat, pheno = pheno)
}

make_age_group <- function(age) {
  cut(age,
      breaks = c(-Inf, 39, 59, Inf),
      labels = c("Young", "Mid", "Old"),
      right = TRUE)
}

# Stratified subsampling: keep approximate distribution for age strata + study
stratified_subsample <- function(pheno, train_frac = 0.8,
                                 strata_cols = c("AgeGroup", "StudyID"),
                                 min_per_stratum_train = 2) {
  if (!all(strata_cols %in% names(pheno))) stop("Missing strata columns in pheno: ", paste(setdiff(strata_cols, names(pheno)), collapse=", "))
  pheno$.strata <- do.call(paste, c(pheno[strata_cols], sep = "||"))
  idx <- seq_len(nrow(pheno))

  train_idx <- unlist(lapply(split(idx, pheno$.strata), function(ii) {
    n <- length(ii)
    k <- max(min_per_stratum_train, floor(n * train_frac))
    k <- min(k, n) 
    sample(ii, size = k, replace = FALSE)
  }))

  train_idx <- sort(unique(train_idx))
  test_idx <- setdiff(idx, train_idx)

  # ensure test has enough samples
  if (length(test_idx) < max(5, floor(nrow(pheno) * (1 - train_frac)))) {
    test_n <- max(5, floor(nrow(pheno) * (1 - train_frac)))
    test_idx <- sample(idx, size = test_n, replace = FALSE)
    train_idx <- setdiff(idx, test_idx)
  }

  list(train = train_idx, test = test_idx)
}

# Spearman prefilter within TRAIN only (reduces noise genes)
spearman_prefilter <- function(dge_train, pheno_train, top_prop = 0.10, min_abs_rho = 0.20, min_keep = 300) {
  expr <- cpm(dge_train, log = TRUE, prior.count = 1)
  age <- pheno_train$Age

  rho <- apply(expr, 1, function(x) suppressWarnings(cor(x, age, method = "spearman")))
  pval <- apply(expr, 1, function(x) cor.test(x, age, method = "spearman")$p.value)
  fdr <- p.adjust(pval, method = "BH")

  df <- data.frame(gene = rownames(expr), rho = rho, fdr = fdr, stringsAsFactors = FALSE) %>%
    filter(!is.na(rho))

  df2 <- df %>% filter(abs(rho) >= min_abs_rho)
  if (nrow(df2) < min_keep) df2 <- df

  keep_n <- max(min_keep, floor(nrow(df2) * top_prop))
  df2 <- df2 %>% arrange(desc(abs(rho)))
  head(df2$gene, keep_n)
}

# Fit ns(Age, df) on TRAIN and select by joint F-test on spline terms
fit_select_ns_age <- function(dge_train, pheno_train, df_ns = 3,
                              use_sex = FALSE,
                              use_cells = FALSE, cellpc_cols = c("CellPC1", "CellPC2"),
                              alpha = 0.05) {
  # Build RHS terms
  terms <- c(sprintf("ns(Age, df=%d)", df_ns), "StudyID")
  if (use_sex) terms <- c(terms, "Sex")
  if (use_cells) {
    miss <- setdiff(cellpc_cols, names(pheno_train))
    if (length(miss) > 0) stop("Missing cell PC columns: ", paste(miss, collapse = ", "))
    terms <- c(terms, paste0("scale(", cellpc_cols, ")"))
  }

  design <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = pheno_train)

  # If singular, try dropping Sex first, then StudyID (last resort)
  if (!limma::is.fullrank(design)) {
    if (use_sex) {
      terms2 <- setdiff(terms, "Sex")
      design <- model.matrix(as.formula(paste("~", paste(terms2, collapse = " + "))), data = pheno_train)
      if (!limma::is.fullrank(design)) {
        # extreme fallback
        terms3 <- setdiff(terms2, "StudyID")
        design <- model.matrix(as.formula(paste("~", paste(terms3, collapse = " + "))), data = pheno_train)
        if (!limma::is.fullrank(design)) stop("Design matrix is not full rank; please check covariates.")
      }
    } else {
      terms3 <- setdiff(terms, "StudyID")
      design <- model.matrix(as.formula(paste("~", paste(terms3, collapse = " + "))), data = pheno_train)
      if (!limma::is.fullrank(design)) stop("Design matrix is not full rank; please check covariates.")
    }
  }

  v <- voom(dge_train, design, plot = FALSE)
  fit <- lmFit(v, design)

  age_cols <- grep("^ns\\(Age", colnames(design))
  if (length(age_cols) < 1) stop("No ns(Age, ...) columns found in design.")

  L <- matrix(0, nrow = ncol(design), ncol = length(age_cols))
  L[age_cols, ] <- diag(length(age_cols))
  rownames(L) <- colnames(design)

  fit2 <- contrasts.fit(fit, L)
  fit2 <- eBayes(fit2)

  tab <- topTable(fit2, number = Inf, sort.by = "F")
  sig <- tab %>% rownames_to_column("gene") %>% filter(adj.P.Val < alpha)
  selected <- sig$gene

  # direction from Spearman in TRAIN (for sign consistency)
  sign_dir <- integer(0)
  if (length(selected) > 0) {
    expr <- cpm(dge_train, log = TRUE, prior.count = 1)
    age <- pheno_train$Age
    rho <- apply(expr[selected, , drop = FALSE], 1, function(x) suppressWarnings(cor(x, age, method = "spearman")))
    sign_dir <- sign(rho)
    names(sign_dir) <- selected
  }

  list(selected_genes = selected, sign_dir = sign_dir)
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

mean_pairwise_jaccard <- function(sets, max_pairs = 2000) {
  n <- length(sets)
  if (n <= 1) return(NA_real_)
  pairs <- combn(n, 2)
  m <- ncol(pairs)
  if (m > max_pairs) {
    sel <- sample(seq_len(m), size = max_pairs, replace = FALSE)
    pairs <- pairs[, sel, drop = FALSE]
  }
  sims <- apply(pairs, 2, function(ii) jaccard(sets[[ii[1]]], sets[[ii[2]]]))
  mean(sims, na.rm = TRUE)
}

# Full-data refit: return joint-F table for ns(Age)
full_refit_table <- function(dge, pheno, df_ns, use_sex = FALSE, use_cells = FALSE, cellpc_cols = c("CellPC1","CellPC2")) {
  terms <- c(sprintf("ns(Age, df=%d)", df_ns), "StudyID")
  if (use_sex) terms <- c(terms, "Sex")
  if (use_cells) terms <- c(terms, paste0("scale(", cellpc_cols, ")"))
  design <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = pheno)

  if (!limma::is.fullrank(design)) {
    if (use_sex) {
      terms2 <- setdiff(terms, "Sex")
      design <- model.matrix(as.formula(paste("~", paste(terms2, collapse = " + "))), data = pheno)
    }
  }
  if (!limma::is.fullrank(design)) stop("Full-data design matrix is not full rank.")

  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  age_cols <- grep("^ns\\(Age", colnames(design))
  L <- matrix(0, nrow = ncol(design), ncol = length(age_cols))
  L[age_cols, ] <- diag(length(age_cols))
  rownames(L) <- colnames(design)
  fit2 <- contrasts.fit(fit, L)
  fit2 <- eBayes(fit2)
  tab <- topTable(fit2, number = Inf, sort.by = "F") %>% rownames_to_column("gene")
  tab
}

# -----------------------------
# Main function
# -----------------------------
run_stable_age_genes_final <- function(
    counts_file = "counts_mat.txt",
    pheno_file  = "pheno.txt",
    outdir      = "stable_age_outputs_final",

    # Resampling settings
    B           = 150,
    train_frac  = 0.8,
    strata_cols = c("AgeGroup", "StudyID"),

    # Candidate filtering (within each resample TRAIN)
    do_prefilter = TRUE,
    prefilter_top_prop = 0.10,
    prefilter_min_abs_rho = 0.20,
    prefilter_min_keep = 300,

    # Model settings
    df_grid     = c(2,3,4,5),
    alpha       = 0.05,
    use_sex     = FALSE,   # default FALSE due to missing Gender in pheno
    use_cells   = FALSE,
    cellpc_cols = c("CellPC1","CellPC2"),

    # Stable gene definition
    freq_thresh = 0.80,
    sign_consistency_thresh = 0.90,
    require_full_refit_fdr = 0.05,

    # Preprocessing
    min_count   = 10,
    seed        = 123
) {
  set.seed(seed)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # --- Read counts ---
  counts_df <- read.table(counts_file)
  counts_mat <- as.matrix(data.matrix(counts_df))

  # --- Read pheno and map columns from your pheno.txt ---
  pheno <- read.table(pheno_file,header = TRUE,
                      sep = "\t",
                      quote = "",
                      comment.char = "",
                      check.names = FALSE,
                      stringsAsFactors = FALSE) %>%
    rename(
      SampleID = `GSM.ID`,
      StudyID  = `GSE.ID`,
      Sex_raw  = Gender
    ) %>%
    mutate(
      Age    = as.numeric(Age),
      Group  = as.character(Group),
      StudyID = factor(StudyID),
      Sex_raw = ifelse(Sex_raw == "" | is.na(Sex_raw), NA_character_, Sex_raw),
      Sex    = factor(Sex_raw, levels = c("F","M"))
    )

  # Add AgeGroup for stratification
  pheno$AgeGroup <- factor(make_age_group(pheno$Age), levels = c("Young", "Mid", "Old"))

  # Drop rows with required covariates missing
  need_cols <- c("SampleID", "Age", "StudyID", "AgeGroup")
  if (use_sex) need_cols <- c(need_cols, "Sex")
  if (use_cells) need_cols <- c(need_cols, cellpc_cols)

  ok <- complete.cases(pheno[, need_cols, drop = FALSE])
  if (any(!ok)) {
    message("Dropping ", sum(!ok), " samples due to missing required covariates: ", paste(need_cols, collapse=", "))
    pheno <- pheno[ok, , drop = FALSE]
  }

  # --- Align ---
  aligned <- align_counts_pheno(counts_mat, pheno)
  counts_mat <- aligned$counts_mat
  pheno <- aligned$pheno

  # --- Build DGE and filter low-expression genes (shared baseline filter) ---
  dge <- DGEList(counts = counts_mat)

  # Design for filterByExpr (conservative; doesn't need ns)
  design0_terms <- c("AgeGroup", "StudyID")
  if (use_sex) design0_terms <- c(design0_terms, "Sex")
  design0 <- model.matrix(as.formula(paste("~", paste(design0_terms, collapse = " + "))), data = pheno)

  keep <- filterByExpr(dge, design = design0, min.count = min_count)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "TMM")

  message("After filterByExpr: genes=", nrow(dge), ", samples=", ncol(dge))

  # --- Evaluate df grid by stability ---
  df_summaries <- list()

  for (df_ns in df_grid) {
    message("\n=== Evaluating df = ", df_ns, " ===")

    selected_sets <- vector("list", B)
    sign_maps <- vector("list", B)

    for (b in seq_len(B)) {
      idx <- stratified_subsample(pheno, train_frac = train_frac, strata_cols = strata_cols)
      tr <- idx$train

      ph_tr <- pheno[tr, , drop = FALSE]
      dge_tr <- dge[, tr]

      # within-resample TRAIN prefilter
      if (do_prefilter) {
        keep_genes <- spearman_prefilter(
          dge_train = dge_tr,
          pheno_train = ph_tr,
          top_prop = prefilter_top_prop,
          min_abs_rho = prefilter_min_abs_rho,
          min_keep = prefilter_min_keep
        )
        dge_tr2 <- dge_tr[rownames(dge_tr) %in% keep_genes, ]
        if (nrow(dge_tr2) < prefilter_min_keep) dge_tr2 <- dge_tr
      } else {
        dge_tr2 <- dge_tr
      }

      fitres <- fit_select_ns_age(
        dge_train = dge_tr2,
        pheno_train = ph_tr,
        df_ns = df_ns,
        use_sex = use_sex,
        use_cells = use_cells,
        cellpc_cols = cellpc_cols,
        alpha = alpha
      )

      selected_sets[[b]] <- fitres$selected_genes
      sign_maps[[b]] <- fitres$sign_dir

      if (b %% 10 == 0) message("  resample ", b, "/", B, " selected=", length(fitres$selected_genes))
    }

    # selection frequency
    all_genes <- unique(unlist(selected_sets))
    if (length(all_genes) == 0) {
      warning("No genes selected for df=", df_ns, ". Consider relaxing alpha or prefilter thresholds.")
      next
    }
    freq <- sapply(all_genes, function(g) mean(vapply(selected_sets, function(S) g %in% S, logical(1))))
    freq <- sort(freq, decreasing = TRUE)

    # sign consistency among times selected (robust to empty sign_dir vectors)
    get_sign <- function(m, g) {
      if (is.null(m) || length(m) == 0) return(NA_real_)
      v <- m[g]                 
      if (length(v) == 0) return(NA_real_)
      as.numeric(v[[1]])
    }
    
    sign_cons <- sapply(names(freq), function(g) {
      s <- vapply(sign_maps, function(m) get_sign(m, g), numeric(1))
      s <- s[!is.na(s)]
      if (length(s) < 5) return(NA_real_)
      pmax(mean(s > 0), mean(s < 0))
    })
    

    stable_genes <- names(freq)[which(freq >= freq_thresh & (is.na(sign_cons) | sign_cons >= sign_consistency_thresh))]

    mean_jac <- mean_pairwise_jaccard(selected_sets, max_pairs = 2000)

    df_summary <- data.frame(
      df_ns = df_ns,
      B = B,
      mean_selected = mean(vapply(selected_sets, length, integer(1))),
      sd_selected   = sd(vapply(selected_sets, length, integer(1))),
      n_unique_selected = length(all_genes),
      mean_pairwise_jaccard = mean_jac,
      n_stable_genes = length(stable_genes),
      stringsAsFactors = FALSE
    )
    df_summaries[[as.character(df_ns)]] <- df_summary

    # Write per-df outputs
    out_freq <- tibble(
      gene = names(freq),
      selection_freq = as.numeric(freq),
      sign_consistency = as.numeric(sign_cons[names(freq)])
    ) %>% arrange(desc(selection_freq))

    write_csv(out_freq, file.path(outdir, sprintf("df%d_selection_frequency.csv", df_ns)))
    write_lines(stable_genes, file.path(outdir, sprintf("df%d_stable_genes.txt", df_ns)))

    # Save resample sets for reproducibility
    sets_lines <- vapply(seq_along(selected_sets), function(i) {
      paste0(i, "\t", paste(selected_sets[[i]], collapse = ","))
    }, character(1))
    write_lines(sets_lines, file.path(outdir, sprintf("df%d_selected_sets.tsv", df_ns)))
  }

  summary_df <- bind_rows(df_summaries) %>%
    arrange(desc(mean_pairwise_jaccard), desc(n_stable_genes), df_ns)

  write_csv(summary_df, file.path(outdir, "df_stability_summary.csv"))

  if (nrow(summary_df) == 0) stop("No df results were produced; please check inputs/thresholds.")

  # Choose df: highest stability; then simplest within delta
  best <- summary_df %>% arrange(desc(mean_pairwise_jaccard), desc(n_stable_genes), df_ns) %>% slice(1)
  best_df <- best$df_ns
  best_score <- best$mean_pairwise_jaccard
  delta <- 0.01
  chosen <- summary_df %>% filter(mean_pairwise_jaccard >= best_score - delta) %>% arrange(df_ns) %>% slice(1)
  chosen_df <- chosen$df_ns

  message("\n=== df selection ===")
  message("Best df by stability: df=", best_df, " (mean_jaccard=", round(best_score, 4), ")")
  message("Chosen simplest df within delta=", delta, ": df=", chosen_df)

  # Load frequency table for chosen df
  freq_tbl <- read_csv(file.path(outdir, sprintf("df%d_selection_frequency.csv", chosen_df)), show_col_types = FALSE)

  # Full-data refit for chosen df
  tab_full <- full_refit_table(dge, pheno, df_ns = chosen_df, use_sex = use_sex, use_cells = use_cells, cellpc_cols = cellpc_cols)
  write_csv(tab_full, file.path(outdir, "final_full_refit_table.csv"))

  # Merge and apply final stable definition (freq + sign + full-data FDR)
  final <- freq_tbl %>%
    left_join(tab_full %>% select(gene, F, P.Value, adj.P.Val), by = "gene") %>%
    mutate(
      pass_freq = selection_freq >= freq_thresh,
      pass_sign = is.na(sign_consistency) | sign_consistency >= sign_consistency_thresh,
      pass_full = is.na(adj.P.Val) | adj.P.Val < require_full_refit_fdr
    ) %>%
    filter(pass_freq, pass_sign, pass_full) %>%
    arrange(desc(selection_freq), adj.P.Val)

  write_csv(final, file.path(outdir, "final_stable_genes.csv"))
  write_lines(final$gene, file.path(outdir, "final_stable_genes.txt"))

  # Direction (full data Spearman) for interpretation
  expr_all <- cpm(dge, log = TRUE, prior.count = 1)
  rho_all <- apply(expr_all[final$gene, , drop = FALSE], 1,
                   function(x) suppressWarnings(cor(x, pheno$Age, method = "spearman")))
  dir_df <- tibble(
    gene = names(rho_all),
    spearman_rho_all = as.numeric(rho_all),
    direction = ifelse(rho_all > 0, "UpWithAge", "DownWithAge")
  ) %>% arrange(desc(abs(spearman_rho_all)))

  write_csv(dir_df, file.path(outdir, "final_stable_genes_direction.csv"))

  message("\nDone. Outputs in: ", normalizePath(outdir))
  invisible(list(
    df_summary = summary_df,
    chosen_df = chosen_df,
    final_stable = final
  ))
}

