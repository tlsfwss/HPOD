#!/usr/bin/env Rscript
# ============================================================
# Visualization script for stable_age_outputs_final
# Reads outputs produced by stable_age_genes_final.R and generates
# publication-friendly plots (PDF + PNG):
#   1) df stability (mean Jaccard vs df, stable gene count)
#   2) selection frequency rank plot (chosen df)
#   3) sign consistency vs frequency (chosen df)
#   4) expression vs age trend plots for top stable genes
#
# Usage (R):
#   source("stable_age_genes_plots.R")
#   make_stability_plots(
#       outdir="stable_age_outputs_final",
#       counts_file="counts_mat.txt",
#       pheno_file="pheno.txt",
#       top_n_trends=12
#   )
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# Optional plotting deps
.has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

# --------- helpers ----------
.make_dir <- function(x) dir.create(x, showWarnings = FALSE, recursive = TRUE)

.read_df_safe <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

# Robust parse chosen df from available files
.infer_chosen_df <- function(outdir) {
  # If final_stable_genes exists, infer from df summary: choose simplest within delta of best
  summ <- .read_df_safe(file.path(outdir, "df_stability_summary.csv"))
  best <- summ %>% arrange(desc(mean_pairwise_jaccard), desc(n_stable_genes), df_ns) %>% slice(1)
  delta <- 0.01
  cand <- summ %>% filter(mean_pairwise_jaccard >= best$mean_pairwise_jaccard - delta) %>% arrange(df_ns) %>% slice(1)
  cand$df_ns
}

# Load pheno as used in final script (mapping)
.load_pheno_mapped <- function(pheno_file) {
  ph <- readr::read_tsv(pheno_file, show_col_types = FALSE) %>%
    rename(
      SampleID = `GSM.ID`,
      StudyID  = `GSE.ID`,
      Sex_raw  = Gender
    ) %>%
    mutate(
      Age = as.numeric(Age),
      Group = as.character(Group),
      StudyID = factor(StudyID),
      Sex_raw = ifelse(Sex_raw == "" | is.na(Sex_raw), NA_character_, Sex_raw),
      Sex = factor(Sex_raw, levels = c("F","M"))
    ) %>%
    filter(Group %in% c("healthy", "normal")) %>%
    mutate(AgeGroup = cut(Age, breaks=c(-Inf,39,59,Inf), labels=c("Young","Mid","Old"), right=TRUE))
  ph
}

.load_counts <- function(counts_file) {
  counts_df <- readr::read_tsv(counts_file, show_col_types = FALSE)
  gene_col <- colnames(counts_df)[1]
  gene_ids <- counts_df[[gene_col]]
  mat <- as.matrix(data.matrix(counts_df[, -1, drop = FALSE]))
  rownames(mat) <- gene_ids
  mat
}

# --------- main plotting function ----------
make_stability_plots <- function(outdir = "stable_age_outputs_final",
                                 counts_file = "counts_mat.txt",
                                 pheno_file  = "pheno.txt",
                                 plots_dir   = file.path(outdir, "plots"),
                                 top_n_trends = 12,
                                 point_alpha = 0.5) {

  .make_dir(plots_dir)

  summ <- .read_df_safe(file.path(outdir, "df_stability_summary.csv"))
  chosen_df <- .infer_chosen_df(outdir)

  freq_path <- file.path(outdir, sprintf("df%d_selection_frequency.csv", chosen_df))
  freq <- .read_df_safe(freq_path)

  final_path <- file.path(outdir, "final_stable_genes.csv")
  final <- .read_df_safe(final_path)

  # ---------- Plot 1: df stability ----------
  if (.has_ggplot2) {
    library(ggplot2)

    p1 <- ggplot(summ, aes(x = df_ns, y = mean_pairwise_jaccard)) +
      geom_line() + geom_point() +
      labs(x="Spline df", y="Mean pairwise Jaccard", title="Stability vs spline df") +
      theme_bw()

    p1b <- ggplot(summ, aes(x = df_ns, y = n_stable_genes)) +
      geom_line() + geom_point() +
      labs(x="Spline df", y=paste0("Stable genes (freq threshold)"),
           title="Stable gene count vs spline df") +
      theme_bw()

    ggsave(file.path(plots_dir, "df_stability_jaccard.pdf"), p1, width=6, height=4)
    ggsave(file.path(plots_dir, "df_stability_jaccard.png"), p1, width=6, height=4, dpi=300)
    ggsave(file.path(plots_dir, "df_stability_stablecount.pdf"), p1b, width=6, height=4)
    ggsave(file.path(plots_dir, "df_stability_stablecount.png"), p1b, width=6, height=4, dpi=300)

    # ---------- Plot 2: selection frequency rank plot ----------
    freq_rank <- freq %>% arrange(desc(selection_freq)) %>% mutate(rank = row_number())
    p2 <- ggplot(freq_rank, aes(x=rank, y=selection_freq)) +
      geom_point(alpha=0.6) +
      geom_hline(yintercept = 0.8, linetype="dashed") +
      labs(x="Gene rank by selection frequency", y="Selection frequency",
           title=paste0("Selection frequency rank plot (df=", chosen_df, ")")) +
      theme_bw()

    ggsave(file.path(plots_dir, sprintf("df%d_selection_frequency_rank.pdf", chosen_df)), p2, width=6, height=4)
    ggsave(file.path(plots_dir, sprintf("df%d_selection_frequency_rank.png", chosen_df)), p2, width=6, height=4, dpi=300)

    # ---------- Plot 3: sign consistency vs frequency ----------
    p3 <- ggplot(freq, aes(x=selection_freq, y=sign_consistency)) +
      geom_point(alpha=0.6) +
      labs(x="Selection frequency", y="Sign consistency",
           title=paste0("Sign consistency vs frequency (df=", chosen_df, ")")) +
      theme_bw()

    ggsave(file.path(plots_dir, sprintf("df%d_signcons_vs_freq.pdf", chosen_df)), p3, width=6, height=4)
    ggsave(file.path(plots_dir, sprintf("df%d_signcons_vs_freq.png", chosen_df)), p3, width=6, height=4, dpi=300)

  } else {
    # base R fallback
    pdf(file.path(plots_dir, "df_stability_jaccard.pdf"), width=6, height=4)
    plot(summ$df_ns, summ$mean_pairwise_jaccard, type="b", xlab="Spline df", ylab="Mean pairwise Jaccard",
         main="Stability vs spline df")
    dev.off()

    pdf(file.path(plots_dir, "df_stability_stablecount.pdf"), width=6, height=4)
    plot(summ$df_ns, summ$n_stable_genes, type="b", xlab="Spline df", ylab="Stable gene count",
         main="Stable gene count vs spline df")
    dev.off()

    freq_rank <- freq[order(-freq$selection_freq),]
    pdf(file.path(plots_dir, sprintf("df%d_selection_frequency_rank.pdf", chosen_df)), width=6, height=4)
    plot(seq_len(nrow(freq_rank)), freq_rank$selection_freq, pch=16, cex=0.6,
         xlab="Gene rank", ylab="Selection frequency",
         main=paste0("Selection frequency rank plot (df=", chosen_df, ")"))
    abline(h=0.8, lty=2)
    dev.off()
  }

  # ---------- Plot 4: expression vs age trend plots for top stable genes ----------
  # Use CPM(log) for visualization
  ph <- .load_pheno_mapped(pheno_file)
  counts <- .load_counts(counts_file)

  # align
  common <- intersect(colnames(counts), ph$SampleID)
  counts <- counts[, common, drop=FALSE]
  ph <- ph[match(common, ph$SampleID), , drop=FALSE]

  # compute logCPM with edgeR 
  suppressPackageStartupMessages(library(edgeR))
  dge <- DGEList(counts=counts)
  dge <- calcNormFactors(dge, method="TMM")
  logcpm <- cpm(dge, log=TRUE, prior.count=1)

  top_genes <- final$gene[1:min(top_n_trends, nrow(final))]
  top_genes <- top_genes[top_genes %in% rownames(logcpm)]
  if (length(top_genes) == 0) {
    warning("No top genes found in logCPM matrix; skipping trend plots.")
    return(invisible(list(chosen_df=chosen_df)))
  }

  if (.has_ggplot2) {
    library(ggplot2)
    trend_df <- lapply(top_genes, function(g) {
      tibble(gene=g, Age=ph$Age, expr=as.numeric(logcpm[g, ]), StudyID=ph$StudyID)
    }) %>% bind_rows()

    p4 <- ggplot(trend_df, aes(x=Age, y=expr)) +
      geom_point(alpha=point_alpha) +
      geom_smooth(method="loess", se=TRUE) +
      facet_wrap(~gene, scales="free_y") +
      labs(x="Age", y="log2-CPM (TMM)", title="Expression trends of top stable genes") +
      theme_bw()

    ggsave(file.path(plots_dir, "top_stable_genes_trends.pdf"), p4, width=10, height=7)
    ggsave(file.path(plots_dir, "top_stable_genes_trends.png"), p4, width=10, height=7, dpi=300)

  } else {
    pdf(file.path(plots_dir, "top_stable_genes_trends.pdf"), width=10, height=7)
    par(mfrow=c(ceiling(length(top_genes)/4), 4), mar=c(3,3,2,1))
    for (g in top_genes) {
      plot(ph$Age, as.numeric(logcpm[g, ]), pch=16, cex=0.6,
           xlab="Age", ylab="log2-CPM", main=g)
      lines(lowess(ph$Age, as.numeric(logcpm[g, ])), lwd=2)
    }
    dev.off()
  }

  message("Plots saved to: ", normalizePath(plots_dir))
  invisible(list(chosen_df=chosen_df, plots_dir=plots_dir))
}


