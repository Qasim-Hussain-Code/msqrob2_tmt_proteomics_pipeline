## Print the numbers the README quotes, each from the results file that
## holds it, so that a reader (or the author) can regenerate them and
## check the README against them.
##
## Usage: scripts/lib/rscript.sh scripts/lib/readme_numbers.R

suppressPackageStartupMessages(library("data.table"))
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."), winslash = "/")
res <- function(f) fread(file.path(root, "results", f))
section <- function(x) cat("\n==== ", x, " ====\n", sep = "")

section("data_checksums.tsv, zenodo_record.tsv")
print(res("data_checksums.tsv")[, .(file, size_bytes, md5, md5_source)])
print(res("zenodo_record.tsv")[, .(record_id, revision, publication_date, license, n_files, total_bytes)])

section("psm_counts_raw_all.tsv")
print(res("psm_counts_raw_all.tsv")[, .(runs = .N, psms = sum(n_psms)), by = dataset])

section("filter logs")
for (ds in c("spikein1", "spikein1_ms2", "mouse")) print(res(paste0(ds, "_filter_log.tsv"))[, .(dataset, step, removed, psms_after)])

section("normalisation_factors: within-run range of channel medians")
for (ds in c("spikein1", "spikein1_ms2", "mouse")) {
    x <- res(paste0(ds, "_normalisation_factors.tsv"))
    r <- x[, diff(range(median_log2_before)), by = run]$V1
    cat(sprintf("%s: median %.3f, max %.3f\n", ds, median(r), max(r)))
}

section("missingness after filters")
for (ds in c("spikein1", "spikein1_ms2", "mouse")) {
    x <- res(paste0(ds, "_missingness.tsv"))
    cat(sprintf("%s: fraction of channel values NA after filters, min %.4f max %.4f\n", ds, min(x$frac_values_na_after), max(x$frac_values_na_after)))
}

section("mds")
for (ds in c("spikein1", "mouse")) { x <- res(paste0(ds, "_mds.tsv")); cat(sprintf("%s: %d complete proteins, %d channels, dim1 %.3f dim2 %.3f\n", ds, x$n_complete_proteins[1], nrow(x), x$var_explained_dim1[1], x$var_explained_dim2[1])) }

section("spikein1_benchmark_metrics.tsv at BH 5%, SPS-MS3, balanced")
m <- res("spikein1_benchmark_metrics.tsv")[alpha == 0.05]
m[, contrast := factor(contrast, levels = c("0.667 - 0.5", "1 - 0.667", "1 - 0.5", "0.5 - 0.125", "0.667 - 0.125", "1 - 0.125"))]
show <- function(d) print(d[, .(variant, contrast, n_tested, tp, fp, fdp = round(fdp, 3), sens = round(sensitivity, 3), bias = round(bias, 3))][order(variant, contrast)])
show(m[acquisition == "ms3" & design == "balanced"])
section("unbalanced")
show(m[design == "unbalanced"])
section("ms2")
show(m[acquisition == "ms2"])

section("spikein1_naive_vs_mixed.tsv")
print(res("spikein1_naive_vs_mixed.tsv")[, .(design, contrast, sig_mixed, fp_mixed, fdp_mixed = round(fdp_mixed, 3), sig_naive, fp_naive, fdp_naive = round(fdp_naive, 3), se_ratio = round(median_se_ratio_naive_over_mixed, 2))])

section("spikein1_benchmark_variants.tsv")
print(res("spikein1_benchmark_variants.tsv")[, .(acquisition, design, variant, n_fit_error, fit_seconds, peak_rss_mb, note)])
if (file.exists(file.path(root, "results", "spikein1_psm_subset.tsv"))) print(res("spikein1_psm_subset.tsv"))

section("spikein1_false_positive_summary.tsv, adopted variants")
print(res("spikein1_false_positive_summary.tsv")[variant %in% c("protein_drop_mixed_robust_noridge_msqrob2", "protein_drop_mixed_ls_noridge_msqrob2", "protein_drop_naive_robust_noridge_msqrob2") & acquisition == "ms3",
      .(design, variant, contrast, n_fp, median_logfc = round(median_logfc, 3), median_abs = round(median_abs_logfc, 3), frac_negative = round(frac_negative, 2), frac_abs_below_0.2 = round(frac_abs_below_0.2, 2))])

section("spikein1_fdp_by_effect_size_floor.tsv, mixed robust, balanced ms3")
print(res("spikein1_fdp_by_effect_size_floor.tsv")[variant == "protein_drop_mixed_robust_noridge_msqrob2" & acquisition == "ms3" & design == "balanced", .(contrast, logfc_floor, tp, fp, fdp = round(fdp, 3))][order(contrast, logfc_floor)])

section("spikein1_lmer_diagnostics_summary.tsv")
print(res("spikein1_lmer_diagnostics_summary.tsv"))

section("compression")
print(res("spikein1_compression_by_contrast.tsv")[, .(acquisition, contrast, expected_log2fc = round(expected_log2fc, 3), ups_median = round(ups_median, 3), ups_iqr = round(ups_iqr, 3), hela_median = round(hela_median, 3), observed_over_expected = round(observed_over_expected, 3))])
print(res("spikein1_compression_slope.tsv")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
print(res("spikein1_compression_by_interference.tsv")[is_ups == TRUE, .(acquisition, interference_bin, n_psms, median_log2fc = round(median_log2fc, 3))])
print(res("spikein1_interference_summary.tsv")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

section("mouse")
print(dcast(res("mouse_significant_counts.tsv")[alpha == 0.05], contrast ~ model, value.var = "n_significant"))
print(dcast(res("mouse_significant_counts.tsv")[alpha == 0.10], contrast ~ model, value.var = "n_significant"))
print(res("mouse_models.tsv")[, .(model, n_proteins, n_fit_error, n_refitted, fit_seconds)])
print(res("mouse_fit_tiers.tsv")[, .(model, tier, attempted, fitted)])
print(res("mouse_lmer_error_causes.tsv"))
print(res("mouse_lmer_diagnostics_summary.tsv")[, .(scope, level, n_proteins, n_error, n_singular, frac_singular = round(frac_singular, 3), median_var_Mixture = signif(median_var_Mixture, 3), median_var_Run = signif(median_var_Run, 3), median_var_BioReplicate = signif(median_var_BioReplicate, 3), median_sigma2 = signif(median_sigma2, 3), frac_zero_var_Mixture = round(frac_zero_var_Mixture, 3))])
print(res("mouse_scope_comparison.tsv")[, .(contrast, n_both_tested, sig_fraction_run, sig_mixture, sig_both, logfc_correlation = round(logfc_correlation, 3), se_ratio = round(median_se_ratio_mixture_over_fraction, 3))])
print(res("mouse_naive_vs_mixed.tsv")[, .(scope, contrast, sig_mixed, sig_naive, sig_both, sig_naive_only, se_ratio = round(median_se_ratio_naive_over_mixed, 2))])
r <- res("mouse_results.tsv.gz")
cat("8wk calls by tier (fraction_run_mixed):\n"); print(table(r[model == "fraction_run_mixed" & contrast == "HF_vs_LF_8wk" & adjPval < 0.05, fit_tier]))
cat("8wk direction (fraction_run_mixed): up", r[model == "fraction_run_mixed" & contrast == "HF_vs_LF_8wk" & adjPval < 0.05 & logFC > 0, .N], "down", r[model == "fraction_run_mixed" & contrast == "HF_vs_LF_8wk" & adjPval < 0.05 & logFC < 0, .N], "\n")
w <- dcast(r[level == "protein" & scope == "fraction_run" & contrast == "HF_vs_LF_18wk"], protein ~ structure, value.var = c("logFC", "se", "adjPval"))
no <- w[adjPval_naive < 0.05 & !(adjPval_mixed < 0.05)]
cat(sprintf("naive-only at 18wk: %d; median |logFC diff| %.2f; median se ratio naive/mixed %.2f; range %.2f to %.2f\n", nrow(no), median(abs(no$logFC_naive - no$logFC_mixed)), median(no$se_naive / no$se_mixed), min(no$se_naive / no$se_mixed), max(no$se_naive / no$se_mixed)))
cat("interaction top 8 (fraction_run_mixed):\n"); print(r[model == "fraction_run_mixed" & contrast == "interaction"][order(pval)][1:8, .(protein, logFC = round(logFC, 2), adjPval = signif(adjPval, 2), fit_tier)])
print(res("mouse_bcaa_kegg_enrichment.tsv")[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
print(res("mouse_bcaa_kegg_by_model.tsv")[, .(model, n_tested, n_pathway_sig_10, n_pathway_sig_5, n_pathway_negative, fisher_p_10 = signif(fisher_p_10, 2))])
print(res("mouse_bcaa_kegg_members.tsv")[, .(gene, logFC = round(logFC, 2), adjPval = signif(adjPval, 2))])
print(res("mouse_bcaa_keyword_check.tsv")[, -"terms"])
print(res("mouse_design_counts_used.tsv"))

section("resources")
print(fread(file.path(root, "logs", "resources_all.tsv")))
