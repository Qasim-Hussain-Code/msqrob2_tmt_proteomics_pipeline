## Stage 08: the publication figures, all from the result tables written
## by stages 04 to 07. Nothing is recomputed here; a figure that needs a
## number the tables do not hold gets the number added to the table in
## the stage that owns it.
##
## Usage: scripts/lib/rscript.sh scripts/08_figures.R [--force]
##
## Figures:
##   figures/fig_benchmark_fdp.png          FDP and sensitivity per contrast, main workflow variants
##   figures/fig_naive_vs_mixed.png         balanced against unbalanced design, naive against mixed
##   figures/fig_reference_treatment.png    the three reference channel treatments
##   figures/fig_variance_components.png    where the variance sits, spike-in and mouse
##   figures/fig_mouse_counts.png           significant proteins per contrast, scope by structure
##   figures/fig_mouse_interaction.png      the interaction result, protein by protein
##   figures/fig_mouse_bcaa.png             branched chain amino acid catabolism proteins at 8 weeks

suppressPackageStartupMessages({
    library("data.table")
    library("ggplot2")
    library("patchwork")
    library("ggrepel")
    library("QFeatures")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
source(file.path(root, "scripts", "lib", "common.R"))
conf <- read_conf(root)
opt <- parse_args()
STAGE <- "08_figures"
if (stage_already_done(STAGE, conf, force = opt$force)) quit(save = "no")
timer <- stage_begin(STAGE, conf)

res_dir <- file.path(root, "results")
fig_dir <- file.path(root, "figures")
expected <- spikein_expected_contrasts()
contrast_levels <- expected[order(expected_log2fc), contrast]
theme_set(theme_bw(base_size = 10))

## ---- Benchmark: FDP and sensitivity ------------------------------------------

metrics <- fread(file.path(res_dir, "spikein1_benchmark_metrics.tsv"))
metrics[, contrast := factor(contrast, levels = contrast_levels)]
label_variant <- function(m) {
    m[, workflow := fifelse(engine == "limma", paste0("limma dupcor, ref ", reference),
                    fifelse(level == "psm", paste0("PSM mixed, ref ", reference, fifelse(robust, ", robust", ", LS"), fifelse(ridge, ", ridge", "")),
                            paste0("protein ", structure, ", ref ", reference, fifelse(robust, ", robust", ", LS"), fifelse(ridge, ", ridge", ""))))]
    m
}
metrics <- label_variant(metrics)

main <- metrics[acquisition == "ms3" & design == "balanced" & alpha == 0.05 & reference == "drop" & !ridge &
                ((level == "protein" & robust & engine == "msqrob2") | engine == "limma" | (level == "psm" & robust))]
main[, workflow := factor(workflow, levels = unique(workflow[order(level != "protein", structure)]))]
p_fdp <- ggplot(main, aes(x = contrast, y = fdp, colour = workflow, group = workflow)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
    geom_line(alpha = 0.6) + geom_point(size = 2.2) +
    labs(x = NULL, y = "realised false discovery proportion at BH 5%", colour = NULL,
         title = "Spike-in, balanced design: what each workflow's 5% FDR list contains") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2))
p_sens <- ggplot(main, aes(x = contrast, y = sensitivity, colour = workflow, group = workflow)) +
    geom_line(alpha = 0.6) + geom_point(size = 2.2) +
    labs(x = "contrast (UPS1 fmol ratio)", y = "sensitivity (UPS1 proteins called)", colour = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
save_fig(p_fdp / p_sens + plot_layout(heights = c(1.2, 1)), file.path(fig_dir, "fig_benchmark_fdp.png"), 8.5, 8)

## ---- Naive against mixed, balanced against unbalanced ---------------------------

h2h <- fread(file.path(res_dir, "spikein1_naive_vs_mixed.tsv"))
h2h[, contrast := factor(contrast, levels = contrast_levels)]
h2h_long <- melt(h2h[, .(design, contrast, mixed = fdp_mixed, naive = fdp_naive)], id.vars = c("design", "contrast"),
                 variable.name = "model", value.name = "fdp")
h2h_sig <- melt(h2h[, .(design, contrast, mixed = sig_mixed, naive = sig_naive)], id.vars = c("design", "contrast"),
                variable.name = "model", value.name = "n_significant")
p_h1 <- ggplot(h2h_long, aes(x = contrast, y = fdp, fill = model)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    facet_wrap(~ design) +
    scale_fill_manual(values = c(mixed = "#2b8cbe", naive = "#d95f0e")) +
    labs(x = NULL, y = "realised FDP at BH 5%", title = "Ignoring the plex: naive against mixed, same proteins, same data") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
p_h2 <- ggplot(h2h_sig, aes(x = contrast, y = n_significant, fill = model)) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    facet_wrap(~ design) +
    scale_fill_manual(values = c(mixed = "#2b8cbe", naive = "#d95f0e")) +
    labs(x = "contrast", y = "proteins at BH 5%") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
p_h3 <- ggplot(h2h, aes(x = contrast, y = median_se_ratio_naive_over_mixed, colour = design, group = design)) +
    geom_hline(yintercept = 1, colour = "grey50") +
    geom_line() + geom_point(size = 2) +
    labs(x = "contrast", y = "median SE ratio, naive / mixed") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
save_fig((p_h1 / p_h2 / p_h3) + plot_layout(heights = c(1.3, 1, 0.9)), file.path(fig_dir, "fig_naive_vs_mixed.png"), 8.5, 10)

## ---- Reference channel treatment ---------------------------------------------

ref <- metrics[acquisition == "ms3" & design == "balanced" & alpha == 0.05 & level == "protein" & robust & !ridge & engine == "msqrob2" & structure == "mixed"]
p_r1 <- ggplot(ref, aes(x = contrast, y = fdp, colour = reference, group = reference)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
    geom_line(alpha = 0.6) + geom_point(size = 2.2) +
    labs(x = NULL, y = "realised FDP at BH 5%", title = "Reference channels: drop, keep as a condition, or divide by them (protein mixed model, robust)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
p_r2 <- ggplot(ref, aes(x = contrast, y = sensitivity, colour = reference, group = reference)) +
    geom_line(alpha = 0.6) + geom_point(size = 2.2) +
    labs(x = "contrast", y = "sensitivity") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
p_r3 <- ggplot(ref, aes(x = contrast, y = bias, colour = reference, group = reference)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_line(alpha = 0.6) + geom_point(size = 2.2) +
    labs(x = "contrast", y = "median UPS1 log2FC minus expected") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
save_fig(p_r1 / (p_r2 | p_r3), file.path(fig_dir, "fig_reference_treatment.png"), 8.5, 7.5)

## ---- Variance components, spike-in and mouse ---------------------------------

vc_files <- list(spikein = "spikein1_lmer_diagnostics.tsv", mouse = "mouse_lmer_diagnostics.tsv")
vc <- rbindlist(lapply(names(vc_files), function(n) {
    d <- fread(file.path(res_dir, vc_files[[n]]))
    if (n == "spikein") d <- d[acquisition == "ms3" & design == "balanced" & level == "protein" & reference == "drop"]
    if (n == "mouse") d <- d[level == "protein"]
    d[, dataset := if (n == "spikein") "spike-in, SPS-MS3, protein level" else paste0("mouse, ", scope, " scope")]
    cols <- grep("^var_|^sigma2$", names(d), value = TRUE)
    melt(d[status != "error", c("protein", "dataset", cols), with = FALSE], id.vars = c("protein", "dataset"),
         variable.name = "component", value.name = "variance")
}), fill = TRUE)[!is.na(variance)]
vc[, component := sub("^var_", "", component)]
vc[component == "sigma2", component := "residual"]
vc[, component := factor(component, levels = c("Mixture", "Run", "BioReplicate", "residual"))]
frac_zero <- vc[, .(frac_zero = mean(variance == 0), n = .N), by = .(dataset, component)]
p_vc <- ggplot(vc, aes(x = component, y = variance + 1e-4)) +
    geom_violin(fill = "grey90", scale = "width") +
    geom_boxplot(width = 0.14, outlier.size = 0.3) +
    geom_text(data = frac_zero, aes(x = component, y = 3, label = sprintf("%.0f%% at zero", 100 * frac_zero)), size = 2.8) +
    scale_y_log10() + facet_wrap(~ dataset, scales = "free_x") +
    labs(x = NULL, y = "variance component per protein (log10, + 1e-4)",
         title = "Variance components from lme4: a component at the floor is a singular fit for that term") +
    theme(legend.position = "none")
save_fig(p_vc, file.path(fig_dir, "fig_variance_components.png"), 10, 4.4)

## ---- Mouse: counts, interaction, BCAA -------------------------------------------

mcounts <- fread(file.path(res_dir, "mouse_significant_counts.tsv"))
mlabels <- fread(file.path(root, "config", "contrasts_mouse.tsv"))$contrast
mcounts[, contrast := factor(contrast, levels = mlabels)]
mcounts[, model := factor(model, levels = c("fraction_run_mixed", "fraction_run_naive", "mixture_mixed", "mixture_naive", "psm_mixed"))]
p_mc <- ggplot(mcounts[alpha == 0.05], aes(x = contrast, y = n_significant, fill = model)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = n_significant), position = position_dodge(0.8), vjust = -0.3, size = 2.6) +
    scale_fill_manual(values = c(fraction_run_mixed = "#2b8cbe", fraction_run_naive = "#d95f0e", mixture_mixed = "#093c57",
                                 mixture_naive = "#fdae6b", psm_mixed = "#7bccc4")) +
    labs(x = NULL, y = "proteins at BH 5%", fill = NULL,
         title = "Mouse adipose: significant proteins per contrast under each model") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1)) +
    guides(fill = guide_legend(nrow = 2))
save_fig(p_mc, file.path(fig_dir, "fig_mouse_counts.png"), 8, 5)

## Interaction: the top proteins by interaction p-value, shown as cell means
## of the mixture-scope values (one per mouse), with each mouse a point.
mres <- fread(file.path(res_dir, "mouse_results.tsv.gz"))
inter <- mres[model == "fraction_run_mixed" & contrast == "interaction" & !is.na(pval)][order(pval)]
top <- head(inter$protein, 6)
qfm <- readRDS(file.path(root, "results", "qfeatures", "mouse_mixture_preprocessed.rds"))
cdm <- as.data.table(as.data.frame(colData(qfm)), keep.rownames = "sample")[Condition != "Norm" & Condition != "Long_M"]
m <- assay(qfm[["proteins"]])[top, cdm$sample, drop = FALSE]
long <- melt(as.data.table(m, keep.rownames = "protein"), id.vars = "protein", variable.name = "sample", value.name = "value")
long <- merge(long[!is.na(value)], cdm, by = "sample")
desc <- as.data.table(as.data.frame(rowData(readRDS(file.path(root, "results", "qfeatures", "mouse_preprocessed.rds"))[["proteins"]])[, "Master.Protein.Descriptions", drop = FALSE]), keep.rownames = "protein")
desc[, gene := sub(".*GN=([^ ]+).*", "\\1", Master.Protein.Descriptions)]
long <- merge(long, desc[, .(protein, gene)], by = "protein")
long <- merge(long, inter[, .(protein, interaction_adjp = adjPval, interaction_logfc = logFC)], by = "protein")
long[, panel := sprintf("%s (%s)\ninteraction log2FC %.2f, adj. p %.2g", gene, protein, interaction_logfc, interaction_adjp)]
long[, panel := factor(panel, levels = unique(panel[order(interaction_adjp)]))]
long[, Duration := factor(Duration, levels = c("Short", "Long"), labels = c("8 weeks", "18 weeks"))]
long[, Diet := factor(Diet, levels = c("LF", "HF"))]
p_int <- ggplot(long, aes(x = Duration, y = value, colour = Diet, shape = mixture_short)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), size = 2) +
    stat_summary(aes(group = Diet), fun = mean, geom = "line", position = position_dodge(0.5), linewidth = 0.5) +
    stat_summary(aes(group = Diet), fun = mean, geom = "point", position = position_dodge(0.5), size = 3, shape = 21, fill = "white") +
    facet_wrap(~ panel, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c(LF = "#1b9e77", HF = "#d95f02")) +
    labs(x = NULL, y = "mixture-scope log2 protein value (centered within mixture)", shape = "mixture",
         title = "Diet by duration interaction: the six smallest interaction p-values, one point per mouse") +
    theme(legend.position = "bottom")
save_fig(p_int, file.path(fig_dir, "fig_mouse_interaction.png"), 10, 7)

## BCAA catabolism proteins at 8 weeks: every KEGG mmu00280 member tested.
kegg_file <- file.path(res_dir, "mouse_bcaa_kegg_members.tsv")
if (file.exists(kegg_file)) {
    km <- fread(kegg_file)
    km[, gene := factor(gene, levels = gene[order(logFC)])]
    km[, significance := fifelse(adjPval < 0.05, "BH < 0.05", fifelse(adjPval < 0.10, "BH < 0.10", "not significant"))]
    p_bcaa <- ggplot(km, aes(x = logFC, y = gene, colour = significance)) +
        geom_vline(xintercept = 0, colour = "grey50") +
        geom_point(size = 2.2) +
        scale_colour_manual(values = c(`BH < 0.05` = "firebrick", `BH < 0.10` = "darkorange", `not significant` = "grey50")) +
        labs(x = "log2 fold change, HF against LF at 8 weeks (fraction-run mixed model)", y = NULL, colour = NULL,
             title = sprintf("KEGG mmu00280 (valine, leucine and isoleucine degradation): %d of its proteins were tested", nrow(km))) +
        theme(legend.position = "bottom")
    save_fig(p_bcaa, file.path(fig_dir, "fig_mouse_bcaa.png"), 7.5, 0.18 * nrow(km) + 2)
}

stage_end(timer)
stage_mark_done(STAGE, conf)
