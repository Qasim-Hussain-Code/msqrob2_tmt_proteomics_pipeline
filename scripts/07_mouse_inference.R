## Stage 07: does a high fat diet change the mouse adipose proteome, and
## does the answer depend on how long the diet ran?
##
## Twenty mice in a diet by duration factorial, five per cell, spread
## unevenly over three TMT 10-plex mixtures, each fractionated into nine
## runs. The fixed effects are diet, duration and their interaction; the
## random effects are the structure the msqrob2TMT mouse vignette uses
## for this design, restricted to what survives summarisation.
##
## Model formulas, stated here, in the fit calls below and in the README:
##   fraction-run scope, mixed : ~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | BioReplicate)
##   mixture scope, mixed      : ~ Diet * Duration + (1 | Mixture)
##   naive (either scope)      : ~ Diet * Duration
##   PSM level, mixed          : ~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | Run:Channel) + (1 | Run:ionID) + (1 | BioReplicate)
##
## Proteins the full model cannot fit are refitted with the next simpler
## random structure (fit_tiered in scripts/lib/msqrob_helpers.R), the
## tier is recorded per protein, and the counts per tier are reported.
## The fixed effects never change between tiers.
##
## In the fraction-run scope each mouse contributes nine protein values
## (one per fraction), so the mouse random effect is identifiable and
## Run carries the fraction. In the mixture scope each mouse contributes
## one value, so a mouse effect would be the residual and only the
## mixture term remains. The naive model treats every value as an
## independent sample; in the fraction-run scope that means it counts
## nine values per mouse as nine mice.
##
## Reference levels: Diet = LF, Duration = Short (8 weeks). DietHF is
## then the high fat effect at 8 weeks and DietHF:DurationLong is how
## much that effect changes by 18 weeks (config/contrasts_mouse.tsv).
##
## Usage: scripts/lib/rscript.sh scripts/07_mouse_inference.R
##          [--reference drop|keep|ratio] [--psm-level yes|no] [--workers N] [--force]
##
## Outputs:
##   results/mouse_results.tsv.gz              every protein, contrast and model
##   results/mouse_significant_counts.tsv      counts at three thresholds
##   results/mouse_scope_comparison.tsv        fraction-run against mixture scope, per contrast
##   results/mouse_naive_vs_mixed.tsv          naive against mixed, per scope and contrast
##   results/mouse_lmer_diagnostics.tsv        singular fits and variance components
##   results/mouse_lmer_diagnostics_summary.tsv
##   results/mouse_bcaa_keyword_check.tsv      keyword hits on protein descriptions, 8 week contrast
##   results/mouse_bcaa_kegg_enrichment.tsv    KEGG mmu00280 against the tested background
##   figures/mouse_design_matrix.png, mouse_volcano.png, mouse_scope_comparison.png, mouse_variance_components.png

suppressPackageStartupMessages({
    library("data.table")
    library("QFeatures")
    library("msqrob2")
    library("BiocParallel")
    library("ggplot2")
    library("ggrepel")
    library("patchwork")
    library("ExploreModelMatrix")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
source(file.path(root, "scripts", "lib", "common.R"))
source(file.path(root, "scripts", "lib", "reference_channels.R"))
source(file.path(root, "scripts", "lib", "msqrob_helpers.R"))
conf <- read_conf(root)
opt <- parse_args()
STAGE <- "07_mouse_inference"
if (stage_already_done(STAGE, conf, force = opt$force)) quit(save = "no")
timer <- stage_begin(STAGE, conf)
set.seed(20260921)
register_parallel(if (is.null(opt$workers)) conf$WORKERS else opt$workers)
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

## The reference channel treatment adopted from the spike-in benchmark.
## The default here is what 05 supported; see the README for the numbers.
reference <- if (is.null(opt$reference)) "drop" else opt$reference
run_psm <- is.null(opt$`psm-level`) || opt$`psm-level` != "no"

res_dir <- file.path(root, "results")
fig_dir <- file.path(root, "figures")
qf_dir <- file.path(root, "results", "qfeatures")

## ---- Data ---------------------------------------------------------------

contr <- fread(file.path(root, "config", "contrasts_mouse.tsv"))
hyp <- contr$hypothesis
labels <- contr$contrast
params <- c("DietHF", "DurationLong", "DietHF:DurationLong")

prep_mouse <- function(qf, sets) {
    ## The four Long_M channels hold the additional aged mice of the
    ## published study; they are outside the two by two and are removed
    ## before the reference treatment.
    qf <- subsetByColData(qf, qf$Condition != "Long_M")
    qf <- apply_reference_treatment(qf, mode = reference, sets = sets, run_col = "Run")$qf
    cd <- colData(qf)
    cd$Diet <- factor(cd$Diet, levels = c("LF", "HF"))
    cd$Duration <- factor(cd$Duration, levels = c("Short", "Long"))
    colData(qf) <- cd
    ## The assay-level colData must not repeat these columns as
    ## character, or getWithColData refuses the factor versions.
    for (i in sets) {
        se <- qf[[i]]
        colData(se) <- colData(se)[, setdiff(colnames(colData(se)), colnames(cd)), drop = FALSE]
        qf <- replaceAssay(qf, se, i)
    }
    qf
}

qf_run <- prep_mouse(readRDS(file.path(qf_dir, "mouse_preprocessed.rds")), c("ions_norm", "proteins"))
qf_mix0 <- readRDS(file.path(qf_dir, "mouse_mixture_preprocessed.rds"))
colData(qf_mix0)$Run <- colData(qf_mix0)$Mixture
qf_mix <- prep_mouse(qf_mix0, "proteins")
message(sprintf("fraction-run scope: %d proteins x %d columns; mixture scope: %d proteins x %d columns",
                nrow(qf_run[["proteins"]]), ncol(qf_run[["proteins"]]), nrow(qf_mix[["proteins"]]), ncol(qf_mix[["proteins"]])))
design_counts <- as.data.table(as.data.frame(colData(qf_mix)))[, .N, by = .(Diet, Duration, Mixture)][order(Diet, Duration, Mixture)]
write_tsv(design_counts, file.path(res_dir, "mouse_design_counts_used.tsv"))

## ---- Design figure ------------------------------------------------------

## VisualizeDesign shows which parameters build each cell mean, which is
## where the four contrasts in config/contrasts_mouse.tsv come from.
vd <- VisualizeDesign(sampleData = as.data.frame(colData(qf_mix)), designFormula = ~ Diet * Duration, textSizeFitted = 4)
p_design <- vd$plotlist[[1]] + ggtitle("Cell means as parameter combinations, ~ Diet * Duration, reference LF and Short")
p_cooc <- vd$cooccurrenceplots[[1]] + ggtitle("Mice per cell (all mixtures)")
save_fig(p_design / p_cooc + plot_layout(heights = c(2, 1.4)), file.path(fig_dir, "mouse_design_matrix.png"), 7.5, 7)

## ---- Models -------------------------------------------------------------

formulas <- list(
    run_mixed = ~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | BioReplicate),
    mix_mixed = ~ Diet * Duration + (1 | Mixture),
    naive     = ~ Diet * Duration
)
formula_psm <- ~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | Run:Channel) + (1 | Run:ionID) + (1 | BioReplicate)

## Refit tiers. A protein with one value per mouse cannot carry a mouse
## effect; a protein seen in one mixture cannot carry a mixture effect.
## Each tier drops the term that the failing proteins cannot support and
## keeps the fixed effects unchanged.
tiers <- list(
    run_mixed = list(full = formulas$run_mixed,
                     no_bioreplicate = ~ Diet * Duration + (1 | Mixture) + (1 | Run),
                     mixture_only = ~ Diet * Duration + (1 | Mixture),
                     no_random = ~ Diet * Duration),
    mix_mixed = list(full = formulas$mix_mixed,
                     no_random = ~ Diet * Duration),
    naive = list(full = formulas$naive),
    psm = list(full = formula_psm,
               no_run_channel = ~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | Run:ionID) + (1 | BioReplicate),
               no_run_channel_no_bioreplicate = ~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | Run:ionID))
)

models <- data.table(
    model = c("fraction_run_mixed", "fraction_run_naive", "mixture_mixed", "mixture_naive"),
    scope = c("fraction_run", "fraction_run", "mixture", "mixture"),
    structure = c("mixed", "naive", "mixed", "naive"),
    formula = c("run_mixed", "naive", "mix_mixed", "naive")
)

all_res <- list(); model_log <- list(); tier_logs <- list()
record <- function(model_name, m_scope, m_structure, m_level, fit, formula, t0) {
    res <- fit$results
    res[, `:=`(model = model_name, scope = m_scope, structure = m_structure, level = m_level)]
    all_res[[model_name]] <<- res
    tl <- copy(fit$tier_log); tl[, model := model_name]
    tier_logs[[model_name]] <<- tl
    model_log[[model_name]] <<- data.table(
        model = model_name, scope = m_scope, structure = m_structure, level = m_level,
        formula = paste(deparse(formula), collapse = ""), n_proteins = fit$n_proteins,
        n_fit_error = sum(res[contrast == labels[1], fit_type == "fitError"]),
        n_refitted = sum(res[contrast == labels[1], !fit_tier %in% c("full", "fitError")]),
        df_prior = round(fit$df_prior, 2), fit_seconds = round(fit$fit_seconds, 1),
        seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
    message(sprintf("   %s: %d proteins, %d refitted on a simpler tier, %d fit errors remain, %.0f s",
                    model_name, fit$n_proteins, model_log[[model_name]]$n_refitted, model_log[[model_name]]$n_fit_error, model_log[[model_name]]$seconds))
}
for (k in seq_len(nrow(models))) {
    m <- models[k]
    qf <- if (m$scope == "fraction_run") qf_run else qf_mix
    message("-- ", m$model, ": ", paste(deparse(formulas[[m$formula]]), collapse = ""))
    t0 <- Sys.time()
    fit <- quiet(fit_tiered(qf, "proteins", tiers[[m$formula]], hyp, params, robust = TRUE, ridge = FALSE, contrast_labels = labels))
    record(m$model, m$scope, m$structure, "protein", fit, formulas[[m$formula]], t0)
}

if (run_psm) {
    message("-- psm_mixed: ", paste(deparse(formula_psm), collapse = ""))
    t0 <- Sys.time()
    fit <- quiet(fit_tiered(qf_run, "ions_norm", tiers$psm, hyp, params, robust = TRUE, ridge = FALSE,
                            psm_level = TRUE, fcol = "Protein.Accessions", contrast_labels = labels))
    record("psm_mixed", "fraction_run", "mixed", "psm", fit, formula_psm, t0)
    mem_checkpoint(timer, "psm-level model done")
}

results <- rbindlist(all_res)
mlog <- rbindlist(model_log)
mlog[, reference := reference]
fwrite(results, file.path(res_dir, "mouse_results.tsv.gz"), sep = "\t", compress = "gzip")
write_tsv(mlog, file.path(res_dir, "mouse_models.tsv"))
tiers_out <- rbindlist(tier_logs)
tiers_out <- merge(tiers_out, results[contrast == labels[1], .N, by = .(model, tier = fit_tier)], by = c("model", "tier"), all.x = TRUE)
setnames(tiers_out, "N", "proteins_with_result_from_tier")
write_tsv(tiers_out, file.path(res_dir, "mouse_fit_tiers.tsv"))
print(tiers_out)

## ---- Counts and comparisons ---------------------------------------------

counts <- rbindlist(lapply(c(0.01, 0.05, 0.10), function(al) results[, .(
    alpha = al, n_tested = sum(!is.na(pval)), n_fit_error = sum(fit_type == "fitError"),
    n_significant = sum(adjPval < al, na.rm = TRUE),
    n_up_in_HF = sum(adjPval < al & logFC > 0, na.rm = TRUE), n_down_in_HF = sum(adjPval < al & logFC < 0, na.rm = TRUE)
), by = .(model, scope, structure, level, contrast)]))
write_tsv(counts, file.path(res_dir, "mouse_significant_counts.tsv"))
print(dcast(counts[alpha == 0.05], contrast ~ model, value.var = "n_significant"))

## Summarisation scope: the same mixed model on the same mice, once with
## nine fraction values per mouse and once with one pooled value.
sc <- dcast(results[model %in% c("fraction_run_mixed", "mixture_mixed")], protein + contrast ~ model,
            value.var = c("logFC", "se", "adjPval", "pval"))
scope_cmp <- sc[, .(
    n_both_tested = sum(!is.na(pval_fraction_run_mixed) & !is.na(pval_mixture_mixed)),
    sig_fraction_run = sum(adjPval_fraction_run_mixed < 0.05, na.rm = TRUE),
    sig_mixture = sum(adjPval_mixture_mixed < 0.05, na.rm = TRUE),
    sig_both = sum(adjPval_fraction_run_mixed < 0.05 & adjPval_mixture_mixed < 0.05, na.rm = TRUE),
    logfc_correlation = cor(logFC_fraction_run_mixed, logFC_mixture_mixed, use = "complete.obs"),
    median_se_ratio_mixture_over_fraction = median(se_mixture_mixed / se_fraction_run_mixed, na.rm = TRUE)
), by = contrast]
write_tsv(scope_cmp, file.path(res_dir, "mouse_scope_comparison.tsv"))
print(scope_cmp)

## Naive against mixed within each scope.
nv <- dcast(results[level == "protein"], scope + protein + contrast ~ structure, value.var = c("logFC", "se", "adjPval", "pval"))
naive_cmp <- nv[, .(
    n_both_tested = sum(!is.na(pval_mixed) & !is.na(pval_naive)),
    sig_mixed = sum(adjPval_mixed < 0.05, na.rm = TRUE), sig_naive = sum(adjPval_naive < 0.05, na.rm = TRUE),
    sig_both = sum(adjPval_mixed < 0.05 & adjPval_naive < 0.05, na.rm = TRUE),
    sig_naive_only = sum(adjPval_naive < 0.05 & !(adjPval_mixed < 0.05), na.rm = TRUE),
    median_se_ratio_naive_over_mixed = median(se_naive / se_mixed, na.rm = TRUE),
    median_abs_logfc_diff = median(abs(logFC_naive - logFC_mixed), na.rm = TRUE)
), by = .(scope, contrast)]
write_tsv(naive_cmp, file.path(res_dir, "mouse_naive_vs_mixed.tsv"))
print(naive_cmp)

## ---- lme4 diagnostics: singular fits and variance components --------------

diag_run <- quiet(lmer_diagnostics(qf_run, "proteins", formulas$run_mixed))
diag_run[, `:=`(scope = "fraction_run", level = "protein")]
diag_mix <- quiet(lmer_diagnostics(qf_mix, "proteins", formulas$mix_mixed))
diag_mix[, `:=`(scope = "mixture", level = "protein")]
diag <- rbind(diag_run, diag_mix, fill = TRUE)
if (run_psm) {
    ## PSM-level diagnostics on a fixed random subset: a full pass would
    ## repeat the most expensive fit of the stage.
    diag_psm <- quiet(lmer_diagnostics(qf_run, "ions_norm", formula_psm, psm_level = TRUE, max_proteins = 400))
    diag_psm[, `:=`(scope = "fraction_run", level = "psm")]
    diag <- rbind(diag, diag_psm, fill = TRUE)
}
write_tsv(diag, file.path(res_dir, "mouse_lmer_diagnostics.tsv"))
## Why the full model fails, protein by protein, from lme4's own messages.
err <- diag[status == "error"]
err[, cause := fifelse(grepl("must have > 1 sampled level", messages), "a grouping factor with a single level (one mixture or one run)",
              fifelse(grepl("number of levels of each grouping factor must be", messages), "a grouping factor with as many levels as observations (one value per mouse)",
              fifelse(grepl("rank deficient|contrasts can be applied|0 \\(non-NA\\) cases|need at least", messages), "fixed effects not estimable (empty factorial cell or too few values)", "other")))]
err_tab <- err[, .N, by = .(scope, level, cause)][order(scope, level, -N)]
write_tsv(err_tab, file.path(res_dir, "mouse_lmer_error_causes.tsv"))
print(err_tab)
vc_cols <- grep("^var_|^sigma2$", names(diag), value = TRUE)
diag_sum <- diag[, c(list(n_proteins = .N, n_error = sum(status == "error"),
                          n_convergence_warning = sum(status == "convergence_warning"),
                          n_singular = sum(singular, na.rm = TRUE), frac_singular = mean(singular, na.rm = TRUE)),
                     lapply(.SD, function(x) median(x, na.rm = TRUE)),
                     lapply(.SD, function(x) mean(x == 0, na.rm = TRUE))),
                by = .(scope, level), .SDcols = vc_cols]
setnames(diag_sum, c(names(diag_sum)[1:7], paste0("median_", vc_cols), paste0("frac_zero_", vc_cols)))
write_tsv(diag_sum, file.path(res_dir, "mouse_lmer_diagnostics_summary.tsv"))
print(diag_sum)

## Variance component figure: where the variance sits, protein by protein.
vc_long <- melt(diag[level == "protein" & status != "error", c("protein", "scope", vc_cols), with = FALSE],
                id.vars = c("protein", "scope"), variable.name = "component", value.name = "variance")[!is.na(variance)]
vc_long[, component := sub("^var_", "", component)]
vc_long[component == "sigma2", component := "residual"]
p_vc <- ggplot(vc_long, aes(x = component, y = variance + 1e-4)) +
    geom_violin(fill = "grey90", scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3) +
    scale_y_log10() + facet_wrap(~ scope, scales = "free_x") +
    labs(x = NULL, y = "variance component (log10 scale, + 1e-4)",
         title = "Variance components per protein from lme4; a value at the floor is a component estimated as zero") +
    theme_bw(base_size = 10)
save_fig(p_vc, file.path(fig_dir, "mouse_variance_components.png"), 8, 4.2)

## ---- Volcano figures ----------------------------------------------------

desc <- as.data.table(as.data.frame(rowData(qf_run[["proteins"]])[, "Master.Protein.Descriptions", drop = FALSE]), keep.rownames = "protein")
desc[, gene := sub(".*GN=([^ ]+).*", "\\1", Master.Protein.Descriptions)]
desc[!grepl("GN=", Master.Protein.Descriptions), gene := protein]
res_v <- merge(results[model == "fraction_run_mixed"], desc, by = "protein", all.x = TRUE)
res_v[, contrast := factor(contrast, levels = labels)]
p_volc <- ggplot(res_v[!is.na(pval)], aes(x = logFC, y = -log10(pval))) +
    geom_point(aes(colour = adjPval < 0.05), size = 0.7, alpha = 0.7) +
    geom_text_repel(data = res_v[adjPval < 0.05], aes(label = gene), size = 2.4, max.overlaps = 25) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), name = "BH adjusted p < 0.05") +
    facet_wrap(~ contrast, nrow = 1) +
    labs(x = "log2 fold change (positive: higher on high fat)", y = "-log10 p",
         title = "Mouse adipose proteome: ~ Diet * Duration + (1 | Mixture) + (1 | Run) + (1 | BioReplicate), fraction-run scope") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
save_fig(p_volc, file.path(fig_dir, "mouse_volcano.png"), 12, 4.2)

## Scope and structure side by side for the 8 week diet contrast.
side <- merge(results[contrast == "HF_vs_LF_8wk" & level == "protein", .(protein, model, logFC, adjPval)],
              results[contrast == "HF_vs_LF_8wk" & model == "fraction_run_mixed", .(protein, logFC_ref = logFC)], by = "protein")
side[, model := factor(model, levels = models$model)]
p_side <- ggplot(side[!is.na(adjPval)], aes(x = logFC_ref, y = logFC, colour = adjPval < 0.05)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey70") +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), name = "significant in this model") +
    facet_wrap(~ model, nrow = 1) +
    labs(x = "log2 fold change, fraction-run mixed model", y = "log2 fold change, this model",
         title = "HF against LF at 8 weeks: the same proteins under four models") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
save_fig(p_side, file.path(fig_dir, "mouse_scope_comparison.png"), 11, 3.6)

## ---- Branched chain amino acid metabolism -------------------------------
##
## Plubell et al. reported that proteins of branched chain amino acid
## catabolism (Mccc1, Pccb, Hibadh among their key drivers) decreased
## under short term high fat feeding. Two checks, both on the fraction-
## run mixed model's 8 week contrast, both with the tested proteins as
## the universe.
##
## 1. A keyword check on the Master Protein Descriptions. This is not an
##    enrichment analysis: the term list is hand written, the matching
##    is a regular expression on free text, and the p-value from
##    Fisher's test is descriptive. It is reported because it needs no
##    external resource.
## 2. KEGG pathway mmu00280 (valine, leucine and isoleucine degradation),
##    fetched from the KEGG REST API and mapped to UniProt accessions,
##    with the background being exactly the proteins tested. That is a
##    properly constructed enrichment test for one pathway, and it is
##    cached in config/ with the date fetched.
ctr8 <- merge(results[model == "fraction_run_mixed" & contrast == "HF_vs_LF_8wk" & !is.na(pval)], desc, by = "protein", all.x = TRUE)
bcaa_terms <- c("branched-chain", "branched chain", "2-oxoisovalerate", "isovaleryl-CoA", "methylcrotonoyl",
                "3-hydroxyisobutyr", "methylmalonate", "methylmalonyl", "propionyl-CoA carboxylase",
                "3-methylglutaconyl", "hydroxymethylglutaryl-CoA lyase", "short/branched", "dihydrolipoyl transacylase",
                "enoyl-CoA hydratase", "3-hydroxyacyl-CoA dehydrogenase", "isobutyryl")
pattern <- paste(bcaa_terms, collapse = "|")
ctr8[, bcaa_keyword := grepl(pattern, Master.Protein.Descriptions, ignore.case = TRUE)]
ctr8[, significant := adjPval < 0.05]
ctr8[, significant_10 := adjPval < 0.10]
kw_tab <- ctr8[, .(n = .N, n_sig_5 = sum(significant), n_sig_10 = sum(significant_10),
                   n_down_sig_10 = sum(significant_10 & logFC < 0), median_logfc = median(logFC)), by = bcaa_keyword]
ft <- fisher.test(table(ctr8$bcaa_keyword, ctr8$significant_10))
kw_hits <- ctr8[bcaa_keyword == TRUE, .(protein, gene, description = Master.Protein.Descriptions, logFC, pval, adjPval)][order(adjPval)]
write_tsv(kw_hits, file.path(res_dir, "mouse_bcaa_keyword_hits.tsv"))
write_tsv(data.table(check = "keyword", contrast = "HF_vs_LF_8wk", model = "fraction_run_mixed", terms = pattern,
                     n_tested = nrow(ctr8), n_keyword = sum(ctr8$bcaa_keyword), n_keyword_sig_10 = kw_tab[bcaa_keyword == TRUE, n_sig_10],
                     n_keyword_down_sig_10 = kw_tab[bcaa_keyword == TRUE, n_down_sig_10],
                     n_keyword_negative_logfc = sum(ctr8$bcaa_keyword & ctr8$logFC < 0),
                     median_logfc_keyword = kw_tab[bcaa_keyword == TRUE, median_logfc],
                     median_logfc_others = kw_tab[bcaa_keyword == FALSE, median_logfc],
                     fisher_odds_ratio_sig_10 = unname(ft$estimate), fisher_p_sig_10 = ft$p.value,
                     wilcoxon_p_logfc = wilcox.test(logFC ~ bcaa_keyword, data = ctr8)$p.value),
          file.path(res_dir, "mouse_bcaa_keyword_check.tsv"))
message(sprintf("keyword check: %d of %d tested proteins match; %d significant at 10%% FDR, %d of them down on HF; median logFC %.3f vs %.3f",
                sum(ctr8$bcaa_keyword), nrow(ctr8), kw_tab[bcaa_keyword == TRUE, n_sig_10], kw_tab[bcaa_keyword == TRUE, n_down_sig_10],
                kw_tab[bcaa_keyword == TRUE, median_logfc], kw_tab[bcaa_keyword == FALSE, median_logfc]))

kegg_file <- file.path(root, "config", "kegg_mmu00280_uniprot.tsv")
kegg <- NULL
if (file.exists(kegg_file)) {
    kegg <- fread(kegg_file)
} else {
    kegg <- tryCatch({
        genes <- readLines("https://rest.kegg.jp/link/mmu/path:mmu00280", warn = FALSE)
        genes <- sub("^.*\t", "", genes)
        up <- character()
        for (chunk in split(genes, ceiling(seq_along(genes) / 50))) {
            conv <- readLines(paste0("https://rest.kegg.jp/conv/uniprot/", paste(chunk, collapse = "+")), warn = FALSE)
            up <- c(up, conv)
        }
        k <- data.table(kegg_gene = sub("\t.*$", "", up), uniprot = sub("^.*up:", "", up))
        k[, `:=`(pathway = "mmu00280", pathway_name = "Valine, leucine and isoleucine degradation",
                 fetched = format(Sys.Date()), source = "https://rest.kegg.jp")]
        fwrite(k, kegg_file, sep = "\t")
        k
    }, error = function(e) { message("KEGG fetch failed: ", conditionMessage(e)); NULL })
}
if (!is.null(kegg)) {
    ctr8[, in_kegg := protein %in% kegg$uniprot]
    tab <- table(kegg = ctr8$in_kegg, sig = ctr8$significant_10)
    ft_k <- fisher.test(tab)
    tab5 <- table(kegg = ctr8$in_kegg, sig = ctr8$significant)
    ft_k5 <- fisher.test(tab5)
    kegg_out <- data.table(
        check = "KEGG mmu00280 enrichment, background = tested proteins", contrast = "HF_vs_LF_8wk", model = "fraction_run_mixed",
        n_pathway_uniprot = length(unique(kegg$uniprot)), n_tested = nrow(ctr8), n_pathway_tested = sum(ctr8$in_kegg),
        n_pathway_sig_10 = sum(ctr8$in_kegg & ctr8$significant_10), n_other_sig_10 = sum(!ctr8$in_kegg & ctr8$significant_10),
        odds_ratio_10 = unname(ft_k$estimate), fisher_p_10 = ft_k$p.value,
        n_pathway_sig_5 = sum(ctr8$in_kegg & ctr8$significant), odds_ratio_5 = unname(ft_k5$estimate), fisher_p_5 = ft_k5$p.value,
        n_pathway_negative_logfc = sum(ctr8$in_kegg & ctr8$logFC < 0),
        median_logfc_pathway = median(ctr8[in_kegg == TRUE, logFC]), median_logfc_others = median(ctr8[in_kegg == FALSE, logFC]),
        wilcoxon_p_logfc = wilcox.test(logFC ~ in_kegg, data = ctr8)$p.value,
        kegg_fetched = kegg$fetched[1])
    write_tsv(kegg_out, file.path(res_dir, "mouse_bcaa_kegg_enrichment.tsv"))
    write_tsv(ctr8[in_kegg == TRUE, .(protein, gene, description = Master.Protein.Descriptions, logFC, pval, adjPval)][order(adjPval)],
              file.path(res_dir, "mouse_bcaa_kegg_members.tsv"))
    message(sprintf("KEGG mmu00280: %d of %d tested proteins in pathway; %d significant at 10%% FDR (others: %d); Fisher p = %.3g; %d of %d with negative logFC",
                    sum(ctr8$in_kegg), nrow(ctr8), kegg_out$n_pathway_sig_10, kegg_out$n_other_sig_10, kegg_out$fisher_p_10,
                    kegg_out$n_pathway_negative_logfc, sum(ctr8$in_kegg)))
    ## Same check under the other models, so that the conclusion's
    ## dependence on the model is visible.
    per_model <- rbindlist(lapply(unique(results$model), function(mm) {
        d <- results[model == mm & contrast == "HF_vs_LF_8wk" & !is.na(pval)]
        d[, in_kegg := protein %in% kegg$uniprot]
        data.table(model = mm, n_tested = nrow(d), n_pathway = sum(d$in_kegg),
                   n_pathway_sig_10 = sum(d$in_kegg & d$adjPval < 0.10), n_pathway_sig_5 = sum(d$in_kegg & d$adjPval < 0.05),
                   n_pathway_negative = sum(d$in_kegg & d$logFC < 0), median_logfc_pathway = median(d[in_kegg == TRUE, logFC]),
                   fisher_p_10 = fisher.test(table(d$in_kegg, d$adjPval < 0.10))$p.value)
    }))
    write_tsv(per_model, file.path(res_dir, "mouse_bcaa_kegg_by_model.tsv"))
    print(per_model)
}

stage_end(timer, note = sprintf("reference=%s; psm_level=%s", reference, run_psm))
stage_mark_done(STAGE, conf)
