## Stage 05: the spike-in benchmark against known truth.
##
## Every UPS1 protein changes between any two dilutions by a ratio the
## design fixed; no HeLa protein does. So for each of the six pairwise
## contrasts every workflow variant can be scored: true positives,
## false positives, realised false discovery proportion at nominal 5 per
## cent, sensitivity, and the bias of the fold change estimate.
##
## Variants are rows of a table crossing
##   level      protein-level model on robustSummary values, or the
##              PSM-level mixed model on all ions of a protein
##   reference  drop, keep or ratio (scripts/lib/reference_channels.R)
##   structure  mixed: (1 | Run) + (1 | Mixture); run_only: (1 | Run);
##              naive: no random effects, every channel measurement
##              treated as an independent sample
##   robust     M-estimation on or off
##   ridge      off by default; one timed run with it on
## plus limma with duplicateCorrelation blocked on mixture as the
## reference arm.
##
## Two designs. "balanced" is the experiment as acquired: every condition
## in every mixture, so the mixture effect cancels inside every contrast.
## "unbalanced" is a defined channel subset in which conditions 1 and
## 0.667 are present only in mixtures 1 to 3 and conditions 0.5 and
## 0.125 only in mixtures 3 to 5, which partially confounds the large
## contrasts with mixture, as the three-plex mouse design does. The
## naive model's behaviour differs between the two, and both are scored
## against the same truth.
##
## The PSM-level fits are the expensive part. The script times a pilot
## batch, projects the full cost, and fits the largest random subset of
## HeLa proteins (always all UPS1 proteins) that fits the time budget.
## What ran at full scale and what did not is written to
## results/spikein1_psm_subset.tsv.
##
## Usage: scripts/lib/rscript.sh scripts/05_benchmark_workflows.R
##          [--acquisition both|ms3|ms2] [--psm-proteins auto|all|N]
##          [--psm-budget-min 90] [--workers N] [--force]
##
## Model formulas, stated here, in the fit calls below and in the README:
##   protein, mixed : ~ 0 + Condition + (1 | Run) + (1 | Mixture)
##   protein, naive : ~ 0 + Condition
##   PSM, mixed     : ~ 0 + Condition + (1 | Run) + (1 | Mixture) + (1 | Run:Channel) + (1 | Run:ionID)

suppressPackageStartupMessages({
    library("data.table")
    library("QFeatures")
    library("msqrob2")
    library("BiocParallel")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
source(file.path(root, "scripts", "lib", "common.R"))
source(file.path(root, "scripts", "lib", "reference_channels.R"))
source(file.path(root, "scripts", "lib", "msqrob_helpers.R"))
conf <- read_conf(root)
opt <- parse_args()
STAGE <- "05_benchmark_workflows"
## --postprocess-only recomputes every summary table from the saved
## per-protein results without refitting anything.
postprocess_only <- isTRUE(opt$`postprocess-only`)
if (!postprocess_only && stage_already_done(STAGE, conf, force = opt$force)) quit(save = "no")
timer <- stage_begin(STAGE, conf)
set.seed(20260921)
register_parallel(if (is.null(opt$workers)) conf$WORKERS else opt$workers)

acq <- if (is.null(opt$acquisition)) "both" else opt$acquisition
acquisitions <- switch(acq, both = c("ms3", "ms2"), ms3 = "ms3", ms2 = "ms2", stop("unknown acquisition"))
psm_proteins <- if (is.null(opt$`psm-proteins`)) "auto" else opt$`psm-proteins`
psm_budget_s <- 60 * as.numeric(if (is.null(opt$`psm-budget-min`)) 120 else opt$`psm-budget-min`)

res_dir <- file.path(root, "results")
qf_dir <- file.path(root, "results", "qfeatures")

## ---- Design and contrasts --------------------------------------------------

contr <- fread(file.path(root, "config", "contrasts_spikein1.tsv"))
expected <- spikein_expected_contrasts()
stopifnot(all(contr$contrast %in% expected$contrast))
cond_params <- paste0("Condition", spikein_conditions()$condition)
hyp <- contr$hypothesis
labels <- contr$contrast

## The three random effect structures. Mixture is a random effect
## because all ten channels of a plex share one labelling reaction, one
## pooling step and every scan of every run; treating its 30 channel
## measurements per condition as 30 samples is the pseudo-replication
## this stage measures. Run within mixture carries the technical
## replicate injections. The naive model has neither.
formulas <- list(
    mixed    = ~ 0 + Condition + (1 | Run) + (1 | Mixture),
    run_only = ~ 0 + Condition + (1 | Run),
    naive    = ~ 0 + Condition
)
formula_psm <- ~ 0 + Condition + (1 | Run) + (1 | Mixture) + (1 | Run:Channel) + (1 | Run:ionID)

## Variant table. ms2 runs the subset needed by the compression stage.
variants <- rbindlist(list(
    data.table(level = "protein", reference = rep(c("drop", "keep", "ratio"), each = 4),
               structure = rep(c("mixed", "naive"), times = 6),
               robust = rep(c(TRUE, TRUE, FALSE, FALSE), times = 3), ridge = FALSE, engine = "msqrob2"),
    data.table(level = "protein", reference = "drop", structure = "run_only", robust = TRUE, ridge = FALSE, engine = "msqrob2"),
    data.table(level = "protein", reference = "drop", structure = "mixed", robust = TRUE, ridge = TRUE, engine = "msqrob2"),
    data.table(level = "protein", reference = c("drop", "keep", "ratio"), structure = "dupcor_mixture", robust = FALSE, ridge = FALSE, engine = "limma"),
    data.table(level = "psm", reference = c("drop", "keep", "drop", "drop"), structure = "mixed",
               robust = c(TRUE, TRUE, FALSE, TRUE), ridge = c(FALSE, FALSE, FALSE, TRUE), engine = "msqrob2")
))
variants[, variant := paste(level, reference, structure, ifelse(robust, "robust", "ls"), ifelse(ridge, "ridge", "noridge"), engine, sep = "_")]
variants[, core := (level == "protein" & reference == "drop" & structure %in% c("mixed", "naive") & robust & !ridge & engine == "msqrob2") |
                   (engine == "limma" & reference == "drop")]
message(nrow(variants), " variants defined")

## The unbalanced design: which mixtures carry which conditions.
unbalanced_rule <- function(qf) {
    cd <- colData(qf)
    keep <- (cd$Condition %in% c("1", "0.667") & cd$Mixture %in% c("Mixture1", "Mixture2", "Mixture3")) |
            (cd$Condition %in% c("0.5", "0.125") & cd$Mixture %in% c("Mixture3", "Mixture4", "Mixture5")) |
            cd$Condition == "Norm"
    subsetByColData(qf, keep)
}
designs <- list(balanced = identity, unbalanced = unbalanced_rule)

## ---- Helpers -------------------------------------------------------------

prepare <- function(qf, reference) {
    r <- apply_reference_treatment(qf, mode = reference, sets = c("ions_norm", "proteins"), run_col = "Run")
    r$qf
}

## PSM-level subset. All UPS1 proteins always; HeLa proteins sampled
## with a fixed seed to whatever count the budget allows.
choose_psm_subset <- function(qf, budget_s) {
    ions <- qf[["ions_norm"]]
    prots <- unique(rowData(ions)$Protein.Accessions)
    ups <- grep("ups", prots, value = TRUE)
    hela <- setdiff(prots, ups)
    if (psm_proteins == "all") return(list(proteins = prots, note = "all proteins requested", per_protein_s = NA, projected_full_s = NA))
    if (psm_proteins != "auto") {
        n <- as.integer(psm_proteins)
        set.seed(20260921)
        return(list(proteins = c(ups, sample(hela, min(n, length(hela)))), note = sprintf("%d HeLa proteins requested", n), per_protein_s = NA, projected_full_s = NA))
    }
    set.seed(20260921)
    pilot <- sample(hela, 40)
    sub <- qf[rowData(ions)$Protein.Accessions %in% pilot, , "ions_norm"]
    t0 <- Sys.time()
    invisible(msqrobAggregate(sub, i = "ions_norm", fcol = "Protein.Accessions", formula = formula_psm,
                              robust = TRUE, ridge = FALSE, name = "pilot", aggregateFun = MsCoreUtils::robustSummary))
    per <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / length(pilot)
    ## Four PSM-level variants run on the subset: the least squares one
    ## costs about half a robust fit, the ridge one about one and a half,
    ## and the reference-kept one has 150 columns instead of 120; the
    ## factor 4.5 is that sum.
    projected_full <- per * length(prots) * 4.5
    n_afford <- floor(budget_s / (per * 4.5))
    if (n_afford >= length(prots)) {
        list(proteins = prots, note = sprintf("pilot %.2f s per protein; all %d proteins fit the %.0f min budget", per, length(prots), budget_s / 60),
             per_protein_s = per, projected_full_s = projected_full)
    } else {
        n_hela <- max(0, n_afford - length(ups))
        set.seed(20260921)
        list(proteins = c(ups, sample(hela, n_hela)),
             note = sprintf("pilot %.2f s per protein; full scale projected at %.0f min for four variants, over the %.0f min budget; %d of %d HeLa proteins sampled",
                            per, projected_full / 60, budget_s / 60, n_hela, length(hela)),
             per_protein_s = per, projected_full_s = projected_full)
    }
}

## ---- Main loop -----------------------------------------------------------

all_results <- list(); all_metrics <- list(); variant_log <- list(); diag_out <- list(); subset_log <- list()
if (postprocess_only) { acquisitions <- character(); designs <- list() }

## lme4 prints "boundary (singular) fit" for every singular protein and
## msqrob2 does not silence it; thousands of those lines would bury the
## log. Singular fits are counted properly by the diagnostics pass.
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

for (a in acquisitions) for (dsg in names(designs)) {
    if (a == "ms2" && dsg == "unbalanced") next
    ds <- if (a == "ms3") "spikein1" else "spikein1_ms2"
    message("==== acquisition ", a, " (", ds, "), design ", dsg, " ====")
    qf0 <- designs[[dsg]](readRDS(file.path(qf_dir, paste0(ds, "_preprocessed.rds"))))
    vset <- if (a == "ms3" && dsg == "balanced") variants else variants[core == TRUE]
    psm_sub <- NULL

    for (v in seq_len(nrow(vset))) {
        vr <- vset[v]
        message(sprintf("-- [%s/%s] variant %d/%d: %s", a, dsg, v, nrow(vset), vr$variant))
        qf <- prepare(qf0, vr$reference)
        cd <- as.data.frame(colData(qf))
        params <- cond_params
        t_start <- Sys.time()

        if (vr$engine == "limma") {
            mat <- assay(qf[["proteins"]])
            fit <- fit_limma_dupcor(mat, cd, ~ 0 + Condition, block = "Mixture", hyp, params, contrast_labels = labels)
            extra <- sprintf("consensus correlation %.3f", fit$consensus_correlation)
        } else if (vr$level == "protein") {
            fit <- quiet(fit_protein_level(qf, "proteins", formulas[[vr$structure]], hyp, params,
                                           robust = vr$robust, ridge = vr$ridge, contrast_labels = labels))
            extra <- ""
        } else {
            if (is.null(psm_sub)) {
                psm_sub <- choose_psm_subset(prepare(qf0, "drop"), psm_budget_s)
                message("PSM subset: ", psm_sub$note)
                subset_log[[a]] <- data.table(acquisition = a, n_proteins_total = length(unique(rowData(qf0[["ions_norm"]])$Protein.Accessions)),
                                              n_proteins_fitted = length(psm_sub$proteins),
                                              n_ups = sum(grepl("ups", psm_sub$proteins)),
                                              pilot_seconds_per_protein = psm_sub$per_protein_s,
                                              projected_full_scale_min = psm_sub$projected_full_s / 60,
                                              budget_min = psm_budget_s / 60, note = psm_sub$note)
            }
            keep <- rowData(qf[["ions_norm"]])$Protein.Accessions %in% psm_sub$proteins
            qf_s <- qf[keep, , "ions_norm"]
            fit <- quiet(fit_psm_level(qf_s, "ions_norm", formula_psm, "Protein.Accessions", hyp, params,
                                       robust = vr$robust, ridge = vr$ridge, contrast_labels = labels))
            extra <- sprintf("%d proteins in subset", length(psm_sub$proteins))
            rm(qf_s)
        }
        elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
        res <- fit$results
        res[, `:=`(acquisition = a, design = dsg, variant = vr$variant, level = vr$level, reference = vr$reference,
                   structure = vr$structure, robust = vr$robust, ridge = vr$ridge, engine = vr$engine)]
        all_results[[length(all_results) + 1L]] <- res
        met <- rbindlist(lapply(c(0.01, 0.05, 0.10), function(al) benchmark_metrics(res, expected, alpha = al)))
        met[, `:=`(acquisition = a, design = dsg, variant = vr$variant, level = vr$level, reference = vr$reference,
                   structure = vr$structure, robust = vr$robust, ridge = vr$ridge, engine = vr$engine)]
        all_metrics[[length(all_metrics) + 1L]] <- met
        variant_log[[length(variant_log) + 1L]] <- data.table(
            acquisition = a, design = dsg, variant = vr$variant, level = vr$level, reference = vr$reference, structure = vr$structure,
            robust = vr$robust, ridge = vr$ridge, engine = vr$engine,
            n_proteins = fit$n_proteins, n_fit_error = sum(res[contrast == labels[1], fit_type == "fitError"]),
            fit_seconds = round(fit$fit_seconds, 1), variant_seconds = round(elapsed, 1),
            n_channels = ncol(qf[["proteins"]]), note = extra, peak_rss_mb = round(peak_rss_mb()))
        m5 <- met[alpha == 0.05]
        message(sprintf("   %.0f s; fit errors %d; FDP at 5%%: %s", elapsed, variant_log[[length(variant_log)]]$n_fit_error,
                        paste(sprintf("%s=%.3f", m5$contrast, m5$fdp), collapse = " ")))
        rm(qf, fit, res); invisible(gc())
    }

    ## Singular fit and variance component diagnostics for the adopted
    ## structures, from a plain lme4 pass (see msqrob_helpers.R).
    message("-- [", a, "/", dsg, "] lme4 diagnostics, protein level, reference dropped")
    qf <- prepare(qf0, "drop")
    d <- quiet(lmer_diagnostics(qf, "proteins", formulas$mixed))
    d[, `:=`(acquisition = a, design = dsg, level = "protein", reference = "drop", structure = "mixed")]
    diag_out[[length(diag_out) + 1L]] <- d
    if (a == "ms3" && dsg == "balanced") {
        qf <- prepare(qf0, "ratio")
        d <- quiet(lmer_diagnostics(qf, "proteins", formulas$mixed))
        d[, `:=`(acquisition = a, design = dsg, level = "protein", reference = "ratio", structure = "mixed")]
        diag_out[[length(diag_out) + 1L]] <- d
        if (!is.null(psm_sub)) {
            message("-- [", a, "] lme4 diagnostics, PSM level subset, reference dropped")
            qf <- prepare(qf0, "drop")
            keep <- rowData(qf[["ions_norm"]])$Protein.Accessions %in% psm_sub$proteins
            d <- quiet(lmer_diagnostics(qf[keep, , "ions_norm"], "ions_norm", formula_psm, psm_level = TRUE))
            d[, `:=`(acquisition = a, design = dsg, level = "psm", reference = "drop", structure = "mixed")]
            diag_out[[length(diag_out) + 1L]] <- d
        }
    }
    rm(qf, qf0); invisible(gc())
    mem_checkpoint(timer, paste("acquisition", a, "design", dsg, "done"))
}

## ---- Write --------------------------------------------------------------

if (postprocess_only) {
    results <- fread(file.path(res_dir, "spikein1_benchmark_results.tsv.gz"))
    vlog <- fread(file.path(res_dir, "spikein1_benchmark_variants.tsv"))
    diag <- fread(file.path(res_dir, "spikein1_lmer_diagnostics.tsv"))
    metrics <- rbindlist(lapply(split(results, by = c("acquisition", "design", "variant")), function(r) {
        met <- rbindlist(lapply(c(0.01, 0.05, 0.10), function(al) benchmark_metrics(r, expected, alpha = al)))
        met[, `:=`(acquisition = r$acquisition[1], design = r$design[1], variant = r$variant[1], level = r$level[1],
                   reference = r$reference[1], structure = r$structure[1], robust = r$robust[1], ridge = r$ridge[1], engine = r$engine[1])]
        met
    }))
    message("postprocess-only: reloaded ", nrow(results), " result rows")
} else {
    results <- rbindlist(all_results)
    metrics <- rbindlist(all_metrics)
    vlog <- rbindlist(variant_log)
    diag <- rbindlist(diag_out, fill = TRUE)
    fwrite(results, file.path(res_dir, "spikein1_benchmark_results.tsv.gz"), sep = "\t", compress = "gzip")
    message("wrote results/spikein1_benchmark_results.tsv.gz (", nrow(results), " rows)")
    write_tsv(vlog, file.path(res_dir, "spikein1_benchmark_variants.tsv"))
    write_tsv(diag, file.path(res_dir, "spikein1_lmer_diagnostics.tsv"))
    if (length(subset_log)) write_tsv(rbindlist(subset_log), file.path(res_dir, "spikein1_psm_subset.tsv"))
}
write_tsv(metrics, file.path(res_dir, "spikein1_benchmark_metrics.tsv"))

## Diagnostics summary: singular fits and the size of each variance component.
vc_cols <- grep("^var_|^sigma2$", names(diag), value = TRUE)
diag_sum <- diag[, c(list(n_proteins = .N, n_error = sum(status == "error"),
                          n_convergence_warning = sum(status == "convergence_warning"),
                          n_singular = sum(singular, na.rm = TRUE),
                          frac_singular = mean(singular, na.rm = TRUE)),
                     lapply(.SD, function(x) median(x, na.rm = TRUE))),
                by = .(acquisition, design, level, reference, structure), .SDcols = vc_cols]
setnames(diag_sum, vc_cols, paste0("median_", vc_cols))
write_tsv(diag_sum, file.path(res_dir, "spikein1_lmer_diagnostics_summary.tsv"))

## Naive against mixed, head to head, on the same proteins.
nm <- results[acquisition == "ms3" & level == "protein" & reference == "drop" & robust & !ridge & engine == "msqrob2" & structure %in% c("mixed", "naive")]
nm <- dcast(nm, design + protein + contrast ~ structure, value.var = c("logFC", "se", "adjPval", "pval"))
nm[, is_ups := grepl("ups", protein)]
head_to_head <- nm[, .(
    n_proteins = sum(!is.na(pval_mixed) & !is.na(pval_naive)),
    median_se_ratio_naive_over_mixed = median(se_naive / se_mixed, na.rm = TRUE),
    median_logfc_diff_naive_minus_mixed = median(logFC_naive - logFC_mixed, na.rm = TRUE),
    sig_mixed = sum(adjPval_mixed < 0.05, na.rm = TRUE), sig_naive = sum(adjPval_naive < 0.05, na.rm = TRUE),
    sig_both = sum(adjPval_mixed < 0.05 & adjPval_naive < 0.05, na.rm = TRUE),
    fp_mixed = sum(adjPval_mixed < 0.05 & !is_ups, na.rm = TRUE), fp_naive = sum(adjPval_naive < 0.05 & !is_ups, na.rm = TRUE),
    tp_mixed = sum(adjPval_mixed < 0.05 & is_ups, na.rm = TRUE), tp_naive = sum(adjPval_naive < 0.05 & is_ups, na.rm = TRUE)
), by = .(design, contrast)]
head_to_head[, fdp_mixed := fp_mixed / sig_mixed]
head_to_head[, fdp_naive := fp_naive / sig_naive]
write_tsv(head_to_head, file.path(res_dir, "spikein1_naive_vs_mixed.tsv"))
print(head_to_head[, .(design, contrast, sig_mixed, fp_mixed, fdp_mixed, sig_naive, fp_naive, fdp_naive, median_se_ratio_naive_over_mixed)])

## What the false positives are. A HeLa protein called significant by
## the mixed model is either a genuinely wrong call or a small,
## run-consistent shift that the design produces (the channel with more
## UPS1 carries less HeLa per unit of total peptide). Their sign and
## size say which. The same table is written for every variant so that
## the pattern can be compared across models.
fp <- results[!is.na(adjPval) & adjPval < 0.05 & !grepl("ups", protein)]
fp_summary <- fp[, .(n_fp = .N, median_logfc = median(logFC), median_abs_logfc = median(abs(logFC)),
                     frac_negative = mean(logFC < 0), frac_abs_below_0.1 = mean(abs(logFC) < 0.1),
                     frac_abs_below_0.2 = mean(abs(logFC) < 0.2), median_se = median(se)),
                 by = .(acquisition, design, variant, contrast)]
write_tsv(fp_summary, file.path(res_dir, "spikein1_false_positive_summary.tsv"))

## Sensitivity of the FDP to an effect size floor, post hoc and reported
## for every floor tried. The floors are arbitrary; none was used to
## pick a result. A floor of 0 is the plain BH list.
floors <- c(0, 0.1, 0.2, 0.3, 0.5)
fdp_floor <- rbindlist(lapply(floors, function(fl) {
    r <- results[!is.na(adjPval)]
    r[, called := adjPval < 0.05 & abs(logFC) >= fl]
    r[, is_ups := grepl("ups", protein)]
    r[, .(logfc_floor = fl, n_called = sum(called), tp = sum(called & is_ups), fp = sum(called & !is_ups),
          fdp = if (sum(called) == 0) NA_real_ else sum(called & !is_ups) / sum(called),
          sensitivity = sum(called & is_ups) / sum(is_ups)),
      by = .(acquisition, design, variant, contrast)]
}))
write_tsv(fdp_floor, file.path(res_dir, "spikein1_fdp_by_effect_size_floor.tsv"))

## Which channels the unbalanced design kept, for the README.
ub <- unbalanced_rule(readRDS(file.path(qf_dir, "spikein1_preprocessed.rds")))
write_tsv(as.data.table(as.data.frame(colData(ub)))[Condition != "Norm", .N, by = .(Mixture, Condition)][order(Mixture, Condition)],
          file.path(res_dir, "spikein1_unbalanced_design.tsv"))

if (postprocess_only) {
    message("postprocess-only: summary tables rewritten")
} else {
    stage_end(timer, note = sprintf("acquisitions=%s; %d variants", paste(acquisitions, collapse = ","), nrow(vlog)))
    stage_mark_done(STAGE, conf)
}
