## Stage 04: from raw reporter intensities to model-ready assays.
##
## For each dataset, in this order: duplicate-spectrum filter (spike-in
## only), missing value encoding, log2, PSM quality filters, per-run
## median centering, robust summarisation to protein within run (and
## within mixture for the mouse data), QC figures, MDS.
##
## The reference channels are kept in the saved objects. Choosing what to
## do with them is the benchmark's job (05), through the functions in
## scripts/lib/reference_channels.R.
##
## Usage: scripts/lib/rscript.sh scripts/04_preprocess.R [--dataset spikein1|spikein1_ms2|mouse|all] [--force]
##
## Outputs:
##   results/qfeatures/<dataset>_preprocessed.rds         run-level assays: ions_norm, proteins
##   results/qfeatures/mouse_mixture_preprocessed.rds     mixture-level assay: proteins (mouse only)
##   results/<dataset>_filter_log.tsv                     every filter with counts and reason
##   results/<dataset>_missingness.tsv                    channel- and run-level missingness
##   results/<dataset>_normalisation_factors.tsv          per-channel medians before centering
##   results/<dataset>_mds.tsv                            MDS coordinates
##   figures/<dataset>_normalisation_densities.png
##   figures/<dataset>_summarisation_example.png
##   figures/<dataset>_mds.png

suppressPackageStartupMessages({
    library("data.table")
    library("QFeatures")
    library("MsCoreUtils")
    library("ggplot2")
    library("patchwork")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
source(file.path(root, "scripts", "lib", "common.R"))
conf <- read_conf(root)
opt <- parse_args()
STAGE <- "04_preprocess"
if (stage_already_done(STAGE, conf, force = opt$force)) quit(save = "no")
timer <- stage_begin(STAGE, conf)
set.seed(20260921)

which_ds <- if (is.null(opt$dataset) || opt$dataset == "all") c("spikein1", "spikein1_ms2", "mouse") else opt$dataset
qf_dir <- file.path(root, "results", "qfeatures")
fig_dir <- file.path(root, "figures")
res_dir <- file.path(root, "results")

## ---- Filter bookkeeping ----------------------------------------------------

## Every filter appends one row: what it removed, the threshold, and why
## the criterion is independent of the test statistic. That independence
## is what keeps filtering from biasing the p-value distribution: a
## filter that looked at fold changes or p-values would select for
## proteins that happen to look differential, and the null distribution
## of the survivors would no longer be uniform. Every criterion here is
## a property of identification or of measurement completeness and can
## be evaluated without knowing which channel is which condition.
n_psms <- function(qf, sets) sum(vapply(sets, function(i) nrow(qf[[i]]), integer(1)))
filter_log <- list()
log_filter <- function(step, before, after, criterion, why) {
    filter_log[[length(filter_log) + 1L]] <<- data.table(
        step = step, psms_before = before, psms_after = after, removed = before - after,
        criterion = criterion, why = why)
    message(sprintf("  %-38s %8d -> %8d (removed %d)", step, before, after, before - after))
    mem_checkpoint(timer, step)
}

## Recompute a rowData column across all run assays from a data.table
## expression on the stacked rowData. Cheaper than rbindRowData for
## the columns we need and it keeps the assay index.
stacked_rowdata <- function(qf, sets, cols) {
    rbindlist(lapply(sets, function(i) {
        rd <- rowData(qf[[i]])[, cols, drop = FALSE]
        dt <- as.data.table(as.data.frame(rd))
        dt[, assay := i]
        dt[, rowname := rownames(rd)]
        dt
    }))
}

set_rowdata_column <- function(qf, sets, dt, column) {
    for (i in sets) {
        d <- dt[assay == i]
        rd <- rowData(qf[[i]])
        rd[[column]] <- d[[column]][match(rownames(rd), d$rowname)]
        rowData(qf[[i]]) <- rd
    }
    qf
}

## ---- MDS ----------------------------------------------------------------

mds_figure <- function(m, cd, dataset, facet = FALSE, label = dataset) {
    keep <- cd[Condition != "Norm" & Condition != "Long_M", sample]
    m <- m[, intersect(colnames(m), keep), drop = FALSE]
    complete <- m[rowSums(is.na(m)) == 0, , drop = FALSE]
    d <- dist(t(complete))
    mds <- cmdscale(d, k = 2, eig = TRUE)
    coords <- data.table(sample = rownames(mds$points), dim1 = mds$points[, 1], dim2 = mds$points[, 2])
    coords <- merge(coords, cd, by = "sample")
    coords[, dataset := dataset]
    coords[, n_complete_proteins := nrow(complete)]
    coords[, var_explained_dim1 := mds$eig[1] / sum(pmax(mds$eig, 0))]
    coords[, var_explained_dim2 := mds$eig[2] / sum(pmax(mds$eig, 0))]
    write_tsv(coords, file.path(res_dir, paste0(dataset, "_mds.tsv")))
    p <- ggplot(coords, aes(x = dim1, y = dim2, colour = Condition, shape = mixture_short)) +
        geom_point(size = 2.4, alpha = 0.9) +
        labs(x = sprintf("MDS 1 (%.0f%%)", 100 * coords$var_explained_dim1[1]),
             y = sprintf("MDS 2 (%.0f%%)", 100 * coords$var_explained_dim2[1]),
             shape = "mixture", colour = "condition",
             title = sprintf("%s: MDS of %d complete proteins, %d channels", label, nrow(complete), ncol(complete))) +
        theme_bw(base_size = 10)
    if (facet) p <- p + facet_wrap(~ mixture_short)
    save_fig(p, file.path(fig_dir, paste0(dataset, "_mds.png")), if (facet) 10 else 7, 4.5)
    invisible(coords)
}

## ---- Main preprocessing ---------------------------------------------------

preprocess_one <- function(dataset) {
    message("== ", dataset, " ==")
    filter_log <<- list()
    is_spike <- startsWith(dataset, "spikein1")
    qf <- readRDS(file.path(qf_dir, paste0(dataset, "_raw.rds")))
    runs <- names(qf)
    n0 <- n_psms(qf, runs)

    ## Memory: every QFeatures operation copies the object, and the raw
    ## rowData carries 34 columns for 300 000 PSMs. Keeping only what the
    ## filters and models need brought the peak from 5.2 GB to under the
    ## 4 GB ceiling on the first run of this stage.
    slim <- c("Annotated.Sequence", "Modifications", "Charge", "Protein.Accessions",
              "Master.Protein.Descriptions", "Marked.as", "PSM.Ambiguity", "Rank",
              "Search.Engine.Rank", "Confidence", "Isolation.Interference....",
              "Average.Reporter.S.N", "First.Scan")
    for (i in runs) rowData(qf[[i]]) <- rowData(qf[[i]])[, intersect(slim, colnames(rowData(qf[[i]]))), drop = FALSE]
    mem_checkpoint(timer, paste(dataset, "loaded and slimmed"))

    ## Ion identifier: sequence, modification string and charge. The
    ## spike-in adds a ups tag so that a light (spike-in) and a heavy
    ## (endogenous SILAC HeLa) version of the same peptide stay apart.
    for (i in runs) {
        rd <- rowData(qf[[i]])
        rd$ionID <- paste(rd$Annotated.Sequence, rd$Modifications, rd$Charge, sep = "_")
        if (is_spike) rd$ionID <- paste0(rd$ionID, ifelse(grepl("ups", rd$Protein.Accessions), "_ups", ""))
        rd$psmID <- paste(i, rd$First.Scan, rd$Rank, sep = "_")
        rowData(qf[[i]]) <- rd
    }

    ## -- Filters that come from the identification, in the vignette's order --

    if (is_spike) {
        ## Two Mascot nodes searched every spectrum, once against SwissProt
        ## with SILAC heavy labels as static modifications and once against
        ## the Sigma UPS sequences without them. A spectrum matched by both
        ## nodes appears twice with byte-identical reporter intensities.
        ## Keeping both would enter one scan's ten reporter ions into two
        ## proteins' summaries, doubling that scan's weight, and keeping
        ## one at random would assign the scan to a database by coin toss.
        ## Both copies go, as in the msqrob2TMT paper.
        n_before <- n_psms(qf, runs)
        for (i in runs) {
            x <- assay(qf[[i]])
            dup <- duplicated(x) | duplicated(x, fromLast = TRUE)
            rowData(qf[[i]])$duplicatedQuant <- as.vector(dup)
        }
        qf <- filterFeatures(qf, ~ !duplicatedQuant, i = runs)
        log_filter("duplicated reporter rows (two search nodes)", n_before, n_psms(qf, runs),
                   "identical reporter intensity vector within a run",
                   "two identifications of one scan; both dropped rather than double counting the scan")

        ## Ambiguous UPS status: the SwissProt accession says ups but the
        ## Sigma UPS node did not mark the PSM, or the reverse. Such a PSM
        ## could be endogenous HeLa protein or spike-in; the ground truth
        ## for it is unknown, so it cannot be scored as true or false.
        n_before <- n_psms(qf, runs)
        for (i in runs) {
            rd <- rowData(qf[[i]])
            is_ups <- grepl("ups", rd$Protein.Accessions)
            marked <- grepl("UPS", rd$Marked.as)
            rd$upsStatus <- ifelse(is_ups & marked, "ups", ifelse(!is_ups & !marked, "hela", "ambiguous"))
            rowData(qf[[i]]) <- rd
        }
        qf <- filterFeatures(qf, ~ upsStatus != "ambiguous", i = runs)
        log_filter("ambiguous UPS status", n_before, n_psms(qf, runs),
                   "accession says ups xor Sigma UPS node marked the PSM",
                   "no ground truth for these PSMs, so they cannot be scored")
    }

    n_before <- n_psms(qf, runs)
    qf <- filterFeatures(qf, ~ Protein.Accessions != "", i = runs)
    log_filter("failed protein inference", n_before, n_psms(qf, runs), "empty Protein.Accessions", "nothing to summarise into")

    n_before <- n_psms(qf, runs)
    qf <- filterFeatures(qf, ~ !grepl(";", Protein.Accessions), i = runs)
    log_filter("shared peptides (protein groups)", n_before, n_psms(qf, runs), "';' in Protein.Accessions",
               "a peptide shared by several proteins cannot be attributed; for the spike-in a group mixing ups and HeLa has no single truth")

    n_before <- n_psms(qf, runs)
    qf <- filterFeatures(qf, ~ PSM.Ambiguity != "Rejected", i = runs)
    log_filter("PSMs rejected by PD consensus", n_before, n_psms(qf, runs), "PSM.Ambiguity == Rejected",
               "PD's own consensus step rejected the match; identification quality, not abundance")

    n_before <- n_psms(qf, runs)
    qf <- filterFeatures(qf, ~ Rank == 1 & Search.Engine.Rank == 1, i = runs)
    log_filter("lower-ranked matches to a spectrum", n_before, n_psms(qf, runs), "Rank == 1 and Search.Engine.Rank == 1",
               "a second-ranked peptide for the same scan shares that scan's reporter ions")

    ## Decoys and contaminants: reported, not filtered, because there is
    ## nothing to filter. The export is post-FDR (Confidence is High for
    ## every row) and no contaminant database was searched.
    conf_tab <- table(unlist(lapply(runs, function(i) rowData(qf[[i]])$Confidence)))
    log_filter("decoys", n_psms(qf, runs), n_psms(qf, runs),
               paste0("Confidence classes present: ", paste(names(conf_tab), conf_tab, sep = "=", collapse = ", ")),
               "export is post-Percolator at q < 0.01; no decoy rows are present to remove")
    log_filter("contaminants", n_psms(qf, runs), n_psms(qf, runs), "no contaminant flag in the export",
               "the searches used SwissProt (plus Sigma UPS) without a contaminant database; nothing to remove")

    ## Ions that map to different proteins in different runs would be
    ## summarised into two proteins.
    n_before <- n_psms(qf, runs)
    rd_all <- stacked_rowdata(qf, runs, c("ionID", "Protein.Accessions"))
    ## unique() then count, not uniqueN() by group: data.table's grouped
    ## uniqueN allocated over 5 GB on this table, the single largest
    ## memory spike of the whole pipeline before it was replaced.
    n_map <- unique(rd_all[, .(ionID, Protein.Accessions)])[, .(nProtsMapped = .N), by = ionID]
    rd_all[n_map, nProtsMapped := i.nProtsMapped, on = "ionID"]
    qf <- set_rowdata_column(qf, runs, rd_all, "nProtsMapped")
    qf <- filterFeatures(qf, ~ nProtsMapped == 1, i = runs)
    log_filter("ions mapped to different proteins across runs", n_before, n_psms(qf, runs), "nProtsMapped == 1",
               "one ion, one protein, in every run")

    if (is_spike) {
        ## A protein seen in one run of fifteen cannot support a mixture
        ## or run effect and is, in the vignette's words, not trustworthy.
        n_before <- n_psms(qf, runs)
        prot_runs <- unique(rd_all[, .(Protein.Accessions, assay)])[, .(n_runs = .N), by = Protein.Accessions]
        one_run <- prot_runs[n_runs == 1, Protein.Accessions]
        qf <- filterFeatures(qf, ~ !Protein.Accessions %in% one_run, i = runs)
        log_filter("proteins seen in a single run", n_before, n_psms(qf, runs),
                   sprintf("%d proteins identified in 1 of %d runs", length(one_run), length(runs)),
                   "no replicate information at all for these proteins")
    }

    ## -- Missing values ---------------------------------------------------
    ##
    ## A zero reporter intensity is not a measured zero. Proteome
    ## Discoverer writes 0 when it integrated nothing at that reporter
    ## m/z, which is a detection failure of one channel in a scan whose
    ## other channels were quantified. Encoding it as NA keeps the
    ## log-transform finite and stops a run's normalisation median from
    ## being pulled by zeros.
    ##
    ## Where missingness sits changes under TMT. In label-free data a
    ## peptide unobserved in one sample is one missing cell, and its
    ## absence is informative because low abundance in that sample is
    ## the likeliest reason. Here every channel of an acquired scan gets
    ## a value in one integration, so a peptide is either quantified for
    ## all ten samples of a run or, when the precursor was never selected
    ## in that run, missing for all ten at once. Missingness moves from
    ## the sample level to the run level: it depends on whether the
    ## instrument picked the precursor in that injection, not on which
    ## sample carried it. The within-run zeros that remain are the
    ## channel-level exceptions and are counted below.
    miss_before <- rbindlist(lapply(runs, function(i) {
        x <- assay(qf[[i]])
        data.table(run = i, n_psms = nrow(x),
                   channel_values_zero = sum(x == 0, na.rm = TRUE),
                   channel_values_na = sum(is.na(x)),
                   psms_all_zero = sum(rowSums(x == 0 | is.na(x)) == ncol(x)),
                   psms_any_zero = sum(rowSums(x == 0 | is.na(x)) > 0))
    }))
    qf <- zeroIsNA(qf, runs)
    mem_checkpoint(timer, "zeroIsNA")

    ## Run-level missingness: in how many runs is each ion observed?
    ion_runs <- unique(stacked_rowdata(qf, runs, "ionID")[, .(ionID, assay)])[, .(n_runs = .N), by = ionID]
    run_level <- ion_runs[, .N, by = n_runs][order(n_runs)]
    setnames(run_level, "N", "n_ions")
    run_level[, dataset := dataset]
    write_tsv(run_level, file.path(res_dir, paste0(dataset, "_ion_runs_observed.tsv")))

    ## One PSM per ion per run: the PSM with the largest summed reporter
    ## intensity, as in the vignette and in MSstatsTMT. Repeated
    ## fragmentation of the same precursor in a run re-measures the same
    ## ten samples; keeping every spectrum would treat them as new
    ## observations of the ion.
    n_before <- n_psms(qf, runs)
    for (i in runs) {
        rd <- as.data.table(as.data.frame(rowData(qf[[i]])[, c("ionID"), drop = FALSE]))
        rd[, psmSum := rowSums(assay(qf[[i]]), na.rm = TRUE)]
        rd[, psmRank := frank(-psmSum, ties.method = "first"), by = ionID]
        rowData(qf[[i]])$psmSum <- rd$psmSum
        rowData(qf[[i]])$psmRank <- rd$psmRank
    }
    qf <- filterFeatures(qf, ~ psmRank == 1, i = runs)
    log_filter("repeat spectra of the same ion in a run", n_before, n_psms(qf, runs), "highest summed reporter intensity per ion per run",
               "one observation of each ion per run; identification-side choice")

    ## Missingness threshold: at most half (spike-in) or 70 per cent
    ## (mouse) of a PSM's channels missing. Both thresholds are the ones
    ## the msqrob2TMT vignettes use for these datasets and both are
    ## arbitrary; the mouse threshold is looser because Huang et al.
    ## kept spectra with reporter ions in at least 3 of 10 channels.
    pna <- if (is_spike) 0.5 else 0.7
    n_before <- n_psms(qf, runs)
    qf <- filterNA(qf, i = runs, pNA = pna)
    log_filter("PSMs with too many missing channels", n_before, n_psms(qf, runs), sprintf("more than %.0f%% of channels NA", 100 * pna),
               "completeness of the measurement, not its value; threshold follows the vignette and is arbitrary")

    miss_after <- rbindlist(lapply(runs, function(i) {
        x <- assay(qf[[i]])
        data.table(run = i, n_psms_after_filters = nrow(x),
                   channel_values_na_after = sum(is.na(x)),
                   frac_values_na_after = round(mean(is.na(x)), 4))
    }))
    mem_checkpoint(timer, paste(dataset, "filters done"))
    miss <- merge(miss_before, miss_after, by = "run")
    miss[, dataset := dataset]
    write_tsv(miss, file.path(res_dir, paste0(dataset, "_missingness.tsv")))

    ## -- Log2 and normalisation ------------------------------------------
    ##
    ## Both transformations are applied in place with replaceAssay rather
    ## than through logTransform() and normalize(), which each add a full
    ## copy of every run assay to the object. The arithmetic is the same
    ## (log2, then subtract each column's median) and the pre-centering
    ## medians are saved so that the un-centred values can be rebuilt
    ## for the QC figure and for the mixture-scope summarisation.
    for (i in runs) {
        se <- qf[[i]]
        assay(se) <- log2(assay(se))
        qf <- replaceAssay(qf, se, i)
    }

    ## Median centering per channel within each run. The assumption is
    ## equal loading: that every channel received the same total peptide
    ## amount, so any difference in the channel medians is technical
    ## (labelling efficiency, pipetting, reporter isotope purity) and not
    ## biological. In the spike-in that holds by construction: 50 ug of
    ## HeLa per channel and the UPS1 proteins are under one per cent of
    ## it. In the mouse data every channel is a different animal, and the
    ## assumption becomes that a high fat diet does not change the
    ## median adipose protein. Most proteins not changing is the usual
    ## justification, and the mixture-level MDS below is the check.
    norm_factors <- rbindlist(lapply(runs, function(i) {
        x <- assay(qf[[i]])
        data.table(run = i, sample = colnames(x),
                   median_log2_before = apply(x, 2, median, na.rm = TRUE))
    }))
    for (i in runs) {
        se <- qf[[i]]
        med <- norm_factors[run == i][match(colnames(se), sample), median_log2_before]
        assay(se) <- sweep(assay(se), 2, med, "-")
        qf <- replaceAssay(qf, se, i)
    }
    norm_sets <- runs
    col_median <- function(samples) norm_factors[match(samples, sample), median_log2_before]
    norm_factors[, dataset := dataset]
    write_tsv(norm_factors, file.path(res_dir, paste0(dataset, "_normalisation_factors.tsv")))

    ## Density figure: every channel of one mixture's runs, before and after.
    show_runs <- if (is_spike) grep("Mixture1_", runs, value = TRUE) else grep("^A-J_", runs, value = TRUE)[1:3]
    dens <- rbindlist(lapply(show_runs, function(r) {
        cd <- colData(qf)[colnames(qf[[r]]), ]
        x <- assay(qf[[r]])
        x_log <- sweep(x, 2, col_median(colnames(x)), "+")
        long <- function(m, stage) data.table(run = r, stage = stage, sample = rep(colnames(m), each = nrow(m)), value = as.vector(m))
        rbind(long(x_log, "log2 intensity"), long(x, "after median centering"))
    }))
    dens <- dens[!is.na(value)]
    dens[, channel := colData(qf)[as.character(sample), "Channel"]]
    dens[, condition := colData(qf)[as.character(sample), "Condition"]]
    dens[, stage := factor(stage, levels = c("log2 intensity", "after median centering"))]
    p_dens <- ggplot(dens, aes(x = value, colour = condition, group = sample)) +
        geom_density(linewidth = 0.4) +
        facet_grid(run ~ stage, scales = "free_x") +
        labs(x = "log2 reporter intensity", y = "density", colour = "condition",
             title = sprintf("%s: per-channel densities before and after median centering", dataset)) +
        theme_bw(base_size = 10) + theme(legend.position = "bottom")
    save_fig(p_dens, file.path(fig_dir, paste0(dataset, "_normalisation_densities.png")), 9, 2 + 1.6 * length(show_runs))
    rm(dens); mem_checkpoint(timer, paste(dataset, "normalised, density figure"))

    ## -- Summarisation within run -----------------------------------------
    ##
    ## robustSummary fits, for one protein in one run, y = feature + sample
    ## by M-estimation, and reports the sample effects as the protein's
    ## per-channel values. Peptides with wildly different intensities
    ## contribute through their offsets, not their raw means, and an
    ## outlying PSM is down-weighted rather than averaged in. The paper
    ## used median polish, which does the same job without weights; the
    ## robust fit is what stage 1 used and what MsCoreUtils recommends.
    message("  summarising ", length(norm_sets), " runs with robustSummary")
    qf <- aggregateFeatures(qf, i = norm_sets, fcol = "Protein.Accessions",
                            name = paste0(runs, "_proteins"), fun = MsCoreUtils::robustSummary)
    prot_sets <- paste0(runs, "_proteins")
    mem_checkpoint(timer, paste(dataset, "summarised"))

    ## Rows of the ion assays are named by ionID so that joining matches
    ## the same ion across runs.
    for (i in norm_sets) {
        se <- qf[[i]]; rownames(se) <- rowData(se)$ionID; qf <- replaceAssay(qf, se, i)
    }
    qf <- joinAssays(qf, i = norm_sets, name = "ions_norm")
    qf <- joinAssays(qf, i = prot_sets, name = "proteins")
    message(sprintf("  joined: %d ions x %d samples; %d proteins x %d samples",
                    nrow(qf[["ions_norm"]]), ncol(qf[["ions_norm"]]), nrow(qf[["proteins"]]), ncol(qf[["proteins"]])))
    mem_checkpoint(timer, paste(dataset, "joined"))

    ## -- Summarisation example figure -------------------------------------
    ## One protein with several ions, shown as its PSM-level values and
    ## as the summarised protein values across one mixture.
    cd <- as.data.table(as.data.frame(colData(qf)), keep.rownames = "sample")
    ex_mixture <- if (is_spike) "Mixture1" else grep("A-J", unique(cd$Mixture), value = TRUE)[1]
    ion_rd <- as.data.table(as.data.frame(rowData(qf[["ions_norm"]])[, c("Protein.Accessions"), drop = FALSE]), keep.rownames = "ionID")
    ion_counts <- ion_rd[, .N, by = Protein.Accessions]
    cand <- ion_counts[N >= 6 & N <= 12]
    if (is_spike) cand <- cand[grepl("ups", Protein.Accessions)]
    ex_prot <- cand$Protein.Accessions[which.max(cand$N)]
    ions_ex <- as.data.table(assay(qf[["ions_norm"]])[ion_rd[Protein.Accessions == ex_prot, ionID], , drop = FALSE], keep.rownames = "ionID") |>
        melt(id.vars = "ionID", variable.name = "sample", value.name = "value")
    ions_ex <- merge(ions_ex[!is.na(value)], cd, by = "sample")[Mixture == ex_mixture]
    prot_ex <- data.table(sample = colnames(qf[["proteins"]]), value = assay(qf[["proteins"]])[ex_prot, ])
    prot_ex <- merge(prot_ex[!is.na(value)], cd, by = "sample")[Mixture == ex_mixture]
    ions_ex[, sample := factor(sample, levels = unique(sample[order(run_short, Condition)]))]
    prot_ex[, sample := factor(sample, levels = levels(ions_ex$sample))]
    p_ex1 <- ggplot(ions_ex, aes(x = sample, y = value)) +
        geom_line(aes(group = ionID), linewidth = 0.2, colour = "grey60") +
        geom_point(aes(colour = Condition), size = 1.2) +
        facet_grid(~ run_short, scales = "free_x") +
        labs(y = "log2 normalised PSM intensity", x = NULL, title = sprintf("%s: %d ions of %s in %s", dataset, length(unique(ions_ex$ionID)), ex_prot, ex_mixture)) +
        theme_bw(base_size = 9) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
    p_ex2 <- ggplot(prot_ex, aes(x = sample, y = value, colour = Condition)) +
        geom_point(size = 2) +
        facet_grid(~ run_short, scales = "free_x") +
        labs(y = "robustSummary protein value", x = "channel", title = "after summarisation: one value per channel per run") +
        theme_bw(base_size = 9) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
    save_fig(p_ex1 / p_ex2, file.path(fig_dir, paste0(dataset, "_summarisation_example.png")), 9, 6.5)
    mem_checkpoint(timer, paste(dataset, "example figure"))

    ## -- MDS on summarised proteins, reference channels dropped for display --
    ## Spike-in: run-level values (150 channels). Mouse: the mixture-scope
    ## values built below, because a fractionated protein is almost never
    ## complete across all 27 fraction-runs (3 were, on the first attempt).
    if (is_spike) mds_figure(assay(qf[["proteins"]]), cd, dataset, facet = FALSE)

    fl <- rbindlist(filter_log)
    fl[, dataset := dataset]
    write_tsv(fl, file.path(res_dir, paste0(dataset, "_filter_log.tsv")))

    ## Drop the per-run intermediate assays before saving: the joined
    ## sets are what the models use, and the raw object is on disk.
    keep_sets <- c("ions_norm", "proteins")
    qf_small <- qf[, , keep_sets]
    saveRDS(qf_small, file.path(qf_dir, paste0(dataset, "_preprocessed.rds")))
    message("  saved ", dataset, "_preprocessed.rds")
    mem_checkpoint(timer, paste(dataset, "saved"))

    ## -- Mouse: summarisation within mixture ----------------------------------
    ##
    ## Fractionation scatters a protein's PSMs for one mixture across nine
    ## runs of the same physical sample pool. Summarising within a run
    ## gives nine protein values per channel per mixture and asks the
    ## model to carry a run effect; summarising within the mixture pools
    ## all fractions first and gives one value per channel, which is what
    ## MSstatsTMT does. Both are built here so that 07 can compare them.
    if (dataset == "mouse") {
        message("  building mixture-scope assays")
        mixes <- unique(cd$mixture_short)
        se_list <- list()
        for (mx in mixes) {
            mruns <- unique(cd[mixture_short == mx, run_short])
            mats <- list(); rds <- list()
            for (r in mruns) {
                se <- qf[[r]]
                x <- sweep(assay(se), 2, col_median(colnames(se)), "+")
                colnames(x) <- colData(qf)[colnames(x), "Channel"]
                rownames(x) <- paste(r, rownames(x), sep = "|")
                rd <- as.data.table(as.data.frame(rowData(se)[, c("ionID", "Protein.Accessions", "Master.Protein.Descriptions", "psmSum")]))
                rd[, run_short := r]
                mats[[r]] <- x; rds[[r]] <- rd
            }
            x <- do.call(rbind, mats)
            rd <- rbindlist(rds)
            ## One PSM per ion per mixture, the highest summed intensity
            ## across fractions, following the paper's mixture script.
            rd[, rank_mix := frank(-psmSum, ties.method = "first"), by = ionID]
            keep_rows <- which(rd$rank_mix == 1)
            x <- x[keep_rows, , drop = FALSE]; rd <- rd[keep_rows]
            cdm <- unique(cd[mixture_short == mx, .(Channel, Condition, Mixture, mixture_short, BioReplicate, Diet, Duration)])
            cdm <- as.data.frame(cdm); rownames(cdm) <- paste(mx, cdm$Channel, sep = "_")
            colnames(x) <- paste(mx, colnames(x), sep = "_")
            cdm <- cdm[colnames(x), ]
            rdf <- as.data.frame(rd); rownames(rdf) <- rownames(x)
            se_list[[mx]] <- SummarizedExperiment(assays = list(x), rowData = rdf, colData = cdm)
        }
        ## QFeatures maps colData rows to assay columns by name; rbind of
        ## named list elements would prefix the rownames, hence unname.
        cd_all <- do.call(rbind, unname(lapply(se_list, function(s) as.data.frame(colData(s)))))
        qfm <- QFeatures(se_list, colData = cd_all)
        names(qfm) <- paste0(mixes, "_log")
        qfm <- normalize(qfm, i = names(qfm), name = paste0(mixes, "_norm"), method = "center.median")
        qfm <- aggregateFeatures(qfm, i = paste0(mixes, "_norm"), fcol = "Protein.Accessions",
                                 name = paste0(mixes, "_proteins"), fun = MsCoreUtils::robustSummary)
        qfm <- joinAssays(qfm, i = paste0(mixes, "_proteins"), name = "proteins")
        message(sprintf("  mixture scope: %d proteins x %d channels", nrow(qfm[["proteins"]]), ncol(qfm[["proteins"]])))
        cdm_all <- as.data.table(as.data.frame(colData(qfm)), keep.rownames = "sample")
        mds_figure(assay(qfm[["proteins"]]), cdm_all, "mouse", facet = FALSE, label = "mouse, mixture-scope summaries")
        saveRDS(qfm[, , "proteins"], file.path(qf_dir, "mouse_mixture_preprocessed.rds"))
        write_tsv(data.table(scope = c("fraction_run", "mixture"),
                             n_proteins = c(nrow(qf[["proteins"]]), nrow(qfm[["proteins"]])),
                             n_columns = c(ncol(qf[["proteins"]]), ncol(qfm[["proteins"]])),
                             psms_used = c(nrow(qf[["ions_norm"]]), sum(vapply(mixes, function(mx) nrow(qfm[[paste0(mx, "_norm")]]), integer(1))))),
                  file.path(res_dir, "mouse_summarisation_scopes.tsv"))
    }
    invisible(NULL)
}

for (d in which_ds) { preprocess_one(d); invisible(gc()) }

stage_end(timer, note = paste("datasets:", paste(which_ds, collapse = ",")))
stage_mark_done(STAGE, conf)
