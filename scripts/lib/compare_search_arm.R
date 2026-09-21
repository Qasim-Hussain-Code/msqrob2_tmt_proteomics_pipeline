## Compare the Sage re-search of one spike-in mixture with the Proteome
## Discoverer export on the same three runs. Called by
## scripts/10_search_arm.sh; not run for the results in this repository.
##
## Comparisons:
##   1. PSMs per run at 1 per cent FDR (Sage spectrum_q; PD Confidence High)
##   2. Identification overlap on scan number, and sequence agreement
##      among scans both engines identified
##   3. Per-channel correlation of log2 reporter intensities on shared PSMs
##   4. Differential abundance from the same protein-level model on both
##      tables: ~ 0 + Condition + (1 | Run), a single mixture so no
##      mixture term, robust M-estimation, six contrasts

suppressPackageStartupMessages({
    library("data.table")
    library("QFeatures")
    library("msqrob2")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."), winslash = "/")
source(file.path(root, "scripts", "lib", "common.R"))
source(file.path(root, "scripts", "lib", "msqrob_helpers.R"))
conf <- read_conf(root)
opt <- parse_args()
mixture <- as.integer(if (is.null(opt$mixture)) 1 else opt$mixture)
out <- file.path(root, "results", "sage")
BiocParallel::register(BiocParallel::SerialParam())

sage <- fread(file.path(out, "results.sage.tsv"))
tmt <- fread(file.path(out, "tmt.tsv"))
sage[, run := sub("\\.mzML(\\.gz)?$", "", basename(filename))]
sage[, scan := as.integer(sub(".*scan=([0-9]+).*", "\\1", scannr))]
sage[, stripped := gsub("\\[[^]]*\\]", "", peptide)]
sage_ok <- sage[spectrum_q <= 0.01 & label == 1]

## PD export for the same runs.
qf <- readRDS(file.path(root, "results", "qfeatures", "spikein1_ms2_raw.rds"))
runs <- grep(sprintf("Mixture%d_", mixture), names(qf), value = TRUE)
pd <- rbindlist(lapply(runs, function(r) {
    rd <- rowData(qf[[r]])
    x <- assay(qf[[r]])
    colnames(x) <- sub("^.*Abundance\\.\\.", "", colnames(x))
    cbind(data.table(run_short = r, scan = rd$First.Scan, stripped = toupper(gsub("^\\[.\\]\\.|\\.\\[.\\]$", "", rd$Annotated.Sequence)),
                     protein = rd$Protein.Accessions, confidence = rd$Confidence, rank = rd$Rank), as.data.table(x))
}))
pd[, run := sub("^.*(Mixture[0-9]+_[0-9]+)$", "161122_SILAC_HeLa_UPS1_TMT10_MS2_\\1", run_short)]
pd_ok <- pd[confidence == "High" & rank == 1]

## 1. counts
counts <- merge(sage_ok[, .(sage_psms = .N), by = run], pd_ok[, .(pd_psms = .N), by = run], by = "run", all = TRUE)
write_tsv(counts, file.path(out, "psm_counts.tsv"))

## 2. overlap on scan number
ov <- merge(sage_ok[, .(run, scan, sage_seq = toupper(stripped))], pd_ok[, .(run, scan, pd_seq = stripped)], by = c("run", "scan"))
overlap <- data.table(
    scans_sage = nrow(sage_ok), scans_pd = nrow(pd_ok), scans_both = nrow(ov),
    sequence_agreement_among_shared = mean(ov$sage_seq == ov$pd_seq))
write_tsv(overlap, file.path(out, "identification_overlap.tsv"))

## 3. reporter correlation on shared PSMs, per channel
tmt[, run := sub("\\.mzML(\\.gz)?$", "", basename(filename))]
tmt[, scan := as.integer(sub(".*scan=([0-9]+).*", "\\1", scannr))]
chan_sage <- grep("^tmt_", names(tmt), value = TRUE)
chan_pd <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131")
sh <- merge(tmt[, c("run", "scan", chan_sage), with = FALSE], pd_ok[, c("run", "scan", chan_pd), with = FALSE], by = c("run", "scan"))
cors <- rbindlist(lapply(seq_along(chan_pd), function(k) {
    a <- log2(sh[[chan_sage[k]]]); b <- log2(sh[[chan_pd[k]]])
    ok <- is.finite(a) & is.finite(b)
    data.table(channel = chan_pd[k], n = sum(ok), pearson_log2 = cor(a[ok], b[ok]))
}))
write_tsv(cors, file.path(out, "reporter_correlation.tsv"))

## 4. differential abundance on both tables with the same model
ann <- fread(manifest_path(conf, "spikein1_ms2_annotations.csv"), colClasses = "character")
ann <- ann[grepl(sprintf("Mixture%d_", mixture), Run)]
build_qf <- function(dt, quant_cols, protein_col, engine) {
    dt <- copy(dt)
    setnames(dt, quant_cols, paste0("Abundance..", chan_pd))
    ann2 <- copy(ann); ann2[, runCol := Run]; ann2[, quantCols := paste0("Abundance..", Channel)]
    dt[, Spectrum.File := paste0(run, ".raw")]
    q <- readQFeatures(as.data.frame(dt), colData = as.data.frame(ann2), quantCols = unique(ann2$quantCols), runCol = "Spectrum.File", name = "psms", verbose = FALSE)
    q <- filterFeatures(q, ~ !grepl(";", .data[[protein_col]]) & .data[[protein_col]] != "")
    q <- zeroIsNA(q, names(q)); q <- filterNA(q, names(q), pNA = 0.5)
    q <- logTransform(q, names(q), name = paste0(names(q), "_log"), base = 2)
    q <- normalize(q, paste0(names(q)[seq_along(runs)], "_log"), name = paste0(names(q)[seq_along(runs)], "_norm"), method = "center.median")
    q <- aggregateFeatures(q, paste0(names(q)[seq_along(runs)], "_norm"), fcol = protein_col, name = paste0(names(q)[seq_along(runs)], "_prot"), fun = MsCoreUtils::robustSummary)
    q <- joinAssays(q, paste0(names(q)[seq_along(runs)], "_prot"), "proteins")
    q <- subsetByColData(q, q$Condition != "Norm")
    q
}
contr <- fread(file.path(root, "config", "contrasts_spikein1.tsv"))
params <- paste0("Condition", spikein_conditions()$condition)
q_pd <- build_qf(pd_ok[, c("run", "protein", chan_pd), with = FALSE], chan_pd, "protein", "pd")
sage_prot <- merge(sage_ok[, .(run, scan, protein = proteins)], tmt[, c("run", "scan", chan_sage), with = FALSE], by = c("run", "scan"))
sage_prot[, protein := sub("^sp\\|([^|]+)\\|.*$", "\\1", protein)]
q_sage <- build_qf(sage_prot, chan_sage, "protein", "sage")
res_pd <- fit_protein_level(q_pd, "proteins", ~ 0 + Condition + (1 | Run), contr$hypothesis, params, contrast_labels = contr$contrast)$results
res_sage <- fit_protein_level(q_sage, "proteins", ~ 0 + Condition + (1 | Run), contr$hypothesis, params, contrast_labels = contr$contrast)$results
expected <- spikein_expected_contrasts()
da <- rbind(benchmark_metrics(res_pd, expected)[, engine := "proteome_discoverer"],
            benchmark_metrics(res_sage, expected)[, engine := "sage"])
write_tsv(da, file.path(out, "differential_abundance_comparison.tsv"))
message("search arm comparison written to results/sage/")
