## Stage 06: how much does ratio compression cost, and does SPS-MS3
## remove it?
##
## Co-isolated precursors contribute reporter ions to every channel of a
## scan, which pulls every measured fold change toward zero. On the
## spike-in the true fold changes are fixed by the pipetting, so the
## compression is a number: the slope of observed against expected
## log2 fold change across the six contrasts and all UPS1 proteins.
##
## The same workflow (protein level, reference channels dropped, mixed
## model, M-estimation) was fitted in 05 on both acquisitions of the same
## five mixtures: the SPS-MS3 export that the msqrob2TMT authors
## timestamped on Zenodo, and the MS2-only export from MassIVE
## MSV000084266. This stage compares them and, at the PSM level, relates
## the per-scan fold change to the isolation interference that Proteome
## Discoverer recorded for that scan.
##
## Usage: scripts/lib/rscript.sh scripts/06_compression.R [--force]
##
## Outputs:
##   results/spikein1_acquisition.tsv                  which table is which acquisition, and the evidence
##   results/spikein1_compression_by_contrast.tsv      medians, IQRs and bias per contrast and acquisition
##   results/spikein1_compression_slope.tsv            attenuation slopes
##   results/spikein1_compression_by_interference.tsv  PSM-level fold change against isolation interference
##   figures/spikein1_compression.png
##   figures/spikein1_compression_by_interference.png

suppressPackageStartupMessages({
    library("data.table")
    library("QFeatures")
    library("ggplot2")
    library("patchwork")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
source(file.path(root, "scripts", "lib", "common.R"))
conf <- read_conf(root)
opt <- parse_args()
STAGE <- "06_compression"
if (stage_already_done(STAGE, conf, force = opt$force)) quit(save = "no")
timer <- stage_begin(STAGE, conf)

res_dir <- file.path(root, "results")
fig_dir <- file.path(root, "figures")
qf_dir <- file.path(root, "results", "qfeatures")
expected <- spikein_expected_contrasts()

## ---- Which acquisition is which ---------------------------------------------
##
## Nothing inside a PD PSM table says MS2 or MS3: the "MS Order" column
## reports the identifying scan (MS2 in both). The assignment rests on
## provenance: the Zenodo file was copied from the MassIVE reanalysis
## container of MSV000084264, whose title reads "acquired using
## SPS-MS3", and its raw file names carry the date 161117; the MassIVE
## dataset MSV000084266, "acquired using MS2-only strategy", carries
## 161122 in its raw file names and the token MS2 in its export name.
acq <- data.table(
    dataset = c("spikein1", "spikein1_ms2"),
    acquisition = c("SPS-MS3", "MS2"),
    source = c("zenodo:14767905 spikein1_psms.txt (from MassIVE RMSV000000265, reanalysis of MSV000084264)",
               "massive:MSV000084266 quant/..._PD22_MS2_Intensity_01-(1)_PSMs.txt"),
    proteomexchange = c("PXD015258", "PXD015261"),
    raw_file_date = c("161117", "161122"),
    evidence = c("MassIVE MSV000084264 title: TMT10 controlled mixtures - SILAC HeLa UPS1, acquired using SPS-MS3; Huang et al. 2020 methods: SpikeIn-5mix-MS3 acquired using SPS",
                 "MassIVE MSV000084266 title: TMT10 controlled mixtures - SILAC HeLa UPS1, acquired using MS2-only strategy; export filename contains MS2")
)
write_tsv(acq, file.path(res_dir, "spikein1_acquisition.tsv"))

## ---- Protein-level compression from the benchmark results ---------------------

bench <- fread(file.path(res_dir, "spikein1_benchmark_results.tsv.gz"))
adopted <- bench[design == "balanced" & variant == "protein_drop_mixed_robust_noridge_msqrob2" & !is.na(logFC)]
adopted[, acquisition := ifelse(acquisition == "ms3", "SPS-MS3", "MS2")]
adopted[, is_ups := grepl("ups", protein)]
adopted <- merge(adopted, expected[, .(contrast, expected_log2fc, ratio)], by = "contrast")

by_contrast <- adopted[, .(
    n_ups = sum(is_ups), n_hela = sum(!is_ups),
    expected_log2fc = expected_log2fc[1],
    ups_median = median(logFC[is_ups]), ups_q25 = quantile(logFC[is_ups], 0.25), ups_q75 = quantile(logFC[is_ups], 0.75),
    ups_iqr = IQR(logFC[is_ups]),
    hela_median = median(logFC[!is_ups]), hela_iqr = IQR(logFC[!is_ups]),
    bias = median(logFC[is_ups]) - expected_log2fc[1],
    observed_over_expected = median(logFC[is_ups]) / expected_log2fc[1]
), by = .(acquisition, contrast)][order(acquisition, expected_log2fc)]
write_tsv(by_contrast, file.path(res_dir, "spikein1_compression_by_contrast.tsv"))
print(by_contrast[, .(acquisition, contrast, expected_log2fc, ups_median, ups_iqr, hela_median, bias)])

## Attenuation slope: every UPS1 protein in every contrast is one point.
## Fitted with an intercept (which should be near zero) and through the
## origin; the with-intercept slope is the headline number. The
## comparison is also made on the UPS1 proteins quantified under both
## acquisitions, so that a different protein set cannot explain it.
ups <- adopted[is_ups == TRUE]
shared <- Reduce(intersect, split(ups$protein, ups$acquisition))
slope_one <- function(d, label) {
    f1 <- lm(logFC ~ expected_log2fc, data = d)
    f0 <- lm(logFC ~ 0 + expected_log2fc, data = d)
    data.table(protein_set = label, n_points = nrow(d), n_proteins = length(unique(d$protein)),
               slope = coef(f1)[["expected_log2fc"]], slope_se = summary(f1)$coefficients["expected_log2fc", 2],
               intercept = coef(f1)[["(Intercept)"]], r_squared = summary(f1)$r.squared,
               slope_through_origin = coef(f0)[["expected_log2fc"]])
}
slopes <- rbind(
    ups[, slope_one(.SD, "all UPS1 proteins quantified"), by = acquisition],
    ups[protein %in% shared, slope_one(.SD, "UPS1 proteins quantified in both acquisitions"), by = acquisition]
)
write_tsv(slopes, file.path(res_dir, "spikein1_compression_slope.tsv"))
print(slopes)

## Per-protein slopes, to show the spread behind the single number.
per_protein <- ups[, .(n_contrasts = .N, slope = if (.N >= 3) coef(lm(logFC ~ 0 + expected_log2fc))[[1]] else NA_real_),
                   by = .(acquisition, protein)]
write_tsv(per_protein, file.path(res_dir, "spikein1_compression_per_protein.tsv"))

## ---- Figure: observed against expected ---------------------------------------

lab <- slopes[protein_set == "all UPS1 proteins quantified"]
lab[, text := sprintf("slope %.2f (n = %d proteins)", slope, n_proteins)]
p1 <- ggplot(ups, aes(x = expected_log2fc, y = logFC)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_jitter(width = 0.04, height = 0, alpha = 0.35, size = 1, colour = "firebrick") +
    stat_summary(fun = median, geom = "point", size = 2.6, shape = 21, fill = "white", colour = "black") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
    geom_text(data = lab, aes(x = 0.2, y = 2.85, label = text), hjust = 0, size = 3.2) +
    facet_wrap(~ acquisition) +
    scale_x_continuous(breaks = round(expected$expected_log2fc, 2)) +
    labs(x = "expected log2 fold change (from fmol)", y = "observed log2 fold change, UPS1 proteins",
         title = "Ratio compression on the spike-in: each dot is one UPS1 protein in one contrast",
         subtitle = "dashed: identity; white points: median per contrast; black line: fitted slope") +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
hela <- adopted[is_ups == FALSE]
p2 <- ggplot(hela, aes(x = factor(round(expected_log2fc, 2)), y = logFC)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_boxplot(outlier.size = 0.4, fill = "grey90") +
    facet_wrap(~ acquisition) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(x = "contrast, as expected UPS1 log2 fold change", y = "observed log2 fold change, HeLa proteins",
         title = "The background: HeLa proteins should sit at zero in every contrast") +
    theme_bw(base_size = 10)
save_fig(p1 / p2 + plot_layout(heights = c(3, 2)), file.path(fig_dir, "spikein1_compression.png"), 9, 8.5)

## ---- PSM-level: fold change against isolation interference -------------------
##
## The mechanism, measured scan by scan. For every UPS1 PSM in every run,
## the log2 ratio of the two 500 fmol channels to the two 62.5 fmol
## channels (expected 3.00) is computed from the normalised reporter
## values, then binned by the isolation interference that Proteome
## Discoverer computed for the precursor's isolation window.
psm_rows <- list()
for (ds in c("spikein1", "spikein1_ms2")) {
    ## The joined ion assay of stage 04 has lost the scan-specific
    ## columns (joinAssays drops rowData that differs between runs), so
    ## this goes back to the raw per-run object and re-applies the
    ## identification filters of stage 04 that matter here: UPS status
    ## confirmed by both search nodes, no duplicated reporter vectors,
    ## rank 1, zeros as missing. The ratio is computed from raw
    ## intensities: within one scan, per-channel normalisation is a
    ## constant offset that median centering removes and the ratio does
    ## not need.
    qf <- readRDS(file.path(qf_dir, paste0(ds, "_raw.rds")))
    cd <- as.data.table(as.data.frame(colData(qf)), keep.rownames = "sample")
    for (r in names(qf)) {
        se <- qf[[r]]
        rd <- rowData(se)
        x <- assay(se); x[x == 0] <- NA
        dup <- duplicated(x) | duplicated(x, fromLast = TRUE)
        is_ups_acc <- grepl("ups", rd$Protein.Accessions)
        marked <- grepl("UPS", rd$Marked.as)
        keep <- !dup & rd$Rank == 1 & !grepl(";", rd$Protein.Accessions) & rd$Protein.Accessions != "" &
            ((is_ups_acc & marked) | (!is_ups_acc & !marked))
        h <- cd[Run == r & Condition == "1", sample]; l <- cd[Run == r & Condition == "0.125", sample]
        fc <- log2(rowMeans(x[, h, drop = FALSE], na.rm = TRUE)) - log2(rowMeans(x[, l, drop = FALSE], na.rm = TRUE))
        psm_rows[[length(psm_rows) + 1L]] <- data.table(
            dataset = ds, run = r, is_ups = is_ups_acc, interference = rd$Isolation.Interference....,
            reporter_sn = rd$Average.Reporter.S.N, log2fc_1_vs_0.125 = fc)[keep]
    }
    rm(qf); invisible(gc())
}
psm <- rbindlist(psm_rows)[is.finite(log2fc_1_vs_0.125) & !is.na(interference)]
psm[, acquisition := ifelse(dataset == "spikein1", "SPS-MS3", "MS2")]
psm[, interference_bin := cut(interference, c(-Inf, 10, 20, 30, 40, 50, 75, 100), right = TRUE,
                              labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-75", "75-100"))]
by_int <- psm[, .(n_psms = .N,
                  median_log2fc = median(log2fc_1_vs_0.125), q25 = quantile(log2fc_1_vs_0.125, 0.25),
                  q75 = quantile(log2fc_1_vs_0.125, 0.75), expected = 3),
              by = .(acquisition, is_ups, interference_bin)][order(acquisition, is_ups, interference_bin)]
write_tsv(by_int, file.path(res_dir, "spikein1_compression_by_interference.tsv"))
print(by_int[is_ups == TRUE])
overall_int <- psm[, .(median_interference = median(interference), mean_interference = mean(interference),
                       frac_over_30 = mean(interference > 30), n = .N), by = .(acquisition, is_ups)]
write_tsv(overall_int, file.path(res_dir, "spikein1_interference_summary.tsv"))

p3 <- ggplot(by_int[is_ups == TRUE], aes(x = interference_bin, y = median_log2fc, colour = acquisition, group = acquisition)) +
    geom_hline(yintercept = 3, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.2, position = position_dodge(0.4)) +
    geom_line(position = position_dodge(0.4)) + geom_point(size = 2.2, position = position_dodge(0.4)) +
    geom_text(aes(label = n_psms, y = q25 - 0.12), size = 2.6, position = position_dodge(0.4), show.legend = FALSE) +
    labs(x = "isolation interference of the identifying scan (%)", y = "median log2 (500 fmol / 62.5 fmol), UPS1 PSMs",
         title = "Compression grows with isolation interference; expected value 3.0 (dashed)",
         subtitle = "median and interquartile range per bin; numbers are PSMs") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
save_fig(p3, file.path(fig_dir, "spikein1_compression_by_interference.png"), 7.5, 4.8)

stage_end(timer, note = sprintf("slopes: %s", paste(sprintf("%s=%.3f", lab$acquisition, lab$slope), collapse = ", ")))
stage_mark_done(STAGE, conf)
