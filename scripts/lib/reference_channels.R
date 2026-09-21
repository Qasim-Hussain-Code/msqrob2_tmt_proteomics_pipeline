## Three treatments of the pooled reference channels, as functions, so
## that the benchmark in 05 selects one with a flag rather than a
## rewrite.
##
## Both spike-in mixtures and the mouse mixtures carry pooled reference
## channels (Condition == "Norm" in the annotation files). There are
## three defensible things to do with them:
##
##   "drop"  Discard them and let the mixed model's run and mixture
##           random effects absorb between-plex shifts. This is what the
##           msqrob2TMT paper does; the authors said they found no
##           benefit from reference normalisation and that dividing by a
##           noisy channel adds that channel's noise to every sample.
##   "keep"  Retain them as ordinary observations with their own
##           Condition level. They then contribute to the run and
##           mixture variance components (they are the only samples
##           present in every plex with the same composition) without
##           being consumed by a division.
##   "ratio" Express every channel as a log ratio to the mean of the
##           run's reference channels, feature by feature, then discard
##           the references. This is the classical TMT normalisation
##           (MSstatsTMT's reference_norm) and it makes intensities
##           comparable across plexes by construction, at the cost of
##           propagating the reference channels' measurement error into
##           every value and of removing the between-plex variation the
##           mixed model would otherwise estimate.
##
## Which one to adopt is an empirical question that the spike-in answers
## in 05_benchmark_workflows.R.

suppressPackageStartupMessages({
    library("QFeatures")
    library("SummarizedExperiment")
})

## Subtract the per-feature mean of the reference channels within each
## run assay. Operates on log2 data, so subtraction is a ratio. Features
## with no observed reference value in a run become NA for that run:
## a ratio to nothing is not a value, and silently keeping the raw
## intensity would mix two scales in one column.
ratio_to_reference <- function(qf, sets, ref_label = "Norm", cond_col = "Condition") {
    for (i in sets) {
        se <- qf[[i]]
        cd <- colData(qf)[colnames(se), , drop = FALSE]
        is_ref <- cd[[cond_col]] == ref_label
        if (!any(is_ref)) stop("set ", i, " has no reference channel labelled ", ref_label)
        m <- assay(se)
        ref_mean <- rowMeans(m[, is_ref, drop = FALSE], na.rm = TRUE)
        ref_mean[is.nan(ref_mean)] <- NA
        assay(se) <- m - ref_mean
        qf <- replaceAssay(qf, se, i)
    }
    qf
}

## Apply a treatment to a QFeatures object. `sets` are the assays on
## which the ratio is computed (the run-level assays the models will
## use); dropping is applied object-wide through the colData.
apply_reference_treatment <- function(qf, mode = c("drop", "keep", "ratio"), sets,
                                      ref_label = "Norm", cond_col = "Condition") {
    mode <- match.arg(mode)
    n_ref <- sum(colData(qf)[[cond_col]] == ref_label)
    if (mode == "keep") {
        return(list(qf = qf, note = sprintf("%d reference channels kept as Condition level '%s'", n_ref, ref_label)))
    }
    if (mode == "ratio") {
        qf <- ratio_to_reference(qf, sets, ref_label, cond_col)
    }
    qf <- subsetByColData(qf, colData(qf)[[cond_col]] != ref_label)
    list(qf = qf,
         note = if (mode == "drop") sprintf("%d reference channels dropped", n_ref)
                else sprintf("%d reference channels consumed as per-feature ratio denominators", n_ref))
}
