## Stage 03: read the PSM tables and the annotation files, join channel to
## sample through the annotations, and build one QFeatures object per
## dataset with one assay per MS run and the reporter channels as columns.
##
## Usage: scripts/lib/rscript.sh scripts/03_build_qfeatures.R [--dataset spikein1|spikein1_ms2|mouse|all] [--force]
##
## Outputs (results/qfeatures/ is not tracked; everything else is):
##   results/qfeatures/<dataset>_raw.rds        QFeatures, raw reporter intensities
##   results/<dataset>_psm_counts_raw.tsv       PSMs, ions and proteins per run
##   results/<dataset>_columns_kept.tsv         PSM table columns retained and why
##   results/mouse_run_structure.tsv            mixtures x fractions derived from the data

suppressPackageStartupMessages({
    library("data.table")
    library("QFeatures")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
source(file.path(root, "scripts", "lib", "common.R"))
conf <- read_conf(root)
opt <- parse_args()
STAGE <- "03_build_qfeatures"
if (stage_already_done(STAGE, conf, force = opt$force)) quit(save = "no")
timer <- stage_begin(STAGE, conf)

which_ds <- if (is.null(opt$dataset) || opt$dataset == "all") c("spikein1", "spikein1_ms2", "mouse") else opt$dataset
out_dir <- file.path(root, "results", "qfeatures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Columns retained from the Proteome Discoverer 2.2 PSM export. The
## tables carry 50 columns; the long free-text description columns and
## the mass columns are not needed downstream and dropping them at read
## time halves the memory footprint of the quarter-gigabyte spike-in
## table. Names are the check.names = TRUE forms. Every column kept is
## either an identifier the model needs, a filter criterion, or a
## quality covariate that the report plots.
keep_columns <- function(available) {
    wanted <- c(
        Confidence = "Proteome Discoverer confidence class; expected High only after FDR control",
        Identifying.Node = "Mascot node (SwissProt or Sigma UPS search) that produced the PSM",
        PSM.Ambiguity = "Unambiguous, Selected or Rejected by the PD consensus step",
        Annotated.Sequence = "peptide sequence, part of the ion identifier",
        Modifications = "modification string, part of the ion identifier",
        Marked.as = "spike-in only: whether the Sigma UPS node marked the PSM as UPS",
        X..Protein.Groups = "number of protein groups the PSM maps to",
        X..Proteins = "number of proteins the PSM maps to",
        Master.Protein.Accessions = "master protein of the group",
        Master.Protein.Descriptions = "description of the master protein, used for the keyword check",
        Protein.Accessions = "all proteins the peptide maps to; ';' separates a group, used as the protein identifier as in the vignette",
        Charge = "precursor charge, part of the ion identifier",
        Rank = "PSM rank for the spectrum",
        Search.Engine.Rank = "Mascot rank",
        MS.Order = "MS level of the identifying scan",
        Isolation.Interference.... = "per cent of the isolation window intensity not from the precursor",
        Average.Reporter.S.N = "mean reporter signal to noise",
        Ion.Inject.Time..ms. = "ion injection time of the quantification scan",
        Spectrum.File = "raw file name; the run identifier joined to the annotation",
        First.Scan = "scan number, part of the PSM identifier",
        Quan.Info = "PD quantification flags",
        Ions.Score = "Mascot ion score",
        Percolator.q.Value = "Percolator PSM q-value",
        Percolator.PEP = "Percolator posterior error probability"
    )
    abund <- grep("^Abundance\\.\\.", available, value = TRUE)
    present <- intersect(names(wanted), available)
    list(cols = c(present, abund),
         table = data.table(column = c(present, abund),
                            reason = c(unname(wanted[present]), rep("reporter ion intensity", length(abund)))))
}

## The annotation files map every channel of every run to a sample. A
## written design table beats parsing conditions or mixtures out of raw
## file names: file names are typed by the person operating the
## instrument, they change between acquisitions of the same design
## (161117 for MS3, 161122 for MS2 here), and nothing in a file name says
## which channel carried which sample. The annotation is the record of
## what was pipetted, and the join through it is the only place that
## information enters the analysis.
read_annotation <- function(path, dataset) {
    ann <- fread(path, colClasses = "character")
    needed <- c("Run", "Channel", "Condition", "Mixture", "TechRepMixture", "BioReplicate", "Fraction")
    miss <- setdiff(needed, names(ann))
    if (length(miss)) stop("annotation ", basename(path), " lacks columns: ", paste(miss, collapse = ", "))
    ann <- ann[, ..needed]
    if (dataset == "mouse") {
        ## Diet and Duration are encoded jointly in Condition (Short_LF and
        ## so on). Split them for the factorial model, and give the
        ## factors the reference levels the contrasts assume: LF and 8
        ## weeks (Short) as baselines, so DietHF is the high fat effect at
        ## 8 weeks and DietHF:DurationLong is how much it changes by 18.
        ann[, Duration := sub("_.*$", "", Condition)]
        ann[, Diet := sub("^.*_", "", Condition)]
        ann[, run_short := paste0(sub("^PAMI-[0-9]+_Mouse_", "", Mixture), "_",
                                 sub("^.*_([0-9]+)pctACN.*$", "F\\1", Fraction))]
        ann[, mixture_short := sub("^PAMI-[0-9]+_Mouse_", "", Mixture)]
    } else {
        ann[, run_short := sub("^.*(Mixture[0-9]+_[0-9]+)\\.raw$", "\\1", Run)]
        ann[, mixture_short := Mixture]
    }
    ann[, runCol := Run]
    ann[, quantCols := paste0("Abundance..", Channel)]
    ann
}

build_one <- function(dataset) {
    message("== ", dataset, " ==")
    psm_file <- manifest_path(conf, paste0(dataset, "_psms.txt"))
    ann_file <- manifest_path(conf, paste0(dataset, "_annotations.csv"))

    ## integer64 = "double": reporter intensities are large integers in
    ## the PD export and fread would otherwise read them as bit64
    ## integer64, which silently mangles arithmetic downstream. Same
    ## reason as in the label-free repository.
    ## fread's select works on the original header, so map the
    ## check.names forms back to the raw column names before selecting.
    header_raw <- names(fread(psm_file, nrows = 0, check.names = FALSE))
    header <- make.names(header_raw, unique = TRUE)
    kc <- keep_columns(header)
    psms <- fread(psm_file, check.names = TRUE, integer64 = "double",
                  select = header_raw[match(kc$cols, header)])
    message(sprintf("read %d PSMs x %d columns (kept %d of %d)", nrow(psms), ncol(psms), length(kc$cols), length(header)))
    write_tsv(kc$table, file.path(root, "results", paste0(dataset, "_columns_kept.tsv")))

    ann <- read_annotation(ann_file, if (startsWith(dataset, "spikein1")) "spikein1" else dataset)

    ## Join check, both directions. A run present in the PSM table but
    ## absent from the annotation would silently drop out of the
    ## QFeatures object; a run in the annotation without PSMs is a sign
    ## of a wrong file pairing.
    runs_psm <- unique(psms$Spectrum.File)
    runs_ann <- unique(ann$Run)
    if (!setequal(runs_psm, runs_ann)) {
        stop(sprintf("run names differ between PSM table and annotation for %s\n  only in PSMs: %s\n  only in annotation: %s",
                     dataset, paste(setdiff(runs_psm, runs_ann), collapse = ", "),
                     paste(setdiff(runs_ann, runs_psm), collapse = ", ")))
    }
    message(sprintf("%d runs in both the PSM table and the annotation", length(runs_psm)))

    ## Run structure derived from the data, not from any paper. For the
    ## mouse study the sources disagree on the fraction count (eight in
    ## the MSstatsTMT main text, nine in its supplement and in the
    ## msqrob2TMT repository); the annotation decides.
    if (dataset == "mouse") {
        rs <- ann[, .(n_channels = .N,
                      n_reference = sum(Condition == "Norm"),
                      conditions = paste(sort(unique(Condition)), collapse = ";")),
                  by = .(Mixture, Run, Fraction, run_short)]
        rs <- merge(rs, psms[, .(n_psms = .N), by = .(Run = Spectrum.File)], by = "Run")
        setorder(rs, Mixture, Fraction)
        fr <- rs[, .(n_fractions = .N, n_psms = sum(n_psms)), by = Mixture]
        message("fractions per mixture: ", paste(sprintf("%s=%d", fr$Mixture, fr$n_fractions), collapse = ", "))
        write_tsv(rs, file.path(root, "results", "mouse_run_structure.tsv"))
        design <- unique(ann[, .(Mixture, mixture_short, Channel, Condition, Diet, Duration, BioReplicate)])
        design_tab <- dcast(design[, .N, by = .(mixture_short, Condition)], Condition ~ mixture_short, value.var = "N", fill = 0L)
        write_tsv(design_tab, file.path(root, "results", "mouse_design_by_mixture.tsv"))
        print(design_tab)
    }

    ## One assay per run. readQFeatures names each sample <run>_<quantCol>,
    ## so sample names stay unique when the runs are joined later.
    qf <- readQFeatures(as.data.frame(psms), colData = as.data.frame(ann),
                        quantCols = unique(ann$quantCols), runCol = "Spectrum.File",
                        name = "psms", verbose = FALSE)
    short <- ann[match(names(qf), Run), run_short]
    names(qf) <- short
    rm(psms); invisible(gc())

    counts <- rbindlist(lapply(names(qf), function(i) {
        rd <- rowData(qf[[i]])
        data.table(run = i, n_psms = nrow(rd),
                   n_ions = length(unique(paste(rd$Annotated.Sequence, rd$Modifications, rd$Charge))),
                   n_proteins = length(unique(rd$Protein.Accessions)),
                   n_channels = ncol(qf[[i]]))
    }))
    write_tsv(counts, file.path(root, "results", paste0(dataset, "_psm_counts_raw.tsv")))
    message(sprintf("%s: %d runs, %d PSMs, %d distinct protein accessions",
                    dataset, length(qf), sum(counts$n_psms),
                    length(unique(unlist(lapply(names(qf), function(i) rowData(qf[[i]])$Protein.Accessions))))))
    saveRDS(qf, file.path(out_dir, paste0(dataset, "_raw.rds")))
    rm(qf); invisible(gc())
    counts
}

all_counts <- rbindlist(lapply(which_ds, function(d) { x <- build_one(d); x[, dataset := d]; x }))
write_tsv(all_counts, file.path(root, "results", "psm_counts_raw_all.tsv"))

stage_end(timer, note = paste("datasets:", paste(which_ds, collapse = ",")))
stage_mark_done(STAGE, conf)
