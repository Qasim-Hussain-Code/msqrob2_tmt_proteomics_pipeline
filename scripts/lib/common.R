## Shared helpers for every R stage: configuration, stage bookkeeping,
## resource logging and a few small utilities. Sourced at the top of
## scripts/03 to 08 and by the report.

suppressPackageStartupMessages({
    library("data.table")
})

## ---- Paths and configuration ---------------------------------------

## Resolve the repository root from the script location so that every
## stage can be run from any working directory.
project_root <- function() {
    env_root <- Sys.getenv("PROJECT_ROOT", unset = "")
    if (nzchar(env_root)) return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
    if (length(file_arg)) {
        return(normalizePath(file.path(dirname(file_arg), ".."), winslash = "/", mustWork = TRUE))
    }
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

## Read project.conf written by scripts/00_configure.sh. The file is a
## flat KEY=value list, deliberately shell-compatible so that bash and R
## stages read the same values.
read_conf <- function(root = project_root()) {
    conf_file <- file.path(root, "project.conf")
    if (!file.exists(conf_file)) stop("project.conf not found; run scripts/00_configure.sh first")
    lines <- readLines(conf_file, warn = FALSE)
    lines <- lines[!grepl("^\\s*#", lines) & grepl("=", lines)]
    keys <- sub("=.*$", "", lines)
    vals <- sub("^[^=]*=", "", lines)
    vals <- gsub('^"|"$', "", vals)
    conf <- as.list(vals)
    names(conf) <- keys
    conf$PROJECT_ROOT <- root
    for (k in c("THREADS", "RAM_GB", "PLAN_RAM_GB", "DISK_GB", "WORKERS")) {
        if (!is.null(conf[[k]])) conf[[k]] <- as.integer(conf[[k]])
    }
    conf
}

## Command line flags of the form --key value or --flag.
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
    out <- list()
    i <- 1L
    while (i <= length(args)) {
        a <- args[i]
        if (startsWith(a, "--")) {
            key <- sub("^--", "", a)
            if (i < length(args) && !startsWith(args[i + 1L], "--")) {
                out[[key]] <- args[i + 1L]
                i <- i + 2L
            } else {
                out[[key]] <- TRUE
                i <- i + 1L
            }
        } else {
            i <- i + 1L
        }
    }
    out
}

## ---- BiocParallel registration ---------------------------------------

## Serial by default. Each forked or socket worker receives its own copy
## of the QFeatures object, and with a 4 GB ceiling two copies of a
## quarter-gigabyte PSM table plus lme4's working memory run out of RAM
## before they win any wall clock. MulticoreParam forks and is
## unavailable on Windows, where SnowParam is the equivalent; both are
## only used when --workers > 1 is passed explicitly.
register_parallel <- function(workers = 1L) {
    suppressPackageStartupMessages(library("BiocParallel"))
    workers <- as.integer(workers)
    if (is.na(workers) || workers <= 1L) {
        register(SerialParam())
        message("BiocParallel: SerialParam (memory-bound default)")
    } else if (.Platform$OS.type == "windows") {
        register(SnowParam(workers = workers, progressbar = FALSE))
        message("BiocParallel: SnowParam with ", workers, " workers (Windows has no fork)")
    } else {
        register(MulticoreParam(workers = workers, progressbar = FALSE))
        message("BiocParallel: MulticoreParam with ", workers, " workers")
    }
    invisible(bpparam())
}

## ---- Resource logging --------------------------------------------------

## Peak resident set size of this R process in MB. Linux exposes it in
## /proc; Windows records it in the process object; macOS has no cheap
## peak figure so the current RSS is returned and flagged.
peak_rss_mb <- function() {
    pid <- Sys.getpid()
    if (file.exists("/proc/self/status")) {
        st <- readLines("/proc/self/status", warn = FALSE)
        hw <- grep("^VmHWM:", st, value = TRUE)
        if (length(hw)) return(as.numeric(gsub("[^0-9]", "", hw)) / 1024)
    }
    if (.Platform$OS.type == "windows") {
        cmd <- sprintf("(Get-Process -Id %d).PeakWorkingSet64", pid)
        out <- tryCatch(system2("powershell", c("-NoProfile", "-Command", shQuote(cmd)),
                                stdout = TRUE, stderr = FALSE), error = function(e) NA)
        val <- suppressWarnings(as.numeric(gsub("[^0-9]", "", out[length(out)])))
        if (length(val) && !is.na(val)) return(val / 1024^2)
    }
    out <- tryCatch(system2("ps", c("-o", "rss=", "-p", pid), stdout = TRUE), error = function(e) NA)
    val <- suppressWarnings(as.numeric(trimws(out[1])))
    if (!is.na(val)) return(val / 1024)
    NA_real_
}

## Stage timers write one TSV row per stage into logs/<stage>.resources.tsv.
## The README's usage table is generated from these files, so a claim
## that a stage runs in N minutes at M MB is always a measurement.
stage_begin <- function(stage, conf) {
    dir.create(file.path(conf$PROJECT_ROOT, "logs"), showWarnings = FALSE, recursive = TRUE)
    list(stage = stage, conf = conf, t0 = Sys.time(), p0 = proc.time())
}

stage_end <- function(timer, note = "") {
    elapsed <- as.numeric(difftime(Sys.time(), timer$t0, units = "secs"))
    rss <- peak_rss_mb()
    row <- data.frame(stage = timer$stage,
                      elapsed_s = round(elapsed, 1),
                      peak_rss_mb = round(rss, 0),
                      note = note,
                      finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                      stringsAsFactors = FALSE)
    out <- file.path(timer$conf$PROJECT_ROOT, "logs", paste0(timer$stage, ".resources.tsv"))
    write.table(row, out, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("%s finished: %.1f s, peak RSS %.0f MB", timer$stage, elapsed, rss))
    invisible(row)
}

## Idempotency: a stage that finds its stamp file skips, and says so.
stage_done_file <- function(stage, conf) file.path(conf$PROJECT_ROOT, "logs", paste0(stage, ".done"))

stage_already_done <- function(stage, conf, force = FALSE) {
    f <- stage_done_file(stage, conf)
    if (file.exists(f) && !isTRUE(force)) {
        message(stage, ": already completed on ", readLines(f, warn = FALSE)[1], "; skipping (pass --force to redo)")
        return(TRUE)
    }
    FALSE
}

stage_mark_done <- function(stage, conf) {
    writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), stage_done_file(stage, conf))
}

## ---- Data manifest ---------------------------------------------------

## The download stage writes data/manifest.tsv; later stages look files
## up by their logical name rather than by path.
read_manifest <- function(conf) {
    f <- file.path(conf$PROJECT_ROOT, "data", "manifest.tsv")
    if (!file.exists(f)) stop("data/manifest.tsv not found; run scripts/02_download_data.sh first")
    m <- fread(f, sep = "\t")
    m[status == "verified"]
}

manifest_path <- function(conf, name) {
    m <- read_manifest(conf)
    p <- m[file == name, path]
    if (!length(p)) stop("file ", name, " is not in the verified manifest")
    p[1]
}

## ---- Small utilities ---------------------------------------------------

write_tsv <- function(x, path) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    ## Unquoted output, so that the tables read plainly; any tab or
    ## newline inside a text field (lme4 messages end with one) is
    ## replaced by a space first, or it would break the row.
    x <- as.data.table(x)
    for (col in names(x)) if (is.character(x[[col]])) set(x, j = col, value = gsub("[\t\r\n]+", " ", x[[col]]))
    fwrite(x, path, sep = "\t", na = "NA", quote = FALSE)
    message("wrote ", path)
}

save_fig <- function(plot, path, width, height, dpi = 150) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white")
    message("wrote ", path)
}

## Design expectations for the UPS1 spike-in, computed from the fmol
## amounts and not from the rounded dilution factors: 333 fmol is not
## exactly two thirds of 500, so log2(500/333) is 0.5865, not 0.5850.
spikein_conditions <- function() {
    data.table(condition = c("1", "0.667", "0.5", "0.125"),
               fmol = c(500, 333, 250, 62.5))
}

spikein_expected_contrasts <- function() {
    cond <- spikein_conditions()
    cmb <- combn(seq_len(nrow(cond)), 2)
    data.table(
        contrast = paste0(cond$condition[cmb[1, ]], " - ", cond$condition[cmb[2, ]]),
        high = cond$condition[cmb[1, ]],
        low = cond$condition[cmb[2, ]],
        fmol_high = cond$fmol[cmb[1, ]],
        fmol_low = cond$fmol[cmb[2, ]],
        ratio = cond$fmol[cmb[1, ]] / cond$fmol[cmb[2, ]],
        expected_log2fc = log2(cond$fmol[cmb[1, ]] / cond$fmol[cmb[2, ]])
    )
}

## Memory checkpoints. Each call appends heap and peak RSS to
## logs/<stage>.memory.tsv so that the step responsible for a stage's
## peak is identifiable after the fact.
mem_checkpoint <- function(timer, label) {
    g <- gc(full = TRUE)
    heap <- sum(g[, 2])
    rss <- peak_rss_mb()
    f <- file.path(timer$conf$PROJECT_ROOT, "logs", paste0(timer$stage, ".memory.tsv"))
    row <- data.frame(stage = timer$stage, checkpoint = label, heap_mb = round(heap),
                      peak_rss_mb = round(rss), elapsed_s = round(as.numeric(difftime(Sys.time(), timer$t0, units = "secs")), 1))
    write.table(row, f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
    message(sprintf("    [mem] %-32s heap %5.0f MB  peak RSS %5.0f MB", label, heap, rss))
    invisible(row)
}
