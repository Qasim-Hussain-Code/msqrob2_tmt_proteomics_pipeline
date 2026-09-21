## Install the R packages the pipeline needs and record what was
## actually installed. Run by scripts/01_install.sh; safe to re-run.
##
## Versions are deliberately not pinned in this file. The brief asks
## for whatever the current Bioconductor release provides at install
## time, with the result recorded (config/r_packages.tsv), so that a
## re-run six months from now reports its own versions rather than
## silently claiming ours.

args <- commandArgs(trailingOnly = TRUE)
out_tsv <- if (length(args) >= 1) args[1] else "config/r_packages.tsv"
out_session <- if (length(args) >= 2) args[2] else "logs/sessionInfo_install.txt"

pkgs <- c(
    "QFeatures", "msqrob2", "MsCoreUtils", "lme4", "limma",
    "data.table", "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork",
    "ComplexHeatmap", "ExploreModelMatrix", "BiocFileCache", "BiocParallel",
    "matrixStats", "jsonlite", "quarto", "rmarkdown", "knitr", "statmod", "R.utils"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

## Let BiocManager pick the Bioconductor release matching this R.
bioc_version <- tryCatch(as.character(BiocManager::version()), error = function(e) NA)
cat("R:", R.version.string, "\n")
cat("Bioconductor:", bioc_version, "\n")

installed <- rownames(installed.packages())
missing <- setdiff(pkgs, installed)
if (length(missing)) {
    cat("Installing:", paste(missing, collapse = ", "), "\n")
    ## update = FALSE keeps a working stack working; ask = FALSE keeps
    ## the script non-interactive. Ncpus speeds up source installs on
    ## Linux and is ignored for binaries.
    BiocManager::install(missing, update = FALSE, ask = FALSE,
                         Ncpus = max(1L, min(4L, parallel::detectCores())))
} else {
    cat("All packages already installed; nothing to do.\n")
}

still_missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(still_missing)) {
    stop("Packages failed to install: ", paste(still_missing, collapse = ", "))
}

ip <- installed.packages()[, c("Package", "Version")]
ip <- ip[ip[, "Package"] %in% pkgs, , drop = FALSE]
ip <- ip[order(ip[, "Package"]), , drop = FALSE]
tab <- data.frame(
    package = ip[, "Package"],
    version = ip[, "Version"],
    r_version = R.version.string,
    bioconductor = bioc_version,
    recorded = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
)
dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)
write.table(tab, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", out_tsv, "\n")

## Load everything once so a broken installation fails here and not
## three stages later inside a two hour model fit.
suppressPackageStartupMessages(for (p in pkgs) library(p, character.only = TRUE))
dir.create(dirname(out_session), showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), out_session)
cat("wrote", out_session, "\n")
