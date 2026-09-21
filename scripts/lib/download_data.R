## Fetch and verify the input tables.
##
## Zenodo record 14767905 is enumerated through its REST API so that the
## filenames, sizes and md5 sums come from the record itself, never from
## this script. Every file is downloaded into a BiocFileCache under
## data/.bfc, verified against the record's md5, and hard-linked into
## data/ under its record filename. A checksum mismatch stops the run;
## a file that cannot be verified is not analysed.
##
## The MS2-only counterpart of the spike-in (MassIVE MSV000084266) is
## enumerated from the MassIVE dataset file listing, which publishes
## sizes but no checksums. Its size is verified against the listing and
## the md5 observed at first download is recorded so that later runs
## detect a changed file.
##
## Called by scripts/02_download_data.sh.

suppressPackageStartupMessages({
    library("BiocFileCache")
    library("jsonlite")
    library("curl")
    library("data.table")
})

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."), winslash = "/")
source(file.path(root, "scripts", "lib", "common.R"))
conf <- read_conf(root)
opt <- parse_args()

datasets <- strsplit(if (is.null(opt$datasets)) "spikein1,mouse" else opt$datasets, ",")[[1]]
want_ms2 <- !isTRUE(opt$`no-ms2`)
force <- isTRUE(opt$force)

ZENODO_RECORD <- "14767905"
ZENODO_API <- sprintf("https://zenodo.org/api/records/%s", ZENODO_RECORD)
MASSIVE_MS2 <- "MSV000084266"

data_dir <- file.path(root, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
bfc <- BiocFileCache(file.path(data_dir, ".bfc"), ask = FALSE)

## ---- helpers ------------------------------------------------------------

md5 <- function(path) unname(tools::md5sum(path))

fetch_json <- function(url) {
    h <- new_handle(followlocation = TRUE, useragent = "msqrob2_tmt_proteomics_pipeline/02_download")
    r <- curl_fetch_memory(url, handle = h)
    if (r$status_code != 200) stop("HTTP ", r$status_code, " for ", url)
    fromJSON(rawToChar(r$content), simplifyVector = FALSE)
}

## Download to a fresh cache slot. On any failure the slot is removed so
## the cache never holds a half-written or wrong file.
cached_download <- function(rname, url, expected_size = NA, expected_md5 = NA) {
    hit <- bfcquery(bfc, rname, field = "rname", exact = TRUE)
    if (nrow(hit) && !force) {
        path <- bfc[[hit$rid[1]]]
        if (file.exists(path)) {
            ok_size <- is.na(expected_size) || file.size(path) == expected_size
            obs <- md5(path)
            ok_md5 <- is.na(expected_md5) || identical(obs, expected_md5)
            if (ok_size && ok_md5) {
                message("cached and verified: ", rname)
                return(list(path = path, md5 = obs, size = file.size(path), fresh = FALSE))
            }
            message("cached copy of ", rname, " fails verification; re-downloading")
        }
        bfcremove(bfc, hit$rid)
    } else if (nrow(hit) && force) {
        bfcremove(bfc, hit$rid)
    }
    path <- bfcnew(bfc, rname = rname, ext = paste0(".", tools::file_ext(rname)))
    message("downloading ", rname, " (", if (is.na(expected_size)) "size unknown" else sprintf("%.1f MB", expected_size / 1e6), ")")
    h <- new_handle(followlocation = TRUE, useragent = "msqrob2_tmt_proteomics_pipeline/02_download",
                    connecttimeout = 60L, low_speed_limit = 1024L, low_speed_time = 120L)
    ok <- tryCatch({ curl_download(url, path, handle = h, quiet = TRUE); TRUE },
                   error = function(e) { message("download failed: ", conditionMessage(e)); FALSE })
    if (!ok) { bfcremove(bfc, names(path)); stop("download failed for ", rname) }
    size <- file.size(path)
    if (!is.na(expected_size) && size != expected_size) {
        bfcremove(bfc, names(path))
        stop(sprintf("size mismatch for %s: expected %d bytes, got %d", rname, expected_size, size))
    }
    obs <- md5(path)
    if (!is.na(expected_md5) && !identical(obs, expected_md5)) {
        bfcremove(bfc, names(path))
        stop(sprintf("md5 mismatch for %s: expected %s, got %s. Refusing to analyse an unverified file.",
                     rname, expected_md5, obs))
    }
    list(path = unname(path), md5 = obs, size = size, fresh = TRUE)
}

## Expose the file under data/<name>. A hard link costs no disk; a copy
## is the fallback on filesystems that refuse links.
link_into_data <- function(src, name) {
    dst <- file.path(data_dir, name)
    if (file.exists(dst)) {
        if (identical(md5(dst), md5(src))) return(dst)
        file.remove(dst)
    }
    ok <- tryCatch(file.link(src, dst), warning = function(w) FALSE, error = function(e) FALSE)
    if (!isTRUE(ok)) file.copy(src, dst, overwrite = TRUE)
    dst
}

manifest <- list()
add_row <- function(dataset, file, src, res, url, exp_size, exp_md5, note = "") {
    manifest[[length(manifest) + 1L]] <<- data.table(
        dataset = dataset, file = file, path = res$path, linked_path = link_into_data(res$path, file),
        source = src, url = url,
        size_expected = exp_size, size_observed = res$size,
        md5_expected = ifelse(is.na(exp_md5), "", exp_md5), md5_observed = res$md5,
        status = "verified", retrieved = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), note = note)
}

## ---- Zenodo -------------------------------------------------------------

message("enumerating Zenodo record ", ZENODO_RECORD)
rec <- fetch_json(ZENODO_API)
## "key" is a reserved argument of data.table(), hence "fname".
zfiles <- rbindlist(lapply(rec$files, function(f) data.table(
    fname = f$key, size = as.numeric(f$size),
    md5 = sub("^md5:", "", f$checksum), url = f$links$self)))
record_info <- data.table(
    record_id = rec$id, doi = rec$doi, concept_doi = rec$conceptdoi,
    revision = rec$revision, publication_date = rec$metadata$publication_date,
    updated = rec$updated, license = rec$metadata$license$id,
    title = rec$metadata$title, n_files = nrow(zfiles), total_bytes = sum(zfiles$size))
message(sprintf("record %s revision %s, %d files, %.1f MB", rec$id, rec$revision, nrow(zfiles), sum(zfiles$size) / 1e6))
print(zfiles[, .(fname, size, md5)])

## Files are selected by dataset prefix, not by name. One of the record's
## description lines spells a filename "spilein2"; the API listing is
## the authority and whatever it says is what gets fetched.
for (ds in datasets) {
    sel <- zfiles[startsWith(fname, paste0(ds, "_")) & !grepl("[.]zip$", fname)]
    if (!nrow(sel)) stop("no Zenodo files with prefix ", ds, "_ in record ", ZENODO_RECORD)
    for (i in seq_len(nrow(sel))) {
        res <- cached_download(sel$fname[i], sel$url[i], sel$size[i], sel$md5[i])
        add_row(ds, sel$fname[i], paste0("zenodo:", ZENODO_RECORD), res, sel$url[i], sel$size[i], sel$md5[i])
    }
}

## ---- MassIVE MS2 counterpart -------------------------------------------

if (want_ms2 && "spikein1" %in% datasets) {
    message("enumerating MassIVE ", MASSIVE_MS2, " file listing")
    h <- new_handle(followlocation = TRUE)
    r <- curl_fetch_memory(sprintf("https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=%s", MASSIVE_MS2), handle = h)
    task <- sub(".*task=([0-9a-f]{32}).*", "\\1", r$url)
    if (!grepl("^[0-9a-f]{32}$", task)) stop("could not resolve MassIVE task id for ", MASSIVE_MS2)
    page <- rawToChar(curl_fetch_memory(sprintf("https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=%s", task), handle = h)$content)
    js <- regmatches(page, regexpr("var dataset_files = \\{.*?\\};", page))
    js <- sub("^var dataset_files = ", "", sub(";$", "", js))
    listing <- fromJSON(js, simplifyVector = FALSE)$row_data
    mfiles <- rbindlist(lapply(listing, function(f) data.table(
        descriptor = f$file_descriptor, name = f$name, collection = f$collection, size = as.numeric(f$size))))
    psm <- mfiles[collection == "quant" & grepl("_PSMs\\.txt$", name)]
    ann <- mfiles[collection == "metadata" & grepl("annotation\\.csv$", name)]
    if (nrow(psm) != 1L || nrow(ann) != 1L) {
        stop("expected exactly one PSM export and one annotation file in ", MASSIVE_MS2,
             "; found ", nrow(psm), " and ", nrow(ann))
    }
    mk_url <- function(desc) sprintf("https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=%s&forceDownload=true",
                                     URLencode(desc, reserved = TRUE))
    ## md5 recorded at first retrieval; later runs verify against it.
    prior_file <- file.path(root, "config", "external_checksums.tsv")
    prior <- if (file.exists(prior_file)) fread(prior_file) else data.table(file = character(), md5 = character())
    for (spec in list(list(name = "spikein1_ms2_psms.txt", row = psm), list(name = "spikein1_ms2_annotations.csv", row = ann))) {
        exp_md5 <- prior[file == spec$name, md5]
        exp_md5 <- if (length(exp_md5)) exp_md5[1] else NA_character_
        res <- cached_download(spec$name, mk_url(spec$row$descriptor), spec$row$size, exp_md5)
        add_row("spikein1_ms2", spec$name, paste0("massive:", MASSIVE_MS2), res, mk_url(spec$row$descriptor),
                spec$row$size, exp_md5, note = paste0("MassIVE path ", spec$row$descriptor,
                                                      if (is.na(exp_md5)) "; md5 recorded at first retrieval, no published checksum" else ""))
    }
    ext <- rbindlist(lapply(manifest, function(m) m[dataset == "spikein1_ms2", .(file, md5 = md5_observed, size = size_observed, source, massive_path = sub("^MassIVE path ", "", sub(";.*$", "", note)))]))
    if (nrow(ext)) fwrite(ext, prior_file, sep = "\t")
}

## ---- Write manifests -----------------------------------------------------

man <- rbindlist(manifest)
fwrite(man, file.path(data_dir, "manifest.tsv"), sep = "\t")
fwrite(record_info, file.path(data_dir, "zenodo_record.tsv"), sep = "\t")
dir.create(file.path(root, "results"), showWarnings = FALSE)
pub <- man[, .(dataset, file, source, size_bytes = size_observed, md5 = md5_observed,
               md5_source = ifelse(nzchar(md5_expected), "record", "observed at first retrieval"),
               status, retrieved)]
fwrite(pub, file.path(root, "results", "data_checksums.tsv"), sep = "\t")
fwrite(record_info, file.path(root, "results", "zenodo_record.tsv"), sep = "\t")
message("verified ", nrow(man), " files; manifest written to data/manifest.tsv and results/data_checksums.tsv")
