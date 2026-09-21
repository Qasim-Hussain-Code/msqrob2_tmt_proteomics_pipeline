# data/

Nothing in this directory is tracked by git except this file. The
tables are fetched by `scripts/02_download_data.sh`, verified against
the md5 sums that Zenodo publishes for record 14767905 (and against the
size and, after first retrieval, the md5 recorded in
`config/external_checksums.tsv` for the MassIVE files), cached in
`.bfc/` by BiocFileCache, and hard-linked here under their record names.

| file | source | size |
|---|---|---|
| `spikein1_psms.txt` | Zenodo 14767905, from MassIVE reanalysis RMSV000000265 of MSV000084264 (SPS-MS3) | 238.8 MB |
| `spikein1_annotations.csv` | Zenodo 14767905, from MassIVE MSV000084264 metadata | 13 kB |
| `spikein1_ms2_psms.txt` | MassIVE MSV000084266 (MS2-only acquisition of the same mixtures) | 301.4 MB |
| `spikein1_ms2_annotations.csv` | MassIVE MSV000084266 metadata | 11 kB |
| `mouse_psms.txt` | Zenodo 14767905, from MassIVE reanalysis RMSV000000264 of MSV000082569 | 80.9 MB |
| `mouse_annotations.csv` | Zenodo 14767905 | 47 kB |
| `spikein2_*.csv` (optional, `--dataset spikein2`) | Zenodo 14767905 | 24 MB |

`manifest.tsv` and `zenodo_record.tsv` are written by the download
stage and copied into `results/` as `data_checksums.tsv` and
`zenodo_record.tsv`, which are tracked.

`search_arm/` holds RAW files, mzML conversions and the search database
of the optional Sage arm (`scripts/10_search_arm.sh`) when it is run.
