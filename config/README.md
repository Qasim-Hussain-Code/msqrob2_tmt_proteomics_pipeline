# config/

- `environment.yml`: conda specification for R and the command line tools; Bioconductor packages are installed by `scripts/lib/install_packages.R` inside that R.
- `r_packages.tsv`: the package versions that were in use when the analysis ran, written by `scripts/01_install.sh`. Not hand-edited.
- `external_checksums.tsv`: md5 and size of the MassIVE files, recorded at first retrieval because MassIVE publishes no checksums.
- `contrasts_spikein1.tsv`: the six pairwise contrasts of the UPS1 dilution series, as msqrob2 hypotheses on the `~ 0 + Condition` parameterisation. Expected log2 fold changes are computed from the fmol columns in `scripts/lib/common.R`, not stored.
- `contrasts_mouse.tsv`: the four contrasts of the diet by duration factorial on the `~ Diet * Duration` parameterisation with LF and Short (8 weeks) as reference levels.
- `sage_spikein_ms2.json`: Sage search parameters for the optional search arm (`scripts/10_search_arm.sh`).
- Sample annotations are not copied here: they are downloaded and checksum-verified by `scripts/02_download_data.sh` and live in `data/`.
