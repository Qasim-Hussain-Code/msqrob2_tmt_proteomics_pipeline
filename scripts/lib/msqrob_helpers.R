## Model fitting, testing and diagnostics shared by 05 and 07.
##
## Everything here wraps msqrob2 so that a workflow variant is a row in a
## table (level, reference treatment, random effect structure, robust,
## ridge) and not a separate script. The functions return plain
## data.tables of per-protein, per-contrast results with the fit type
## attached, because the benchmark needs to count fit errors and
## singular fits as carefully as it counts significant proteins.

suppressPackageStartupMessages({
    library("data.table")
    library("QFeatures")
    library("msqrob2")
    library("lme4")
    library("limma")
})

## ---- Contrast handling ---------------------------------------------------

## Hypotheses are written once with plain parameter names. Ridge
## regression in msqrob2 prefixes every fixed effect with "ridge", so
## the same hypotheses need rewriting when ridge is on. Longest names
## first, whole tokens only, so that DurationLong inside
## DietHF:DurationLong is left alone.
ridge_names <- function(hypotheses, params) {
    params <- params[order(-nchar(params))]
    for (p in params) {
        pat <- paste0("(?<![A-Za-z0-9_.:])", gsub("([.:()])", "\\\\\\1", p), "(?![A-Za-z0-9_.:])")
        hypotheses <- gsub(pat, paste0("ridge", p), hypotheses, perl = TRUE)
    }
    hypotheses
}

build_contrasts <- function(hypotheses, params, ridge = FALSE) {
    if (ridge) {
        hypotheses <- ridge_names(hypotheses, params)
        params <- paste0("ridge", params)
    }
    makeContrast(hypotheses, parameterNames = params)
}

## ---- Result extraction -----------------------------------------------------

## msqrob2 stores a StatModel per protein; its @type is "lm", "rlm",
## "lmer" or "fitError". The hypothesis test writes one DataFrame per
## contrast into the rowData; both are gathered into one long table.
collect_results <- function(se, L, model_col = "msqrobModels", contrast_labels = NULL) {
    rd <- rowData(se)
    models <- rd[[model_col]]
    types <- vapply(models, function(m) m@type, character(1))
    dfs <- vapply(models, function(m) if (m@type == "fitError") NA_real_ else getDF(m), numeric(1))
    df_post <- vapply(models, function(m) if (m@type == "fitError") NA_real_ else m@dfPosterior, numeric(1))
    out <- rbindlist(lapply(seq_len(ncol(L)), function(k) {
        cn <- colnames(L)[k]
        res <- rd[[cn]]
        data.table(protein = rownames(rd),
                   contrast = if (is.null(contrast_labels)) cn else contrast_labels[k],
                   logFC = res$logFC, se = res$se, df = res$df, t = res$t,
                   pval = res$pval, adjPval = res$adjPval,
                   fit_type = types, df_residual = dfs, df_posterior = df_post)
    }))
    out
}

## ---- Fitting -------------------------------------------------------------

## Protein-level model on a summarised assay. Returns results and the
## wall clock of the fit alone.
fit_protein_level <- function(qf, i, formula, hypotheses, params, robust = TRUE, ridge = FALSE,
                              contrast_labels = NULL) {
    L <- build_contrasts(hypotheses, params, ridge)
    t0 <- Sys.time()
    qf <- msqrob(qf, i = i, formula = formula, robust = robust, ridge = ridge,
                 modelColumnName = "msqrobModels", overwrite = TRUE)
    fit_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    qf <- hypothesisTest(qf, i = i, contrast = L, modelColumn = "msqrobModels", overwrite = TRUE)
    res <- collect_results(qf[[i]], L, "msqrobModels", contrast_labels)
    list(results = res, fit_seconds = fit_s, n_proteins = nrow(qf[[i]]))
}

## PSM-level model: msqrobAggregate fits one mixed model per protein on
## all its ions at once and creates a summarised assay as a by-product.
fit_psm_level <- function(qf, i, formula, fcol, hypotheses, params, robust = TRUE, ridge = FALSE,
                          contrast_labels = NULL, name = "proteins_msqrob") {
    L <- build_contrasts(hypotheses, params, ridge)
    t0 <- Sys.time()
    qf <- msqrobAggregate(qf, i = i, fcol = fcol, formula = formula, robust = robust, ridge = ridge,
                          name = name, modelColumnName = "msqrobModels",
                          aggregateFun = MsCoreUtils::robustSummary)
    fit_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    qf <- hypothesisTest(qf, i = name, contrast = L, modelColumn = "msqrobModels", overwrite = TRUE)
    res <- collect_results(qf[[name]], L, "msqrobModels", contrast_labels)
    list(results = res, fit_seconds = fit_s, n_proteins = nrow(qf[[name]]))
}

## ---- Diagnostics: singular fits and variance components ------------------
##
## msqrob2 keeps only the coefficients and their covariance, not the
## lme4 object, so whether a fit was singular has to come from a second
## pass with lme4 itself. This pass fits the same formula (without
## robust weights or ridge) protein by protein and records the variance
## component of every random effect, the residual variance, and lme4's
## singularity and convergence messages. A singular fit is one where at
## least one variance component was estimated at exactly zero: the
## data for that protein carry no evidence of, say, a mixture effect
## beyond what the residuals already explain. With three mixtures that
## is expected to be common and is reported, not hidden.
lmer_diagnostics <- function(qf, i, formula, psm_level = FALSE, fcol = "Protein.Accessions",
                             max_proteins = Inf, seed = 20260921) {
    se <- getWithColData(qf, i)
    cd <- as.data.frame(colData(se))
    vars <- all.vars(formula)
    y <- assay(se)
    groups <- if (psm_level) rowData(se)[[fcol]] else rownames(se)
    prots <- unique(groups)
    if (is.finite(max_proteins) && length(prots) > max_proteins) {
        set.seed(seed); prots <- sample(prots, max_proteins)
    }
    rowvars <- intersect(vars, colnames(rowData(se)))
    colvars <- intersect(vars, colnames(cd))
    fml <- update.formula(formula, y ~ .)
    re_names <- vapply(findbars(formula), function(b) deparse(b[[3]]), character(1))
    out <- vector("list", length(prots))
    for (k in seq_along(prots)) {
        p <- prots[k]
        rows <- which(groups == p)
        ymat <- y[rows, , drop = FALSE]
        d <- data.frame(y = as.vector(ymat),
                        cd[rep(seq_len(ncol(ymat)), each = nrow(ymat)), colvars, drop = FALSE],
                        stringsAsFactors = FALSE)
        for (rv in rowvars) d[[rv]] <- rep(rowData(se)[rows, rv], times = ncol(ymat))
        d <- d[!is.na(d$y), , drop = FALSE]
        msgs <- character()
        fit <- tryCatch(withCallingHandlers(
            lmer(fml, data = d, control = lmerControl(calc.derivs = FALSE)),
            warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
            message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }),
            error = function(e) { msgs <<- c(msgs, conditionMessage(e)); NULL })
        if (is.null(fit)) {
            out[[k]] <- data.table(protein = p, n_obs = nrow(d), status = "error", singular = NA,
                                   sigma2 = NA_real_, messages = paste(unique(msgs), collapse = " | "))
            next
        }
        vc <- as.data.frame(VarCorr(fit))
        vc <- vc[is.na(vc$var2), ]
        comps <- setNames(vc$vcov, ifelse(vc$grp == "Residual", "sigma2", paste0("var_", vc$grp)))
        row <- data.table(protein = p, n_obs = nrow(d),
                          status = if (any(grepl("converge", msgs))) "convergence_warning" else "ok",
                          singular = isSingular(fit),
                          messages = paste(unique(msgs), collapse = " | "))
        for (nm in names(comps)) row[[nm]] <- comps[[nm]]
        out[[k]] <- row
    }
    res <- rbindlist(out, fill = TRUE)
    setcolorder(res, c("protein", "n_obs", "status", "singular"))
    res
}

## ---- limma reference arm -------------------------------------------------
##
## limma with duplicateCorrelation is the transcriptomics reader's
## calibration point: a fixed-effects design plus one consensus
## intra-block correlation shared by every protein, with the block set
## to the mixture. It cannot carry a second level (run within mixture),
## which is the point of including it.
fit_limma_dupcor <- function(mat, cd, design_formula, block, hypotheses, params, contrast_labels = NULL) {
    design <- model.matrix(design_formula, data = cd)
    colnames(design) <- sub("^Condition", "Condition", colnames(design))
    ## Rows for every design column, zero for parameters no hypothesis
    ## uses (the reference channel level when it is kept).
    L0 <- makeContrast(hypotheses, parameterNames = params)
    L <- matrix(0, nrow = ncol(design), ncol = ncol(L0), dimnames = list(colnames(design), colnames(L0)))
    L[rownames(L0), ] <- L0
    t0 <- Sys.time()
    dc <- duplicateCorrelation(mat, design, block = cd[[block]])
    fit <- lmFit(mat, design, block = cd[[block]], correlation = dc$consensus)
    fit2 <- eBayes(contrasts.fit(fit, L))
    fit_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    res <- rbindlist(lapply(seq_len(ncol(L)), function(k) {
        tt <- topTable(fit2, coef = k, number = Inf, sort.by = "none", adjust.method = "BH")
        data.table(protein = rownames(tt),
                   contrast = if (is.null(contrast_labels)) colnames(L)[k] else contrast_labels[k],
                   logFC = tt$logFC, se = tt$logFC / tt$t, df = fit2$df.total[1], t = tt$t,
                   pval = tt$P.Value, adjPval = tt$adj.P.Val, fit_type = "limma_dupcor",
                   df_residual = fit2$df.residual[1], df_posterior = fit2$df.total[1])
    }))
    list(results = res, fit_seconds = fit_s, n_proteins = nrow(mat), consensus_correlation = dc$consensus)
}

## ---- Benchmark metrics -----------------------------------------------------

## Realised false discovery proportion, sensitivity and fold change
## bias per contrast, against the design's known truth. A UPS1 protein
## is differential in every contrast; a HeLa protein never is.
benchmark_metrics <- function(res, expected, alpha = 0.05) {
    res <- copy(res)
    res[, is_ups := grepl("ups", protein)]
    res <- merge(res, expected[, .(contrast, expected_log2fc)], by = "contrast", all.x = TRUE)
    res[, .(
        n_tested = sum(!is.na(pval)),
        n_ups_tested = sum(!is.na(pval) & is_ups),
        n_hela_tested = sum(!is.na(pval) & !is_ups),
        n_fit_error = sum(fit_type == "fitError"),
        n_significant = sum(adjPval < alpha, na.rm = TRUE),
        tp = sum(adjPval < alpha & is_ups, na.rm = TRUE),
        fp = sum(adjPval < alpha & !is_ups, na.rm = TRUE),
        fdp = {
            s <- sum(adjPval < alpha, na.rm = TRUE)
            if (s == 0) NA_real_ else sum(adjPval < alpha & !is_ups, na.rm = TRUE) / s
        },
        sensitivity = sum(adjPval < alpha & is_ups, na.rm = TRUE) / sum(!is.na(pval) & is_ups),
        ups_median_log2fc = median(logFC[is_ups], na.rm = TRUE),
        ups_iqr_log2fc = IQR(logFC[is_ups], na.rm = TRUE),
        hela_median_log2fc = median(logFC[!is_ups], na.rm = TRUE),
        expected_log2fc = expected_log2fc[1],
        bias = median(logFC[is_ups], na.rm = TRUE) - expected_log2fc[1],
        alpha = alpha
    ), by = contrast]
}

## ---- Tiered fitting: refit proteins that the full model cannot fit ---------
##
## msqrob2 returns a fitError when the fixed effects are not estimable
## from a protein's data (an empty factorial cell) or when lme4 refuses
## the random effects (a grouping factor with as many levels as
## observations, or a single sampled level). The second kind is a
## property of the model, not of the protein: a protein quantified in
## one fraction per mouse gives one value per mouse, and a mouse random
## effect then has nothing to estimate beyond the residual. Following
## the msqrob2TMT vignette's refit of one-hit wonders, such proteins are
## refitted with the next simpler random structure in `tiers`, the
## moderated variance is re-estimated across all proteins together (as
## the vignette's msqrobRefit does), and the tier that produced each
## protein's result is recorded. Nothing is dropped silently: a protein
## that fails every tier stays a fitError in the results.
fit_tiered <- function(qf, i, tiers, hypotheses, params, robust = TRUE, ridge = FALSE,
                       psm_level = FALSE, fcol = "Protein.Accessions", name = "proteins_msqrob",
                       contrast_labels = NULL) {
    stopifnot(is.list(tiers), !is.null(names(tiers)))
    L <- build_contrasts(hypotheses, params, ridge)
    t0 <- Sys.time()
    fit_one <- function(obj, formula, set_name) {
        if (psm_level) {
            obj <- msqrobAggregate(obj, i = i, fcol = fcol, formula = formula, robust = robust, ridge = ridge,
                                   name = set_name, modelColumnName = "msqrobModels",
                                   aggregateFun = MsCoreUtils::robustSummary)
        } else {
            obj <- msqrob(obj, i = i, formula = formula, robust = robust, ridge = ridge,
                          modelColumnName = "msqrobModels", overwrite = TRUE)
        }
        obj
    }
    set <- if (psm_level) name else i
    ## msqrob2 builds the fixed-effect model matrix before its own error
    ## guard, so a protein whose observed channels cover a single level
    ## of a factor (all high fat, say) stops the whole call with R's
    ## "contrasts can be applied only to factors with 2 or more levels".
    ## Such proteins cannot support the fixed effects under any random
    ## structure; they are screened out here and reported as fitError.
    se0 <- getWithColData(qf, i)
    fixed_vars <- intersect(all.vars(nobars(tiers[[1]])), colnames(colData(se0)))
    groups0 <- if (psm_level) rowData(se0)[[fcol]] else rownames(se0)
    obs <- !is.na(assay(se0))
    cd0 <- as.data.frame(colData(se0))
    level_ok <- vapply(split(seq_len(nrow(obs)), groups0), function(rows) {
        cols <- colSums(obs[rows, , drop = FALSE]) > 0
        all(vapply(fixed_vars, function(v) length(unique(cd0[[v]][cols])) >= 2, logical(1)))
    }, logical(1))
    absent <- names(level_ok)[!level_ok]
    if (length(absent)) {
        keep_rows <- !groups0 %in% absent
        qf <- qf[keep_rows, , i]
    }
    qf <- fit_one(qf, tiers[[1]], set)
    models <- as.list(rowData(qf[[set]])[["msqrobModels"]])
    names(models) <- rownames(qf[[set]])
    is_err <- function(ms) vapply(ms, function(m) m@type == "fitError", logical(1))
    tier <- setNames(rep(names(tiers)[1], length(models)), names(models))
    tier[is_err(models)] <- "fitError"
    tier_log <- data.table(tier = names(tiers)[1], attempted = length(models) + length(absent), fitted = sum(!is_err(models)))
    ## msqrob2 moderates the variance inside every call with
    ## limma::squeezeVar, which errors when a subset's variances are
    ## degenerate (a tier where almost nothing fits). A tier is therefore
    ## tried whole and, on failure, in chunks of 25 proteins, so that only
    ## the chunk that cannot be moderated stays unfitted. The moderation
    ## is redone across all proteins afterwards in any case.
    fit_subset <- function(todo, formula) {
        sub <- if (psm_level) qf[rowData(qf[[i]])[[fcol]] %in% todo, , i] else qf[todo, , i]
        sub <- fit_one(sub, formula, "refit")
        newm <- as.list(rowData(sub[[if (psm_level) "refit" else i]])[["msqrobModels"]])
        names(newm) <- rownames(sub[[if (psm_level) "refit" else i]])
        newm
    }
    for (k in seq_along(tiers)[-1]) {
        todo <- names(tier)[tier == "fitError"]
        if (!length(todo)) break
        newm <- tryCatch(fit_subset(todo, tiers[[k]]), error = function(e) NULL)
        n_failed_chunks <- 0L
        if (is.null(newm)) {
            newm <- list()
            for (chunk in split(todo, ceiling(seq_along(todo) / 25))) {
                got <- tryCatch(fit_subset(chunk, tiers[[k]]), error = function(e) { n_failed_chunks <<- n_failed_chunks + 1L; NULL })
                if (!is.null(got)) newm <- c(newm, got)
            }
        }
        ok <- if (length(newm)) !is_err(newm) else logical()
        models[names(newm)[ok]] <- newm[ok]
        tier[names(newm)[ok]] <- names(tiers)[k]
        tier_log <- rbind(tier_log, data.table(tier = names(tiers)[k], attempted = length(todo), fitted = sum(ok)))
        if (n_failed_chunks) message("   tier ", names(tiers)[k], ": ", n_failed_chunks, " chunk(s) of 25 could not be variance-moderated and stay unfitted")
    }
    ## Re-estimate the moderated variance across every protein together;
    ## the tiers were squeezed separately by their own msqrob calls.
    vars <- vapply(models, function(m) if (m@type == "fitError") NA_real_ else getVar(m), numeric(1))
    dfs <- vapply(models, function(m) if (m@type == "fitError") NA_real_ else getDF(m), numeric(1))
    okv <- is.finite(vars) & vars > 0 & is.finite(dfs) & dfs > 0
    hlp <- limma::squeezeVar(var = vars[okv], df = dfs[okv])
    idx <- which(okv)
    for (j in seq_along(idx)) {
        models[[idx[j]]]@varPosterior <- as.numeric(hlp$var.post[j])
        models[[idx[j]]]@dfPosterior <- as.numeric(hlp$df.prior + dfs[idx[j]])
    }
    rowData(qf[[set]])[["msqrobModels"]] <- models[rownames(qf[[set]])]
    fit_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    qf <- hypothesisTest(qf, i = set, contrast = L, modelColumn = "msqrobModels", overwrite = TRUE)
    res <- collect_results(qf[[set]], L, "msqrobModels", contrast_labels)
    res[, fit_tier := tier[protein]]
    if (length(absent)) {
        res <- rbind(res, CJ(protein = absent, contrast = res[, unique(contrast)])[
            , `:=`(logFC = NA_real_, se = NA_real_, df = NA_real_, t = NA_real_, pval = NA_real_, adjPval = NA_real_,
                   fit_type = "fitError", df_residual = NA_real_, df_posterior = NA_real_, fit_tier = "fitError")], fill = TRUE)
        tier_log <- rbind(tier_log, data.table(tier = "absent_factor_level", attempted = length(absent), fitted = 0L))
    }
    tier_log[, formula := vapply(names(tiers), function(n) paste(deparse(tiers[[n]]), collapse = ""), character(1))[tier]]
    list(results = res, fit_seconds = fit_s, n_proteins = nrow(qf[[set]]) + length(absent), tier_log = tier_log,
         df_prior = hlp$df.prior, var_prior = hlp$var.prior, n_absent_level = length(absent))
}
