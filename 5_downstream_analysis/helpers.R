## helpers.R — consolidated helper library for the downstream-analysis QMDs.
##
## Replaces the former split into `global_analysis_helpers.R` and
## `rt_analysis_function.R`. Only functions actually used by a QMD are
## kept; orphaned dead code from the old `rt_analysis_function.R`
## (compute_bin_ratio, compute_bin_count, SavitzkyGolayParam,
## IsolatePeakParam, flagBadEICs, BaselineParam, calculatePeakMetrics,
## plot_signal_comparison — all unused as of the consolidation) was
## dropped.
##
## Public surface (each function lists the QMD(s) that source it):
##
##   load_data()                      — QMD 1 (sp_full + sp_full_detect builds)
##   detect_signal()                  — QMD 1 (xcms-detected per-lab object build)
##   run_for_polarities()             — QMD 7 (polarity-iteration wrapper)
##   run_for_polarities_names()       — QMD 7 (variant with pol name)
##   check_metric_quality()           — QMD 7 §1.2.1 (variance decomposition via lm)
##   compute_icc()                    — QMDs 5 §P3.5/§P3.6 and 7 §1.2.4
##                                      (variance-component ICC on log2 values)
##   plot_violin()                    — QMD 7 §1.1.1 (per-(lab, polarity, data_level)
##                                      violin plot of bpc/tic intensities)

# ── 1. xcms / Spectra source loaders (QMD 1) ───────────────────────

#' Load the lab × study_group MS-Experiment from disk into a Spectra
#' (or full MsExperiment) object, with `mixture` and `lab`
#' spectraVariables populated and MS1-only filtering applied.
#' @param lab        character — "afekta" / "cembio" / "hmgu" / "icl"
#' @param study_group character — "HE" (currently the only study group)
#' @param return     "sp" (default; returns Spectra) or "mse" (returns
#'                   the MsExperiment for chromPeaks() etc.)
#' Used by: 5_downstream_analysis/1_full_detected_objects.qmd
load_data <- function(lab, study_group, return = "sp") {
  dr <- here::here("1_preprocessing", lab, study_group)
  mse <- readMsObject(
    XcmsExperiment(),
    AlabasterParam(path = file.path(dr, "mse")),
    spectraPath = file.path(dr, "mzml"))
  sampleData(mse)$mixture <- sub(".*_", "", sampleData(mse)$Sample.Name)
  sampleData(mse)$mixture <- gsub("\\.", "_", sampleData(mse)$mixture)
  spectra(mse)$lab <- lab
  spectra(mse)$mixture <- sampleData(mse)[
    match(spectra(mse)$dataOrigin, sampleData(mse)$spectraOrigin),
    "mixture"]
  sp <- spectra(mse)
  if (length(unique(sp$msLevel)) > 1) {
    sp <- filterMsLevel(sp, 1L)
  }
  if (return == "sp") {
    return(sp)
  } else if (length(unique(spectra(mse)$msLevel)) > 1) {
    mse <- filterMsLevel(mse, 1L)
  } else {
    return(mse)
  }
}

#' Build the detected-signal Spectra object for one lab by intersecting
#' the lab's raw spectra with its xcms chromPeaks. With `annotated =
#' TRUE`, restricts to the per-lab ring-trial library subset; with
#' `annotated = FALSE`, also writes `detected_peaks_<lab>_HE.csv` to
#' `5_downstream_analysis/object/`.
#' Used by: 5_downstream_analysis/1_full_detected_objects.qmd
detect_signal <- function(study_group, annotated = TRUE,
                           lab = character(), bpparam) {

  a <- load_data(lab = lab, study_group = study_group, return = "mse")

  if (annotated) {
    res <- read.csv(here::here("4_library_generation", lab, study_group,
                                "ring_trial_library_HE.csv"))
    cpks <- chromPeaks(a)[res$chrom_peak_id,
      c("rtmin", "rtmax", "mzmin", "mzmax", "sample")]
  } else {
    cpks <- chromPeaks(a)[, c("rtmin", "rtmax", "mzmin", "mzmax",
                              "rt", "mz", "sample", "into", "intb", "maxo")]
    cpks <- as.data.frame(cpks)
    cpks$polarity <- sampleData(a)$polarity[
      match(cpks$sample, seq_len(nrow(sampleData(a))))]
    cpks$mixture <- sampleData(a)$mixture[
      match(cpks$sample, seq_len(nrow(sampleData(a))))]
    cpks$chrom_peak_id <- row.names(cpks)
    write.csv(cpks,
      file = here::here("5_downstream_analysis", "object",
        paste0("detected_peaks_", lab, "_", study_group, ".csv")),
      row.names = FALSE)
  }

  spectra(a) <- setBackend(spectra(a), MsBackendMemory())
  cpk_split <- split(as.data.frame(cpks), cpks[, "sample"])
  cpk_split <- lapply(cpk_split, function(df)
    as.matrix(df[, c("rtmin", "rtmax", "mzmin", "mzmax")]))

  bg <- bpmapply(
    FUN = function(s, pks) {
      s <- Spectra::filterPeaksRanges(s,
        mz = pks[, c("mzmin", "mzmax")],
        rtime = pks[, c("rtmin", "rtmax")],
        keep = TRUE)
      Spectra::applyProcessing(s)
    },
    split(spectra(a), spectraSampleIndex(a)),
    cpk_split,
    BPPARAM = bpparam)

  bg_full <- concatenateSpectra(bg)
  bg_full <- filterEmptySpectra(bg_full)
  return(bg_full)
}

# ── 2. Polarity-iteration wrappers (QMD 7) ─────────────────────────

#' `bind_rows` over (positive, negative) calling
#' `fn(pol_int, pol_name, ...)` for each — collapses the boilerplate
#' two-polarity loop used across QMD 7's metric builders.
#' Used by: 5_downstream_analysis/7_reproducibility_results.qmd
run_for_polarities <- function(fn, ...) {
  dplyr::bind_rows(
    fn(1L, "Positive", ...),
    fn(0L, "Negative", ...))
}

#' Variant of `run_for_polarities()` that doesn't forward extra args.
#' Used by: 5_downstream_analysis/7_reproducibility_results.qmd
run_for_polarities_names <- function(fn) {
  dplyr::bind_rows(
    fn(1L, "Positive"),
    fn(0L, "Negative"))
}

# ── 3. Variance decomposition (lm + ANOVA) — QMD 7 §1.2 ───────────

#' Fit `lm(log2(metric) ~ lab + mixture)` and partition the ANOVA Sum
#' Sq into Lab / Mixture / Residual variance percentages. Plus
#' coefficients flagged as Lab vs Mixture effects.
#' Used by: 5_downstream_analysis/7_reproducibility_results.qmd §1.2.1
#' @return list(model, r_squared, variance_table, coefficients)
check_metric_quality <- function(data, metric_name) {
  if (metric_name %in% c("tic", "bpc")) {
    model <- lm(log2(get(metric_name)) ~ lab + mixture, data = data)
  } else {
    model <- lm(get(metric_name) ~ lab + mixture, data = data)
  }

  r_squared <- summary(model)$r.squared
  an <- anova(model)
  ss_total <- sum(an[, "Sum Sq"])

  variance_data <- data.frame(
    Source     = c("Lab Variance", "Mixture Variance", "Residual Variance"),
    Percentage = c(
      an["lab", "Sum Sq"]       / ss_total * 100,
      an["mixture", "Sum Sq"]   / ss_total * 100,
      an["Residuals", "Sum Sq"] / ss_total * 100))

  coefs <- coef(model)
  coef_df <- data.frame(
    coefficient = names(coefs),
    estimate    = coefs,
    type        = dplyr::case_when(
      grepl("^lab", names(coefs))     ~ "Lab Effect",
      grepl("^mixture", names(coefs)) ~ "Mixture Effect",
      TRUE                            ~ "Intercept")) |>
    dplyr::filter(type != "Intercept")

  list(model         = model,
       r_squared     = r_squared,
       variance_table = variance_data,
       coefficients  = coef_df)
}

# ── 4. Variance-component ICC on log2 values ──────────────────────

#' Variance-component ICC. log2-transforms `metric_name` (unless its
#' name contains "ratio") and computes
#'   ICC = Var(between-mixture-means)
#'         / [Var(between) + mean(Var(within-mixture))]
#' across the `(lab, mixture)` rows in `data`. Equivalent to ICC(1)
#' from a one-way random-effects model with mixture as the random
#' factor — computed directly from the variance components so the
#' formula is transparent (no `lme4`).
#' Used by:
#'   - 5_downstream_analysis/5_normalization.qmd §P3.5 / §P3.6
#'   - 5_downstream_analysis/7_reproducibility_results.qmd §1.2.4
#' @return data.frame(Metric, ICC, Between_mixture_var, Within_mixture_var)
compute_icc <- function(data, metric_name) {
  if (grepl("ratio", metric_name, ignore.case = TRUE)) {
    data$value <- data[[metric_name]]
  } else {
    data$value <- log2(data[[metric_name]])
  }

  between_var <- data |>
    dplyr::group_by(mixture) |>
    dplyr::summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::pull(mean_val) |>
    var(na.rm = TRUE)

  within_var <- data |>
    dplyr::group_by(mixture) |>
    dplyr::summarise(var_val = var(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::pull(var_val) |>
    mean(na.rm = TRUE)

  data.frame(
    Metric              = metric_name,
    ICC                 = between_var / (between_var + within_var),
    Between_mixture_var = between_var,
    Within_mixture_var  = within_var)
}

# ── 5. Plot helpers ───────────────────────────────────────────────

#' Per-(lab, polarity, mixture) violin plot of a bulk metric, faceted
#' by `data_level` (single row, left-to-right Full → Detected →
#' Consensus → Annotated).
#' Used by: 5_downstream_analysis/7_reproducibility_results.qmd §1.1.1
plot_violin <- function(data, metric) {
  if (metric %in% c("bpc", "tic")) {
    data <- data |> dplyr::mutate(dplyr::across(dplyr::all_of(metric), log2))
    ylab  <- paste(metric, "(log2)")
    title <- paste("Violin plot of", metric, "(log2 scale)")
  } else {
    ylab  <- metric
    title <- paste("Violin plot of", metric)
  }

  ggplot2::ggplot(data,
                  ggplot2::aes(x = lab, y = .data[[metric]], fill = polarity)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.6) +
    ggplot2::geom_jitter(ggplot2::aes(color = as.factor(mixture)),
                          width = 0.2, size = 1.5, alpha = 0.7) +
    ggplot2::scale_color_hue() +
    ggplot2::facet_wrap(~ data_level, nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = "Laboratory", y = ylab, color = "Mixture")
}
