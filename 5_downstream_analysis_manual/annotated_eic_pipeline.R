## annotated_eic_pipeline.R
##
## Stage 3 of the MS1-anchored annotation pipeline. Sourced by
## `3_rt_alignment.qmd` §P2.10 and by `ms1_anchor.R`.
##
## Inputs:
##   - feat_rts_ms1_anchored.csv — per-row rt_sec snapped to the Stage-2
##     consensus candidate (produced by ms1_anchor.R Stage 1 + 2)
##   - annotation_identity.csv  — SIRIUS identifications + cosmic scores
## Outputs:
##   - annotation_chrom_metrics.csv — per (lab, mixture, polarity,
##     compound, adduct) Stage-3 picker result
##
## Stage 2's cross-lab consensus + isobaric resolution have already
## run upstream; this file's old greedy-cluster + isobar-check pass
## remains as a defensive audit but is expected to be a no-op (each
## (mixture, compound) collapses to a single cluster, 0 collisions).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(here)
  library(Spectra); library(Chromatograms); library(MsExperiment); library(xcms)
})

# ───────────────────────────────────────────────────────────────────
# 1. Helpers
# ───────────────────────────────────────────────────────────────────

# Per (lab, polarity, mixture), the MS1 mzML path.
build_lab_files <- function(labs, mixtures, study_group = "HE") {
  mix_us <- function(m) gsub("\\.", "_", m)

  ms1_files_for_lab <- function(lab_name, study_group, mixture) {
    study_path <- here("1_preprocessing", lab_name, study_group)
    mzml_dir   <- file.path(study_path, "mzml")
    seq_pos    <- read_xlsx(file.path(study_path, "seq_pos.xlsx")) |>
      as.data.frame()
    seq_neg    <- read_xlsx(file.path(study_path, "seq_neg.xlsx")) |>
      as.data.frame()
    seq_pos$polarity <- "pos"; seq_neg$polarity <- "neg"
    seq <- rbind(seq_pos, seq_neg)
    seq$filename <- paste0(seq$`Data File`, ".mzML")
    seq <- seq[seq$filename %in% list.files(mzml_dir), ]
    seq$mixture <- gsub("\\.", "_",
                        sub("target$", "",
                            sub(".*_", "",
                                sub("_MS2$", "", seq$`Data File`))))
    # MS1 files only (exclude *_MS2*.mzML)
    seq <- seq[!grepl("_MS2", seq$filename), ]
    separate <- any(grepl("_MS2", list.files(mzml_dir), fixed = TRUE))
    seq <- if (separate)
      seq[grepl("MS2", seq$filename, fixed = TRUE) == FALSE &
            seq$mixture == mix_us(mixture), ]
    else
      seq[seq$mixture == mix_us(mixture), ]
    tibble(lab = lab_name, mixture = mixture, polarity = seq$polarity,
           path = file.path(mzml_dir, seq$filename))
  }

  bind_rows(lapply(mixtures, function(m)
    bind_rows(lapply(labs, ms1_files_for_lab,
                     study_group = study_group, mixture = m))))
}

# Greedy 1-D clustering: sort, walk, start new cluster on gap > window.
cluster_rts_greedy <- function(rt, window = 30) {
  if (length(rt) == 0L) return(integer(0))
  ord <- order(rt); s <- rt[ord]
  out <- integer(length(rt))
  cl <- 1L; out[ord[1L]] <- cl
  if (length(s) >= 2L) for (i in 2:length(s)) {
    if (s[i] - s[i - 1L] > window) cl <- cl + 1L
    out[ord[i]] <- cl
  }
  out
}

# Forward apply: rt_sec → rt_aligned (consensus space). Identity
# fallback for out-of-anchor-range values.
apply_fwd <- function(rt_sec, lab, polarity, splines, anchor_table) {
  out <- rt_sec
  key <- paste(lab, polarity, sep = "_")
  if (is.null(splines[[key]])) return(out)
  rng <- anchor_table |>
    dplyr::filter(lab == !!lab, polarity == !!polarity) |>
    dplyr::pull(rt_local) |> range(na.rm = TRUE)
  in_range <- !is.na(rt_sec) & rt_sec >= rng[1] & rt_sec <= rng[2]
  out[in_range] <- as.numeric(splines[[key]](rt_sec[in_range]))
  out
}

# Build inverse splines (consensus → lab-local) once.
build_inverse_splines <- function(splines, anchor_table) {
  inv <- list()
  for (k in names(splines)) {
    parts <- strsplit(k, "_")[[1]]
    sub <- anchor_table |>
      dplyr::filter(lab == parts[1], polarity == parts[2]) |>
      dplyr::arrange(consensus_rt)
    if (nrow(sub) < 4L) next
    yf <- isoreg(sub$consensus_rt, sub$rt_local)$yf
    df <- aggregate(yf ~ consensus_rt, data = data.frame(
      consensus_rt = sub$consensus_rt, yf = yf), FUN = median)
    df <- df[order(df$consensus_rt), ]
    inv[[k]] <- tryCatch(
      splinefun(df$consensus_rt, df$yf, method = "hyman"),
      error = function(e) approxfun(df$consensus_rt, df$yf, rule = 2))
  }
  inv
}

apply_inv <- function(rt_aligned, lab, polarity, inv_splines, anchor_table) {
  out <- rt_aligned
  key <- paste(lab, polarity, sep = "_")
  if (is.null(inv_splines[[key]])) return(out)
  rng <- anchor_table |>
    dplyr::filter(lab == !!lab, polarity == !!polarity) |>
    dplyr::pull(consensus_rt) |> range(na.rm = TRUE)
  in_range <- !is.na(rt_aligned) &
              rt_aligned >= rng[1] & rt_aligned <= rng[2]
  out[in_range] <- as.numeric(inv_splines[[key]](rt_aligned[in_range]))
  out
}

# SNR-aware peak picker. Linear interior NA-fill of the EIC so sparse
# scan grids (notably icl, and some late-RT hmgu m/z) don't
# under-sample a real peak to a single apex scan that the
# n_scans_in_shape guard then rejects. Imputation only rescues high-SNR
# peaks; it cannot fabricate a detection on noise because the SNR
# guard (snr >= 5) still applies.
extract_one <- function(rt, intensity, rt_prior,
                        apex_search_pad   = 10,
                        tail_pad          = 25,
                        snr_threshold     = 5,
                        abs_min_counts    = 30,
                        min_abs_intensity = 100,
                        strong_snr        = 10,
                        impute            = TRUE) {
  # Linear-impute INTERIOR NA gaps on the RT grid before dropping NAs.
  # rule = 1 (interior only): leading/trailing NAs stay NA and get
  # dropped, so the tail-noise / baseline estimate is not contaminated
  # by flat edge-extrapolation (rule = 2 caused false below_detection
  # regressions on sparse consensus EICs — see impute impact test).
  # `real0` tracks which scans are genuine (non-imputed) so the
  # detection guard can require real-scan support and not certify a
  # peak that is purely interpolated between 2-3 isolated points.
  real0 <- !is.na(intensity)
  isok <- !is.na(rt) & !is.na(intensity)
  if (impute && sum(isok) >= 2L) {
    fill <- !is.na(rt)
    intensity[fill] <- stats::approx(rt[isok], intensity[isok],
                                     xout = rt[fill], rule = 1)$y
  }
  ok <- !is.na(intensity); rt <- rt[ok]; intensity <- intensity[ok]
  real0 <- real0[ok]
  if (length(rt) < 5L)
    return(list(apex_rt = NA_real_, maxIntensity = NA_real_,
                areaUnderTic = NA_real_, xicFwhm = NA_real_,
                snr = NA_real_, baseline = NA_real_, noise = NA_real_,
                n_scans_in_shape = 0L, detection = "no_data"))
  sm <- as.numeric(stats::filter(intensity, rep(1, 3) / 3, sides = 2))
  sm[is.na(sm)] <- intensity[is.na(sm)]
  search_idx <- which(rt >= rt_prior - apex_search_pad &
                      rt <= rt_prior + apex_search_pad)
  if (!length(search_idx))
    return(list(apex_rt = NA_real_, maxIntensity = NA_real_,
                areaUnderTic = NA_real_, xicFwhm = NA_real_,
                snr = NA_real_, baseline = NA_real_, noise = NA_real_,
                n_scans_in_shape = 0L, detection = "no_data"))
  apex_local <- search_idx[which.max(intensity[search_idx])]
  apex_rt    <- rt[apex_local]; apex_int <- intensity[apex_local]

  tail_idx <- which(abs(rt - apex_rt) > tail_pad)
  if (length(tail_idx) < 10L) {
    n <- length(rt); tail_idx <- c(seq_len(round(0.25 * n)),
                                    (n - round(0.25 * n) + 1L):n)
  }
  tail_int <- intensity[tail_idx]
  tail_low <- tail_int[tail_int <= quantile(tail_int, 0.5, na.rm = TRUE)]
  baseline <- median(tail_low, na.rm = TRUE)
  q15 <- as.numeric(quantile(tail_low, 0.15, na.rm = TRUE))
  q85 <- as.numeric(quantile(tail_low, 0.85, na.rm = TRUE))
  raw_noise <- max(mad(tail_low, na.rm = TRUE), (q85 - q15) / 2, na.rm = TRUE)
  noise_at_floor <- !is.finite(raw_noise) || raw_noise <= 0
  noise <- max(raw_noise, 1, na.rm = TRUE)
  snr   <- (apex_int - baseline) / noise

  v_idx <- tryCatch(MsCoreUtils::valleys(sm), error = function(e) integer(0))
  if (length(v_idx) >= 2L) {
    left  <- suppressWarnings(max(v_idx[v_idx <= apex_local], na.rm = TRUE))
    right <- suppressWarnings(min(v_idx[v_idx >= apex_local], na.rm = TRUE))
    if (!is.finite(left))  left  <- 1L
    if (!is.finite(right)) right <- length(rt)
  } else {
    left <- apex_local; right <- apex_local
    while (left  > 1L         && sm[left  - 1L] > baseline + noise) left  <- left  - 1L
    while (right < length(rt) && sm[right + 1L] > baseline + noise) right <- right + 1L
  }
  pk_rt <- rt[left:right]; pk_int <- intensity[left:right]
  area  <- if (length(pk_rt) >= 2L)
             sum(diff(pk_rt) * (pk_int[-1] + pk_int[-length(pk_int)]) / 2)
           else NA_real_
  hm    <- (apex_int - baseline) / 2 + baseline
  above <- which(pk_int >= hm)
  fwhm  <- if (length(above) >= 2L) diff(range(pk_rt[above])) else NA_real_

  nb_lo  <- max(1L, apex_local - 2L); nb_hi <- min(length(rt), apex_local + 2L)
  n_shape <- sum(intensity[nb_lo:nb_hi] >= 0.3 * apex_int, na.rm = TRUE)
  # Real-scan support: count genuine (non-imputed) scans within ±10 s of
  # the apex. A real peak is sampled by several scans across its width;
  # an interpolation artifact (2-3 isolated real points joined by a line)
  # has only 1 real scan near the apex and is rejected here.
  n_real_near_apex <- sum(real0 & abs(rt - apex_rt) <= 10)
  # Signal concentration: fraction of above-baseline intensity within
  # ±12 s of the detected apex. A real chromatographic peak concentrates
  # its signal at the apex (~>0.5); broad noise scatter spreads it out,
  # and a marginal bump at a wrong RT prior (the real peak eluting
  # elsewhere) has most of its above-baseline mass away from the apex —
  # both give a low fraction and are rejected.
  ab_bl <- pmax(intensity - baseline, 0)
  near_apex <- abs(rt - apex_rt) <= 12
  concentration <- sum(ab_bl[near_apex], na.rm = TRUE) /
                   max(sum(ab_bl, na.rm = TRUE), 1)

  # Quality guards first (genuine non-detections), THEN the absolute-intensity
  # floor as a BACKSTOP. A fixed count floor is vendor-miscalibrated (Waters/Sciex
  # run on a lower scale, so a 1000-count floor wrongly kills real low-abundance
  # peaks), so it must not veto a peak that already cleared SNR + shape +
  # concentration + real-scan support with strong SNR. The floor therefore applies
  # only to marginal peaks (snr < strong_snr) or when noise sat at the floor (so
  # SNR is unreliable). This rescues strong, well-shaped, well-placed low-abundance
  # peaks that DDA correctly triggered MS2 on. See QMD 8 "MS2 without MS1".
  reject_reason <-
        if (apex_int <= baseline) "apex<=baseline"
        else if (snr < snr_threshold) "snr<thresh"
        else if (n_shape < 2L) "n_shape<2"
        else if (n_real_near_apex < 3L) "n_real_near_apex<3"
        else if (concentration < 0.45) "concentration<0.45"
        else if ((snr < strong_snr || noise_at_floor) &&
                 apex_int < min_abs_intensity) "apex_int<min_abs"
        else if (noise_at_floor && apex_int < abs_min_counts) "noise_floor"
        else NA_character_
  detection <- if (is.na(reject_reason)) "detected" else "below_detection"

  list(apex_rt = apex_rt,
       maxIntensity = if (detection == "detected") apex_int else NA_real_,
       areaUnderTic = if (detection == "detected") area     else NA_real_,
       xicFwhm     = if (detection == "detected") fwhm     else NA_real_,
       snr = snr, baseline = baseline, noise = noise,
       n_scans_in_shape = n_shape,
       n_real_near_apex = n_real_near_apex, concentration = concentration,
       reject_reason = reject_reason, detection = detection)
}

# ───────────────────────────────────────────────────────────────────
# 2. Main pipeline
# ───────────────────────────────────────────────────────────────────

run_annotated_eic_pipeline <- function(feat_rts_path,
                                       identity_path,
                                       splines, anchor_table,
                                       output_path,
                                       labs       = c("afekta", "cembio", "hmgu", "icl"),
                                       mixtures,
                                       ppm_tol    = 20,
                                       pad_wide   = 60,
                                       rt_cluster_window = 30,
                                       eic_dump_path = NULL,
                                       # Per-lab absolute intensity floor (counts). Derived from
                                       # the per-lab clean-vs-scatter apex distributions: Agilent
                                       # labs (afekta, cembio) and icl separate around ~1000
                                       # counts; hmgu (Sciex, lower intensity scale) at ~100.
                                       lab_intensity_floor = c(afekta = 1000, cembio = 1000,
                                                               hmgu = 100, icl = 1000)) {

  message("[1/8] Loading inputs...")
  eic_candidates <- read.csv(feat_rts_path, stringsAsFactors = FALSE,
                             colClasses = c(mixture = "character"))
  identity_df    <- read.csv(identity_path, stringsAsFactors = FALSE,
                             colClasses = c(mixture = "character"))
  # Normalize mixture IDs to underscore form. Defensive: ms1_anchor.R
  # already does this on its output, but annotation_identity.csv is
  # the SIRIUS-producer's CSV and is still in dot form ("2.1", "2.10").
  eic_candidates$mixture <- gsub("\\.", "_", eic_candidates$mixture)
  identity_df$mixture    <- gsub("\\.", "_", identity_df$mixture)
  # Same for the requested mixtures parameter so the per-lab file lookup
  # works regardless of which form the caller passed.
  mixtures <- gsub("\\.", "_", as.character(mixtures))
  inv_splines    <- build_inverse_splines(splines, anchor_table)

  message("[2/8] Forward-aligning MS2 RTs...")
  eic_candidates <- eic_candidates |>
    dplyr::rowwise() |>
    dplyr::mutate(rt_aligned = apply_fwd(rt_sec, lab, polarity,
                                          splines, anchor_table)) |>
    dplyr::ungroup()

  # [3] Cross-lab consensus clustering — historically a 30 s greedy
  # 1-D pass over rt_aligned per (mixture, compound). After the
  # ms1_anchor.R Stage 2 refactor, every row's rt_sec was already
  # snapped to its inverse-aligned consensus candidate, so rt_aligned
  # values for one (mixture, compound) now converge; this step is now
  # a no-op (always produces one cluster per compound). Kept as an
  # audit pass — a non-trivial output here would indicate Stage 2
  # missed a case.
  message("[3/8] Consensus clustering (audit pass — expected to be no-op post-refactor)...")
  eic_candidates <- eic_candidates |>
    dplyr::group_by(mixture, compound) |>
    dplyr::mutate(rt_cluster = cluster_rts_greedy(rt_aligned,
                                                  window = rt_cluster_window)) |>
    dplyr::ungroup()

  cluster_scores <- eic_candidates |>
    dplyr::group_by(mixture, compound, rt_cluster) |>
    dplyr::summarise(n_labs       = dplyr::n_distinct(lab),
                     n_polarities = dplyr::n_distinct(polarity),
                     n_features   = dplyr::n(),
                     sum_cosmic   = sum(cosmic, na.rm = TRUE),
                     rt_center_aligned = median(rt_aligned),
                     rt_center    = median(rt_sec),
                     .groups      = "drop") |>
    dplyr::mutate(score = n_labs * 10 + n_polarities * 5 + sum_cosmic)

  winning_cluster <- cluster_scores |>
    dplyr::group_by(mixture, compound) |>
    dplyr::slice_max(score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(mixture, compound,
                  winning_cluster      = rt_cluster,
                  consensus_rt_aligned = rt_center_aligned,
                  consensus_rt         = rt_center,
                  consensus_n_labs     = n_labs,
                  consensus_n_pol      = n_polarities,
                  consensus_score      = score)

  # [3.1] Isobaric-collision resolution — AUDIT PASS.
  # After the ms1_anchor.R Stage 2 refactor, cross-compound isobaric
  # collisions are resolved upstream (per (mixture) pairs with any
  # adduct m/z within ppm_tol AND consensus RT within isobar_rt_window
  # → drop the lower-scoring compound). The lower-scoring side's rows
  # never reach this pipeline, so this duplicate check should always
  # find 0 collisions to resolve. A non-zero count would indicate
  # Stage 2 missed something — investigate.
  isobar_rt_window <- 5  # seconds — must match ms1_anchor.R::anchor_rts
  per_adduct_mz <- eic_candidates |>
    dplyr::group_by(mixture, compound, adduct) |>
    dplyr::summarise(ion_mz = stats::median(ionMass, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(is.finite(ion_mz), ion_mz > 0)
  compound_aux <- eic_candidates |>
    dplyr::group_by(mixture, compound) |>
    dplyr::summarise(sum_cosmic_all = sum(cosmic, na.rm = TRUE),
                     .groups = "drop")
  wc_keyed <- winning_cluster |>
    dplyr::inner_join(compound_aux, by = c("mixture", "compound"))

  # Pairwise collision detection within each mixture: two compounds collide
  # if ANY of their adducts' m/z fall within ppm_tol AND their consensus_rt
  # within isobar_rt_window.
  collide_check <- function(wcx) {
    n <- nrow(wcx)
    if (n < 2L) return(rep(TRUE, n))
    mz_lookup <- split(per_adduct_mz$ion_mz[per_adduct_mz$mixture == wcx$mixture[1]],
                       per_adduct_mz$compound[per_adduct_mz$mixture == wcx$mixture[1]])
    keep <- rep(TRUE, n)
    ord  <- order(wcx$consensus_score, wcx$sum_cosmic_all, decreasing = TRUE)
    for (i in seq_len(n - 1)) {
      ii <- ord[i]; if (!keep[ii]) next
      mz_i <- mz_lookup[[wcx$compound[ii]]]
      if (is.null(mz_i) || !length(mz_i)) next
      for (j in (i + 1):n) {
        jj <- ord[j]; if (!keep[jj]) next
        if (abs(wcx$consensus_rt[ii] - wcx$consensus_rt[jj]) > isobar_rt_window)
          next
        mz_j <- mz_lookup[[wcx$compound[jj]]]
        if (is.null(mz_j) || !length(mz_j)) next
        # Any-adduct any-adduct m/z match within ppm_tol?
        ppm <- outer(mz_i, mz_j, function(a, b) abs(a - b) / a * 1e6)
        if (any(is.finite(ppm) & ppm <= ppm_tol)) {
          keep[jj] <- FALSE   # lower-scoring compound loses
        }
      }
    }
    keep
  }
  wc_keyed$.keep <- unlist(by(wc_keyed, wc_keyed$mixture, collide_check),
                            use.names = FALSE)
  isobar_diag <- wc_keyed |>
    dplyr::group_by(mixture) |>
    dplyr::summarise(n_collisions_resolved = sum(!.keep),
                     .groups = "drop") |>
    dplyr::filter(n_collisions_resolved > 0)
  isobar_dropped <- wc_keyed |> dplyr::filter(!.keep) |>
    dplyr::select(mixture, compound, consensus_rt, consensus_score)
  n_collisions <- sum(!wc_keyed$.keep)
  message(sprintf("[3.1] Isobaric collisions audit — found at downstream pass: %d (expected 0; real resolution happens in ms1_anchor.R Stage 2)",
                  n_collisions))
  winning_cluster <- wc_keyed |> dplyr::filter(.keep) |>
    dplyr::select(-.keep, -sum_cosmic_all, -consensus_score)

  message("[4/8] Building per-row priors (inverse-aligned to lab-local)...")
  eic_params <- eic_candidates |>
    dplyr::inner_join(winning_cluster |>
                      dplyr::select(mixture, compound, winning_cluster),
                      by = c("mixture", "compound")) |>
    dplyr::filter(rt_cluster == winning_cluster) |>
    dplyr::group_by(mixture, lab, polarity, compound, adduct) |>
    dplyr::summarise(ionMass   = median(ionMass),
                     rt_sec    = median(rt_sec),
                     rt_spread = max(rt_sec) - min(rt_sec),
                     .groups   = "drop")

  adduct_inv <- eic_candidates |>
    dplyr::distinct(mixture, polarity, compound, adduct, ionMass) |>
    dplyr::group_by(mixture, polarity, compound, adduct) |>
    dplyr::summarise(ionMass = median(ionMass), .groups = "drop")

  eic_consensus <- adduct_inv |>
    dplyr::inner_join(winning_cluster |>
                      dplyr::select(mixture, compound, consensus_rt_aligned),
                      by = c("mixture", "compound")) |>
    dplyr::rename(rt_sec_aligned = consensus_rt_aligned)

  peak_table <- tidyr::expand_grid(
      lab             = labs,
      mix_pol_cpd_add = eic_consensus |>
                          dplyr::distinct(mixture, polarity, compound, adduct)) |>
    tidyr::unnest(mix_pol_cpd_add) |>
    dplyr::left_join(eic_params,
                     by = c("lab", "mixture", "polarity", "compound", "adduct")) |>
    dplyr::left_join(eic_consensus,
                     by = c("mixture", "polarity", "compound", "adduct"),
                     suffix = c("", "_cons")) |>
    dplyr::rowwise() |>
    dplyr::mutate(rt_sec_cons_local =
                    apply_inv(rt_sec_aligned, lab, polarity,
                              inv_splines, anchor_table)) |>
    dplyr::ungroup() |>
    dplyr::mutate(source = dplyr::if_else(is.na(ionMass), "consensus", "own"),
                  rt = dplyr::coalesce(rt_sec, rt_sec_cons_local),
                  mz = dplyr::coalesce(ionMass, ionMass_cons))

  message("[5/8] Building lab_files (mzML paths)...")
  lab_files <- build_lab_files(labs, mixtures)
  peak_table <- peak_table |>
    dplyr::left_join(lab_files |>
                       dplyr::select(lab, mixture, polarity, path),
                     by = c("lab", "mixture", "polarity")) |>
    dplyr::filter(!is.na(path), !is.na(rt), !is.na(mz))

  # Cosmic for the per-row prior demotion guard.
  cosmic_per_row <- identity_df |>
    dplyr::group_by(mixture, lab, compound, polarity) |>
    dplyr::summarise(cosmic = suppressWarnings(
                       max(confidenceExactMatch, na.rm = TRUE)),
                     .groups = "drop") |>
    dplyr::mutate(cosmic = ifelse(is.infinite(cosmic), NA_real_, cosmic))

  rt_consensus_per_cpd <- winning_cluster |>
    dplyr::select(mixture, compound, rt_aligned = consensus_rt_aligned) |>
    tidyr::expand_grid(polarity = c("pos", "neg"), lab = labs) |>
    dplyr::rowwise() |>
    dplyr::mutate(rt_consensus = apply_inv(rt_aligned, lab, polarity,
                                            inv_splines, anchor_table)) |>
    dplyr::ungroup() |>
    dplyr::select(mixture, compound, polarity, lab, rt_consensus)

  prior_per_row <- peak_table |>
    dplyr::select(mixture, lab, polarity, compound, adduct, source) |>
    dplyr::left_join(eic_params |>
                       dplyr::select(mixture, lab, polarity, compound, adduct,
                                     rt_own = rt_sec, rt_spread_own = rt_spread),
                     by = c("mixture","lab","polarity","compound","adduct")) |>
    dplyr::left_join(rt_consensus_per_cpd,
                     by = c("mixture","lab","polarity","compound")) |>
    dplyr::left_join(cosmic_per_row,
                     by = c("mixture","lab","compound","polarity")) |>
    dplyr::mutate(prior_kind = dplyr::case_when(
                    source == "own" & !is.na(rt_own) & !is.na(cosmic) &
                      cosmic >= 0.1 & !is.na(rt_consensus) &
                      abs(rt_own - rt_consensus) <= 30 ~ "own",
                    TRUE ~ "consensus"),
                  rt_prior  = ifelse(prior_kind == "own", rt_own, rt_consensus),
                  apex_pad  = ifelse(prior_kind == "own",
                                     pmax(8, tidyr::replace_na(rt_spread_own, 0) / 2 + 5),
                                     15))

  message("[6/8] Loading Spectra + chromExtract...")
  sps_all <- Spectra(unique(peak_table$path), backend = MsBackendMzR())
  sps_all <- sps_all[msLevel(sps_all) == 1L]
  canon   <- unique(dataOrigin(sps_all))
  peak_table$dataOrigin <- canon[match(basename(peak_table$path),
                                       basename(canon))]
  stopifnot(all(!is.na(peak_table$dataOrigin)))

  pt_for_extract <- peak_table |>
    dplyr::transmute(mzMin = mz * (1 - ppm_tol / 1e6),
                     mzMax = mz * (1 + ppm_tol / 1e6),
                     rtMin = pmax(0, rt - pad_wide),
                     rtMax = rt + pad_wide,
                     msLevel = 1L,
                     dataOrigin, lab, mixture, polarity,
                     compound, adduct, source,
                     mz_center = mz, rt_center = rt)
  chr_base <- Chromatograms(sps_all)
  # peak.table intentionally carries `polarity` (we want it in chromData);
  # chromExtract warns that it replaces the existing column. That
  # replacement is the desired behaviour, so muffle ONLY that warning.
  chr_all  <- withCallingHandlers(
    chromExtract(chr_base, peak.table = pt_for_extract,
                 by = c("msLevel", "dataOrigin")),
    warning = function(w) {
      if (grepl("already exist in", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  cd_all <- as.data.frame(chromData(chr_all))
  pkData <- peaksData(chr_all)

  message("[7/8] Running SNR-aware peak picker per row...")
  stopifnot(nrow(cd_all) == nrow(prior_per_row))
  metrics_rows <- vector("list", nrow(cd_all))
  for (i in seq_len(nrow(cd_all))) {
    pd <- pkData[[i]]
    prior_i <- prior_per_row[i, ]
    if (is.null(pd) || !nrow(pd)) {
      metrics_rows[[i]] <- list(apex_rt = NA_real_, maxIntensity = NA_real_,
                                areaUnderTic = NA_real_, xicFwhm = NA_real_,
                                snr = NA_real_, baseline = NA_real_,
                                noise = NA_real_, n_scans_in_shape = 0L,
                                detection = "no_data")
      next
    }
    lab_i  <- cd_all$lab[i]
    floor_i <- if (!is.na(lab_intensity_floor[lab_i]))
                 lab_intensity_floor[[lab_i]] else 100
    metrics_rows[[i]] <- extract_one(pd[, "rtime"], pd[, "intensity"],
                                      rt_prior = cd_all$rt_center[i],
                                      apex_search_pad = prior_i$apex_pad,
                                      min_abs_intensity = floor_i)
  }

  metrics_df <- dplyr::bind_cols(
    cd_all |> dplyr::select(lab, mixture, polarity, compound, adduct, source),
    dplyr::bind_rows(lapply(metrics_rows, as.data.frame)),
    prior_per_row |> dplyr::select(prior_kind))

  # Optional: dump per-row EICs + priors + metrics for offline diagnostics
  # (e.g. plotting rescued peaks). Not written in normal renders.
  if (!is.null(eic_dump_path)) {
    saveRDS(list(cd_all = cd_all, pkData = pkData,
                 prior_per_row = prior_per_row, metrics_df = metrics_df),
            eic_dump_path)
    message("  [debug] EIC dump written to ", eic_dump_path)
  }

  message("[8/8] Writing annotation_chrom_metrics.csv...")
  annotation_chrom_metrics <- metrics_df |>
    dplyr::select(mixture, lab, polarity, compound, adduct, source,
                  areaUnderTic, xicFwhm, maxIntensity, apex_rt,
                  snr, baseline, noise, n_scans_in_shape,
                  n_real_near_apex, concentration, reject_reason,
                  prior_kind, detection)
  write.csv(annotation_chrom_metrics, output_path, row.names = FALSE)
  message(sprintf("  Wrote %d rows to %s", nrow(annotation_chrom_metrics),
                  output_path))
  invisible(list(annotation_chrom_metrics = annotation_chrom_metrics,
                 winning_cluster          = winning_cluster,
                 cluster_scores           = cluster_scores,
                 isobar_diagnostics       = list(
                   n_collisions_resolved = n_collisions,
                   per_mixture           = isobar_diag,
                   dropped               = isobar_dropped)))
}
