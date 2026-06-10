## ms1_anchor.R — function library, no side effects.
##
## Replaces each lab's MS2-trigger RT (rt_sec) in feat_rts with the MS1
## chromatographic peak apex matching the cross-lab consensus, using
## the SAME MS1 file the downstream EIC extraction reads from.
##
## Architecture (3 stages):
##   STAGE 1 — within-lab candidate finding.
##     For each MS2 trigger, extract the MS1 EIC at the precursor m/z,
##     run find_all_peaks (SNR ≥ 5), validate each peak with the SAME
##     stacked guards Stage 3 uses (per-lab intensity floor, peak shape
##     ≥ 2 scans, real-scan support, signal concentration ≥ 0.45), then
##     keep peaks ≥ candidate_frac of the surviving max as candidates.
##     No cross-lab look. Multiple candidates per row are allowed.
##
##   STAGE 2 — cross-lab consensus in aligned space.
##     Forward-align every candidate to the consensus axis (via NAPS
##     splines), then per (mixture, compound) apply the similarity rule:
##     for each candidate, count distinct (lab, polarity) within
##     ±agreement_window of it; among ALL candidates that reach the max
##     support level, take the median of their aligned RTs as the
##     centre; consensus_rt_aligned = median of all candidates within
##     ±agreement_window of that centre. Then cross-compound isobaric
##     resolution: pairs of compounds in the same mixture whose
##     consensus RTs are within isobar_rt_window AND ANY of their
##     adducts' m/z fall within ppm_tol → drop the lower-scoring.
##
##   STAGE 2 SNAP — inverse-align consensus per row, snap to nearest
##     candidate. consensus_rt_aligned → row's lab-local target RT via
##     inverse spline; each row's rt_sec is replaced by the candidate
##     nearest to that target. Rows whose compound was isobar-dropped
##     are dropped entirely. Rows whose compound had no cross-lab
##     consensus fall back to the strongest candidate.
##
## Sourced by `5_downstream_analysis/3_rt_alignment.qmd` §P2.9.
## Public entry point: `anchor_rts(feat, splines, anchor_table, ...)`.
##
## Stage 3 (per-lab EIC extraction with the same stacked guards) lives
## in `annotated_eic_pipeline.R::extract_one` and is sourced here so
## Stage 1 reuses the same guard logic — keeping "what counts as a
## peak" consistent across the whole pipeline.
suppressPackageStartupMessages({
  library(dplyr); library(here); library(readxl); library(tidyr)
  library(Spectra); library(Chromatograms)
})

# Null-coalescing operator (defined for R < 4.4 compatibility).
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# Sourced for: extract_one (Stage 3's peak-quality guards, reused as
# Stage 1's candidate-validation pass), apply_fwd, apply_inv,
# build_inverse_splines.
source(here::here("5_downstream_analysis", "annotated_eic_pipeline.R"))

# ── helpers ─────────────────────────────────────────────────────────
.mix_us <- function(m) gsub("\\.", "_", m)

# Find the MS1 mzML file for a given (lab, mixture, polarity).
# For separate-MS2 labs (afekta, cembio, icl) this returns the MS1-only
# injection (no _MS2 in filename). For hmgu it returns the combined
# MS1+MS2 file.
ms1_path_for <- function(lab, mixture, polarity, study_group = "HE") {
  study_path <- here::here("1_preprocessing", lab, study_group)
  mzml_dir   <- file.path(study_path, "mzml")
  if (!dir.exists(mzml_dir)) return(NA_character_)
  seq_xlsx <- file.path(study_path, sprintf("seq_%s.xlsx", polarity))
  if (!file.exists(seq_xlsx)) return(NA_character_)
  seq <- as.data.frame(read_xlsx(seq_xlsx))
  seq$filename <- paste0(seq$`Data File`, ".mzML")
  seq <- seq[seq$filename %in% list.files(mzml_dir), , drop = FALSE]
  if (!nrow(seq)) return(NA_character_)
  has_separate_ms2 <- any(grepl("_MS2", list.files(mzml_dir), fixed = TRUE))
  if (has_separate_ms2)
    seq <- seq[!grepl("_MS2", seq$filename, fixed = TRUE), , drop = FALSE]
  seq$mixture_us <- gsub("\\.", "_",
                         sub("target$", "",
                             sub(".*_", "",
                                 sub("_MS2$", "", seq$`Data File`))))
  hit <- seq[seq$mixture_us == .mix_us(mixture), , drop = FALSE]
  if (!nrow(hit)) return(NA_character_)
  file.path(mzml_dir, hit$filename[1])
}

# Linear NA-fill of an EIC's INTERIOR gaps along the RT axis.
impute_eic_linear <- function(rt, intensity) {
  ok <- !is.na(rt) & !is.na(intensity)
  if (sum(ok) < 2L) return(intensity)
  out <- intensity
  out[!is.na(rt)] <- stats::approx(rt[ok], intensity[ok],
                                   xout = rt[!is.na(rt)], rule = 1)$y
  out
}

# SNR-aware multi-peak finder over a chromatogram.
find_all_peaks <- function(eic, snr_threshold = 5) {
  empty <- data.frame(apex_rt = numeric(0), maxIntensity = numeric(0),
                      snr = numeric(0))
  if (is.null(eic) || nrow(eic) < 10L) return(empty)
  rt <- eic$rt; intensity <- eic$intensity
  intensity <- impute_eic_linear(rt, intensity)
  ok <- !is.na(rt) & !is.na(intensity)
  rt <- rt[ok]; intensity <- intensity[ok]
  if (length(rt) < 10L) return(empty)
  sm <- as.numeric(stats::filter(intensity, rep(1, 3) / 3, sides = 2))
  sm[is.na(sm)] <- intensity[is.na(sm)]
  is_lm <- c(FALSE,
             sm[-c(1, length(sm))] > sm[-c(length(sm) - 1, length(sm))] &
             sm[-c(1, length(sm))] > sm[-c(1, 2)],
             FALSE)
  cands <- which(is_lm)
  if (!length(cands)) return(empty)
  out <- lapply(cands, function(a) {
    tail_idx <- which(abs(rt - rt[a]) > 25)
    if (length(tail_idx) < 10L) {
      n <- length(rt); tail_idx <- c(seq_len(round(0.25 * n)),
                                      (n - round(0.25 * n) + 1L):n)
    }
    ti <- intensity[tail_idx]
    tl <- ti[ti <= quantile(ti, 0.5, na.rm = TRUE)]
    bl <- median(tl, na.rm = TRUE)
    q15 <- as.numeric(quantile(tl, 0.15, na.rm = TRUE))
    q85 <- as.numeric(quantile(tl, 0.85, na.rm = TRUE))
    raw_noise <- max(mad(tl, na.rm = TRUE), (q85 - q15) / 2, na.rm = TRUE)
    nu <- max(raw_noise, 1, na.rm = TRUE)
    list(apex_rt = rt[a], maxIntensity = intensity[a],
         snr = (intensity[a] - bl) / nu)
  })
  r <- dplyr::bind_rows(lapply(out, as.data.frame))
  r <- r[r$snr >= snr_threshold & r$maxIntensity > 0, , drop = FALSE]
  r[order(-r$maxIntensity), ]
}

# Validate one peak with Stage 3's stacked guards: drops peaks the
# downstream extract_one would refuse (intensity floor, peak shape ≥
# 2 scans, real-scan support, signal concentration ≥ 0.45). Returns
# TRUE iff extract_one with a tight apex_search_pad lands "detected"
# on this candidate's apex.
validate_peak <- function(apex_rt, eic, lab_floor) {
  res <- extract_one(eic$rt, eic$intensity,
                     rt_prior          = apex_rt,
                     apex_search_pad   = 3,
                     min_abs_intensity = lab_floor)
  isTRUE(res$detection == "detected")
}

# ── public entry point ─────────────────────────────────────────────
anchor_rts <- function(feat,
                       splines, anchor_table,
                       labs            = c("afekta", "cembio", "hmgu", "icl"),
                       ppm_tol         = 20,
                       search_window   = 60,
                       snr_threshold   = 5,
                       gradient_max_rt = 10000,
                       # Per-lab absolute intensity floor (counts). Same values
                       # Stage 3 uses; tuned from the per-lab clean-vs-scatter
                       # apex distributions (Agilent / Waters ≈ 1000, Sciex ≈ 100).
                       lab_intensity_floor = c(afekta = 1000, cembio = 1000,
                                                hmgu = 100, icl = 1000),
                       # A peak counts as a candidate if it reaches this fraction
                       # of the row's strongest peak that survived the Stage 3 guards.
                       candidate_frac  = 0.5,
                       # Window (seconds in the ALIGNED consensus axis) for the
                       # similarity rule that defines cross-lab agreement.
                       agreement_window = 8,
                       # Window (aligned seconds) for cross-compound isobaric
                       # collision detection.
                       isobar_rt_window = 5,
                       # Labs whose MS2 fires far from the real peak (cembio
                       # inclusion-list schedule); search the full gradient.
                       wide_search_labs = "cembio",
                       verbose         = TRUE) {

  # Normalize mixture IDs to underscore form ("2_1", "2_10") — protects
  # against silent dot/underscore mixing downstream (numeric coercion of
  # "2.10" → 2.1 collides with "2.1"). Doing it once here propagates to
  # feat_rts_ms1_anchored.csv and on to annotation_chrom_metrics.csv.
  if ("mixture" %in% names(feat))
    feat$mixture <- gsub("\\.", "_", as.character(feat$mixture))

  in_scope <- feat |> dplyr::filter(lab %in% labs)
  pass_through <- feat |> dplyr::filter(!(lab %in% labs))
  if (verbose)
    message(sprintf("[anchor] in-scope rows: %d  (labs: %s) | pass-through: %d",
                    nrow(in_scope), paste(labs, collapse = ","),
                    nrow(pass_through)))

  if (!nrow(in_scope)) {
    result <- feat
    attr(result, "diagnostics") <- list(n_input = 0L, n_kept = 0L,
                                         n_dropped = 0L)
    return(result)
  }

  # SIRIUS's per-feature ionMass drifts with mass accuracy; use the
  # group median so all features of one ion share a single EIC.
  if (all(c("compound", "adduct") %in% names(in_scope))) {
    in_scope <- in_scope |>
      dplyr::group_by(lab, mixture, polarity, compound, adduct) |>
      dplyr::mutate(ionMass = stats::median(ionMass, na.rm = TRUE)) |>
      dplyr::ungroup()
  }

  win_for <- function(lab_name) {
    if (length(search_window) == 1L) return(search_window)
    if (!is.null(names(search_window)) &&
        lab_name %in% names(search_window))
      return(search_window[[lab_name]])
    60
  }

  # 1. Resolve MS1 paths per (lab, mixture, polarity).
  groups <- in_scope |> dplyr::distinct(lab, mixture, polarity) |>
    dplyr::arrange(lab, mixture, polarity)
  groups$path <- mapply(ms1_path_for, groups$lab,
                        groups$mixture, groups$polarity)
  groups <- groups |> dplyr::filter(!is.na(path), file.exists(path))
  if (verbose)
    message(sprintf("[anchor] resolved MS1 paths: %d / %d (lab, mixture, polarity)",
                    nrow(groups),
                    nrow(dplyr::distinct(in_scope, lab, mixture, polarity))))

  in_scope <- in_scope |>
    dplyr::left_join(groups |> dplyr::select(lab, mixture, polarity, path),
                     by = c("lab", "mixture", "polarity"))
  bad_path <- is.na(in_scope$path)
  bad_mz   <- is.na(in_scope$ionMass) |
              !is.finite(in_scope$ionMass) | in_scope$ionMass <= 0
  if (verbose)
    message(sprintf("[anchor] pre-extract drops: %d (no_path=%d, no_mz=%d)",
                    sum(bad_path | bad_mz), sum(bad_path), sum(bad_mz)))
  workable <- in_scope |> dplyr::filter(!bad_path, !bad_mz)

  # 2. Build peak.table: one row per UNIQUE (lab, mixture, polarity, ionMass).
  unique_eics <- workable |>
    dplyr::group_by(lab, mixture, polarity, path, ionMass) |>
    dplyr::summarise(trigger_min = min(rt_sec, na.rm = TRUE),
                     trigger_max = max(rt_sec, na.rm = TRUE),
                     .groups     = "drop") |>
    dplyr::mutate(
      win    = vapply(lab, win_for, numeric(1)),
      wide   = lab %in% wide_search_labs,
      rtMin  = ifelse(wide, 0, pmax(0, trigger_min - win - 30)),
      rtMax  = ifelse(wide, gradient_max_rt,
                      pmin(gradient_max_rt, trigger_max + win + 30)),
      mz_key = sprintf("%.4f", ionMass))

  if (verbose)
    message(sprintf("[anchor] unique EICs to extract: %d  (vs %d rows → %.1fx dedupe)",
                    nrow(unique_eics), nrow(workable),
                    nrow(workable) / pmax(1, nrow(unique_eics))))

  # 3. Load Spectra (MS1 only), canonicalize paths.
  sps_all <- Spectra(unique(unique_eics$path), backend = MsBackendMzR())
  sps_all <- sps_all[msLevel(sps_all) == 1L]
  canon <- unique(dataOrigin(sps_all))
  unique_eics$dataOrigin <- canon[match(basename(unique_eics$path),
                                         basename(canon))]
  stopifnot(all(!is.na(unique_eics$dataOrigin)))

  # 4. chromExtract over all (mz, file) windows in one call.
  pt <- unique_eics |>
    dplyr::transmute(mzMin     = ionMass * (1 - ppm_tol / 1e6),
                     mzMax     = ionMass * (1 + ppm_tol / 1e6),
                     rtMin, rtMax,
                     msLevel   = 1L,
                     dataOrigin,
                     lab, mixture, polarity, ionMass, mz_key)
  if (verbose) message("[anchor] running chromExtract over ",
                       nrow(pt), " (mz, file) windows...")
  chr_all <- withCallingHandlers(
    chromExtract(Chromatograms(sps_all), peak.table = pt,
                 by = c("msLevel", "dataOrigin")),
    warning = function(w) {
      if (grepl("already exist in", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  cd <- as.data.frame(chromData(chr_all))
  pkData <- peaksData(chr_all)

  cd$eic_idx <- seq_len(nrow(cd))
  eic_lookup <- cd |>
    dplyr::select(lab, mixture, polarity, mz_key, eic_idx)

  workable <- workable |>
    dplyr::mutate(mz_key = sprintf("%.4f", ionMass)) |>
    dplyr::left_join(eic_lookup,
                     by = c("lab", "mixture", "polarity", "mz_key"))

  # 5. STAGE 1 — per-row candidate finding with Stage 3 guards.
  if (verbose) message("[anchor] Stage 1: finding candidates per row...")
  cand_rts     <- vector("list", nrow(workable))
  strongest_rt <- rep(NA_real_, nrow(workable))
  status_v     <- rep("", nrow(workable))
  for (ri in seq_len(nrow(workable))) {
    idx <- workable$eic_idx[ri]
    if (is.na(idx)) { status_v[ri] <- "no_eic"; next }
    pd <- pkData[[idx]]
    if (is.null(pd) || !nrow(pd)) { status_v[ri] <- "no_data"; next }
    eic <- data.frame(rt = pd[, "rtime"], intensity = pd[, "intensity"])
    pks <- find_all_peaks(eic, snr_threshold = snr_threshold)
    if (!nrow(pks)) { status_v[ri] <- "no_peak"; next }

    floor_i <- if (!is.na(lab_intensity_floor[workable$lab[ri]]))
                 lab_intensity_floor[[workable$lab[ri]]] else 100
    keep_validated <- vapply(seq_len(nrow(pks)), function(j)
      validate_peak(pks$apex_rt[j], eic, floor_i), logical(1))
    pks_val <- pks[keep_validated, , drop = FALSE]
    if (!nrow(pks_val)) { status_v[ri] <- "no_peak_passes_guards"; next }
    keep <- pks_val$maxIntensity >= candidate_frac * max(pks_val$maxIntensity)
    cand_rts[[ri]]   <- pks_val$apex_rt[keep]
    strongest_rt[ri] <- pks_val$apex_rt[which.max(pks_val$maxIntensity)]
    status_v[ri]     <- "candidates_found"
  }
  workable$candidate_rts <- vapply(cand_rts, function(x)
    if (length(x)) paste(sprintf("%.3f", x), collapse = ";") else "",
    character(1))
  workable$n_candidates <- lengths(cand_rts)
  workable$status       <- status_v
  workable$row_id       <- seq_len(nrow(workable))

  # 6. STAGE 2 — cross-lab consensus in ALIGNED space.
  if (verbose) message("[anchor] Stage 2: aligned-space consensus per (mixture, compound)...")
  have_cpd <- all(c("mixture", "compound") %in% names(workable))
  consensus_aligned   <- new.env(parent = emptyenv())
  consensus_n_labs    <- new.env(parent = emptyenv())
  consensus_n_pol     <- new.env(parent = emptyenv())
  consensus_score     <- new.env(parent = emptyenv())
  if (have_cpd) {
    cand_rows <- which(status_v == "candidates_found")
    # Build ballots: one row per candidate. Forward-align to aligned space.
    ballots_lab    <- character(0)
    ballots_pol    <- character(0)
    ballots_mix    <- character(0)
    ballots_cpd    <- character(0)
    ballots_addct  <- character(0)
    ballots_rti    <- numeric(0)
    ballots_rta    <- numeric(0)
    ballots_cosmic <- numeric(0)
    ballots_mz     <- numeric(0)
    ballots_rowid  <- integer(0)
    for (ri in cand_rows) {
      cs <- cand_rts[[ri]]; if (!length(cs)) next
      al <- apply_fwd(cs, workable$lab[ri], workable$polarity[ri],
                      splines, anchor_table)
      ok <- !is.na(al)
      if (!any(ok)) next
      n <- sum(ok)
      ballots_lab    <- c(ballots_lab,    rep(workable$lab[ri], n))
      ballots_pol    <- c(ballots_pol,    rep(workable$polarity[ri], n))
      ballots_mix    <- c(ballots_mix,    rep(workable$mixture[ri], n))
      ballots_cpd    <- c(ballots_cpd,    rep(workable$compound[ri], n))
      ballots_addct  <- c(ballots_addct,
                          rep(workable$adduct[ri] %||% NA_character_, n))
      ballots_rti    <- c(ballots_rti,    cs[ok])
      ballots_rta    <- c(ballots_rta,    al[ok])
      cosmic_i <- if ("cosmic" %in% names(workable)) workable$cosmic[ri] else NA_real_
      ballots_cosmic <- c(ballots_cosmic, rep(cosmic_i, n))
      ballots_mz     <- c(ballots_mz,     rep(workable$ionMass[ri], n))
      ballots_rowid  <- c(ballots_rowid,  rep(ri, n))
    }
    ballots <- data.frame(
      lab        = ballots_lab,
      polarity   = ballots_pol,
      mixture    = ballots_mix,
      compound   = ballots_cpd,
      adduct     = ballots_addct,
      rt_cand    = ballots_rti,
      rt_aligned = ballots_rta,
      cosmic     = ballots_cosmic,
      ionMass    = ballots_mz,
      row_id     = ballots_rowid,
      stringsAsFactors = FALSE)

    # Similarity rule per (mixture, compound), in aligned space.
    if (nrow(ballots)) {
      ck <- paste(ballots$mixture, ballots$compound, sep = "\r")
      for (k in unique(ck)) {
        sub <- ballots[ck == k, , drop = FALSE]
        rts <- sub$rt_aligned
        lp  <- paste(sub$lab, sub$polarity, sep = "/")
        supp <- vapply(rts, function(r)
          length(unique(lp[abs(rts - r) <= agreement_window])), integer(1))
        if (!length(supp) || max(supp) == 0L) next
        max_supp_idx <- which(supp == max(supp))
        center <- stats::median(rts[max_supp_idx])
        in_w <- abs(rts - center) <= agreement_window
        supporting <- sub[in_w, , drop = FALSE]
        supporting <- supporting[!duplicated(supporting$row_id), , drop = FALSE]
        assign(k, stats::median(rts[in_w]), envir = consensus_aligned)
        assign(k, length(unique(supporting$lab)),
               envir = consensus_n_labs)
        assign(k, length(unique(supporting$polarity)),
               envir = consensus_n_pol)
        sc <- length(unique(supporting$lab)) * 10 +
              length(unique(supporting$polarity)) * 5 +
              sum(supporting$cosmic, na.rm = TRUE)
        assign(k, sc, envir = consensus_score)
      }
    }
  }

  # 7. STAGE 2.1 — cross-compound isobaric collision resolution.
  # Per mixture: pairs of compounds with consensus RTs within
  # isobar_rt_window AND any-adduct-pair m/z within ppm_tol → drop the
  # lower-scoring compound.
  isobar_dropped_set <- character(0)
  if (have_cpd) {
    # Compound → list of adduct m/z (per mixture).
    adduct_mz_per <- workable |>
      dplyr::filter(!is.na(ionMass), ionMass > 0,
                    status == "candidates_found") |>
      dplyr::distinct(mixture, compound, adduct, ionMass)
    by_mix <- split(adduct_mz_per, adduct_mz_per$mixture)

    keys_all <- ls(consensus_aligned)
    by_mix_k <- split(keys_all,
                       vapply(keys_all, function(k)
                         strsplit(k, "\r", fixed = TRUE)[[1]][1],
                         character(1)))
    for (mx in names(by_mix_k)) {
      ks <- by_mix_k[[mx]]
      if (length(ks) < 2L) next
      cpd_for_k <- vapply(ks, function(k)
        strsplit(k, "\r", fixed = TRUE)[[1]][2], character(1))
      rt_for_k  <- vapply(ks, function(k)
        get(k, envir = consensus_aligned), numeric(1))
      sc_for_k  <- vapply(ks, function(k)
        get0(k, envir = consensus_score, ifnotfound = 0), numeric(1))
      mz_map <- split(by_mix[[mx]]$ionMass, by_mix[[mx]]$compound)
      ord <- order(sc_for_k, decreasing = TRUE)
      keep <- rep(TRUE, length(ks))
      for (i in seq_len(length(ks) - 1)) {
        ii <- ord[i]; if (!keep[ii]) next
        if (is.na(rt_for_k[ii])) next
        mzi <- mz_map[[cpd_for_k[ii]]]; if (is.null(mzi)) next
        for (j in (i + 1):length(ks)) {
          jj <- ord[j]; if (!keep[jj]) next
          if (is.na(rt_for_k[jj])) next
          if (abs(rt_for_k[ii] - rt_for_k[jj]) > isobar_rt_window) next
          mzj <- mz_map[[cpd_for_k[jj]]]; if (is.null(mzj)) next
          ppm <- outer(mzi, mzj, function(a, b) abs(a - b) / a * 1e6)
          if (any(is.finite(ppm) & ppm <= ppm_tol)) {
            keep[jj] <- FALSE
            isobar_dropped_set <- c(isobar_dropped_set, ks[jj])
          }
        }
      }
    }
  }

  # 8. STAGE 2 SNAP — inverse-align consensus per row, snap to nearest
  #    candidate. Rows whose compound was isobar-dropped are dropped.
  inv_splines <- build_inverse_splines(splines, anchor_table)
  workable$rt_sec_anchored <- NA_real_
  workable$anchor_status   <- workable$status
  workable$rt_sec_original <- workable$rt_sec
  workable$consensus_rt_aligned <- NA_real_

  for (ri in seq_len(nrow(workable))) {
    cs <- cand_rts[[ri]]
    if (is.null(cs) || !length(cs)) next  # no_eic / no_peak / no_peak_passes_guards
    k <- paste(workable$mixture[ri], workable$compound[ri], sep = "\r")
    if (have_cpd && exists(k, envir = consensus_aligned, inherits = FALSE)) {
      if (k %in% isobar_dropped_set) {
        workable$anchor_status[ri] <- "isobar_dropped"
        next
      }
      cons_a <- get(k, envir = consensus_aligned)
      cons_l <- apply_inv(cons_a, workable$lab[ri], workable$polarity[ri],
                          inv_splines, anchor_table)
      if (!is.finite(cons_l)) {
        workable$rt_sec_anchored[ri] <- strongest_rt[ri]
        workable$anchor_status[ri]   <- "anchored_strongest_no_inv"
      } else {
        workable$rt_sec_anchored[ri] <- cs[which.min(abs(cs - cons_l))]
        workable$anchor_status[ri]   <- "anchored_consensus"
      }
      workable$consensus_rt_aligned[ri] <- cons_a
    } else {
      workable$rt_sec_anchored[ri] <- strongest_rt[ri]
      workable$anchor_status[ri]   <- "anchored_strongest"
    }
  }

  # 9. Re-assemble. Keep `consensus_rt_aligned` and `anchor_status` as
  # extra columns so downstream code (annotated_eic_pipeline.R, QMDs)
  # can read the aligned-space consensus directly without recomputing.
  anchored_keep <- c("anchored_consensus", "anchored_strongest",
                     "anchored_strongest_no_inv")
  anchored_out <- workable |>
    dplyr::mutate(rt_sec = ifelse(is.na(rt_sec_anchored), rt_sec,
                                   rt_sec_anchored)) |>
    dplyr::filter(anchor_status %in% anchored_keep) |>
    dplyr::select(dplyr::any_of(c(names(feat),
                                   "consensus_rt_aligned",
                                   "anchor_status")))

  # Per-lab diagnostics.
  per_lab_shifts <- workable |>
    dplyr::filter(anchor_status %in% anchored_keep) |>
    dplyr::mutate(rt_shift = rt_sec_anchored - rt_sec_original) |>
    dplyr::group_by(lab) |>
    dplyr::summarise(n = dplyr::n(),
                     median_shift_s = stats::median(rt_shift),
                     p25_shift_s    = stats::quantile(rt_shift, 0.25),
                     p75_shift_s    = stats::quantile(rt_shift, 0.75),
                     abs_p95_s      = stats::quantile(abs(rt_shift), 0.95),
                     .groups = "drop")

  per_lab_status <- workable |>
    dplyr::count(lab, anchor_status, name = "n") |>
    tidyr::pivot_wider(names_from = anchor_status, values_from = n,
                       values_fill = 0L)

  isobar_summary <- workable |>
    dplyr::filter(anchor_status == "isobar_dropped") |>
    dplyr::distinct(mixture, compound) |>
    dplyr::arrange(mixture, compound)

  diagnostics <- list(
    n_input        = nrow(in_scope),
    n_workable     = nrow(workable),
    n_kept         = nrow(anchored_out),
    n_dropped      = nrow(in_scope) - nrow(anchored_out),
    n_unique_eics  = nrow(unique_eics),
    n_isobar_dropped_rows = sum(workable$anchor_status == "isobar_dropped"),
    n_isobar_dropped_cpds = nrow(isobar_summary),
    per_lab_shifts = per_lab_shifts,
    per_lab_status = per_lab_status,
    isobar_summary = isobar_summary,
    overall_status = as.data.frame(
                       table(c(rep("no_path",  sum(bad_path)),
                               rep("no_mz",    sum(bad_mz & !bad_path)),
                               workable$anchor_status), useNA = "ifany"))
  )

  if (verbose) {
    message(sprintf("[anchor] kept %d / %d rows (dropped %d, isobar %d compounds)",
                    diagnostics$n_kept, diagnostics$n_input,
                    diagnostics$n_dropped, diagnostics$n_isobar_dropped_cpds))
    message("[anchor] per-lab RT shift summary:")
    print(per_lab_shifts)
  }

  result <- dplyr::bind_rows(pass_through, anchored_out)
  attr(result, "diagnostics") <- diagnostics
  result
}

