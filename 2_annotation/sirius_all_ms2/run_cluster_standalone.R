# SIRIUS all-MS2 annotation runner (cluster). Edit the CONFIG block, then:
#   export SIRIUS_PASSWORD='your_password'
#   Rscript run_sirius_cluster.R
# Each parallel job: give it its own `sirius_port` and a disjoint `mixtures`.
#
# ===================== FEATURE-DEFINITION STRATEGY =====================
# The recovery of a spiked-in panel compound depends only on how we turn raw
# MS2 scans into the ONE spectrum-per-feature that SIRIUS annotates. Two steps:
#
# (1) GROUPING (identical for every lab): precursor-tolerant MERGE.
#     Group MS2 by precursor m/z within 20 ppm, then split each m/z group
#     wherever the RT gap exceeds 15 s (single-linkage) -> one feature per real
#     precursor peak. This replaces the old exact-precursor "rle" split. Paired
#     with the richest-scan step below it loses no compound on ANY lab (incl.
#     icl/Waters) and yields far fewer features than rle (icl 946->205, afekta
#     2715->1548, hmgu ~1094, cembio 194->183) => much less SIRIUS compute.
#
# (2) MS2 -> spectrum(s) per feature: RICHEST-SCAN representative (all labs except
#     afekta -- see afekta note below). We do NOT import every scan and let SIRIUS
#     build a consensus. We keep the SINGLE scan with the most fragment peaks and set
#     its precursor to the feature's MEDIAN precursor m/z. Why this beats SIRIUS's own
#     consensus merge: a feature's scans vary in quality (weak scans, over-fragmented
#     CE-ladder rungs, isobaric contamination) and in precursor accuracy (Waters
#     drift). A blind consensus dilutes the best scan and can inherit a drifted
#     precursor. One clean, peak-rich scan + an accurate median precursor recovered
#     MORE compounds (icl UMP, cembio glucosyl-mannose). Confirmed panel-wide: an
#     all-scans consensus was net -5 for afekta (regresses N-heterocycles/nucleotides).
#
# Per-lab differences are only in the PRE-STEPS, not the grouping/representative:
#   * icl  (Waters): mzML precursor m/z is unreliable/empty -> re-derive it with
#            estimatePrecursorMz() before anything else (without it recovery=0). EXCEPTION to
#            step (2) -- uses ALL-SCANS (like cembio): richest picks the noisiest co-isolated scan
#            of over-split features. Modest win (+1..+3/mix, beats semi-auto everywhere, a couple
#            regressions) -> richest UNION all-scans is the true ceiling. Branch = reduce_allscans.
#   * hmgu (DDA, noisy): a precursor-abundance filter (prec_thr=300) drops weak-
#            precursor DDA-noise features. Intensity is the MAX precursor peak
#            over the feature's scans, NOT one median-RT scan (that single-scan
#            estimate sampled off-apex dips and wrongly dropped real compounds).
#   * cembio: MS2 injection is MS2-only (CE ladder 10/20/40 eV), no prec_thr. EXCEPTION to
#            step (2) -- uses ALL-SCANS (keep every MS2 scan/feature, capped at 40 richest,
#            precursor=median; SIRIUS builds its own consensus). richest-scan collapses to
#            2-3/20 on the group-2/3 mixtures; all-scans rebuilds it (+53 over all 19).
#            richest-per-CE was decisively WORSE (-48) and is rejected; do not retry it.
#            Branch = `reduce_allscans`; batch_size 30.
#            MS1 (adopted 2026-07-18): cembio DOES get an isotope pattern, taken from its
#            SEPARATE MS1 INJECTION rather than from the MS2 injection (whose own MS1 holds
#            the precursor only ~74% of the time vs ~100% for the MS1 injection). See the
#            ms1_injection block below for the pairing and the three safety rules. Worth
#            +13 (185 -> 198), which beats semi-auto's 193 for the first time.
#   * afekta: clean, prec_thr=0, MS1 attached. EXCEPTION to step (2) -- afekta
#            acquires a 10/20/40 eV CE ladder, so it uses RICHEST-PER-CE (richest scan
#            at EACH collision energy, ~3 energy-resolved spectra/feature) + a cleaner
#            prep (5%-only filter, no cap, >=2 peaks). Full-panel over 16 mixtures:
#            net +6 vs single richest (+18 small-acid gains fumaric/taurine/glutaric/
#            mandelic/malonic..., -12). Branch = `is_afekta` / `reduce_richest_perCE`.
# mix-1_1 recovery: hmgu 17, cembio 14 (top-1; a GOOD mixture, method-insensitive), icl 6; afekta richest-per-CE
# net +6 panel-wide (vs 14 single-richest). cembio's all-scans win shows only on weak group-2/3 mixtures (see above).
# =======================================================================
#
## ===================== CONFIG =====================
sirius_bin    <- "/home/plouail/tools/sirius/bin"  # dir containing the sirius executable
DATA_DIR      <- "~/ringtrial_data"   # ONE folder holding <lab>/<class>/{mzml, seq_pos.xlsx, seq_neg.xlsx}; results go to <DATA_DIR>/results/
sirius_port   <- 8080                 # base port, UNIQUE per concurrent job. Prefer >= 1024 (ports
                                                # below 1024 are privileged on Linux and may fail to bind). If
                                                # the service won't start, the port is often held by a STALE
                                                # daemon from a prior failed attempt (orphaned java) -> kill it
                                                # or pick a fresh port. Space concurrent jobs >= 20 apart; the
                                                # per-batch retry auto-bumps the port on failure.
# ---------------------------------------------------------------------------------------------
# CURRENT JOB (2026-07-20): regenerate the 12 hmgu PRIMARY mixtures whose sirius_summary_*.csv
# was lost. sirius_full exists for all 19 hmgu mixtures, but the cache test needs BOTH files, so
# these 12 re-run and write summary + full + extraction_counts fresh.
# Paste back ALL THREE files per mixture (summary, full, extraction_counts): re-running also
# regenerates `full`, and SIRIUS is not perfectly deterministic, so keeping the NEW summary next
# to the OLD full would leave that mixture internally inconsistent.
# The other three labs are already complete at 19/19 summary + full and must NOT be re-run.
mixtures      <- NULL          # NULL = auto-discover every mixture in the class; or set c("4_1","4_2",...) to subset
labs          <- c("icl")      # EDIT: ONE lab per job (afekta / hmgu / icl). Not cembio.
# ---------------------------------------------------------------------------------------------
representative <- "primary"   # "primary" = the finalized best method per lab (afekta richest-per-CE,
                              #             cembio + icl ALL-SCANS, hmgu richest-scan)  -> results/primary
                              #             cembio additionally attaches MS1 from its separate MS1
                              #             injection (see the ms1_injection block below).
                              # "richest"  = force SINGLE-RICHEST for every lab (+ the production prep, so it
                              #             matches the historical baseline) -> results/richest
                              # Union only ACROSS representatives (richest vs per-CE/all-scans), never a
                              # no-MS1 pass with its +MS1 counterpart: MS1 is independent evidence, so an
                              # annotation the isotope pattern contradicts should not be kept anyway.
provenance_only <- FALSE      # TRUE = ONLY extract + write provenance_<lab>_mix_<m>.csv (which scans of which
                              # file were fed to SIRIUS for each feature), then SKIP SIRIUS entirely. The
                              # extraction is deterministic, so this regenerates provenance for results that
                              # already exist WITHOUT re-running any annotation. Needed for spectral-library
                              # building: sirius_full gives feature -> structure, provenance gives feature -> scans.
study_group   <- "FE"          # EDIT: "FE" (Food Exposome) or "Others" (Chemical+Eukaryotic+Plant)
sirius_user   <- "plouail@eurac.edu"
sirius_pass   <- Sys.getenv("SIRIUS_PASSWORD")   # run: export SIRIUS_PASSWORD=your_pw
hmgu_prec_thr <- 300                               # hmgu precursor-intensity filter (others = 0/off)
mz_range      <- c(30, 1050)                        # precursor-m/z window applied to every lab. Spans the
                                                    # full mass range ever RECORDED in the semi-auto curated
                                                    # annotation (formamide ~44 to Glycerol Trieicosanoate
                                                    # ~1020 m/z) with margin; drops only higher-mass matrix
                                                    # and makes labs cross-comparable (cf. lab_setup_phase1.csv).
sirius_batch_size <- 100                            # features per SIRIUS batch. Each batch is cached to CSV +
                                                    # retried (fresh daemon) on failure, so a crash / internet
                                                    # blip resumes mid-lab from the first un-cached batch, not
                                                    # from scratch. Set larger to cut per-batch overhead if the
                                                    # connection is stable.
sirius_xmx_gb <- 32                                 # HARD CAP on the SIRIUS JVM heap. Unset, the JVM sizes its
                                                    # heap from the NODE's RAM (not the SLURM --mem cgroup) and
                                                    # grows past your limit -> SLURM OOM-kills it. Keep it well
                                                    # under --mem: leave room for R (~8 GB) + JVM off-heap
                                                    # (~heap*0.2) + margin. Rule of thumb: 2*xmx + 16 < --mem
                                                    # (32 -> needs ~80 GB, safe at --mem=100G). Small batches
                                                    # need far less than 32 GB, so lower it if --mem is tighter.
## ==================================================

if (!nzchar(sirius_pass)) stop("SIRIUS_PASSWORD env var not set.")
Sys.setenv(PATH = paste(sirius_bin, Sys.getenv("PATH"), sep = ":"))
# Bound the SIRIUS service JVM heap so it respects the SLURM --mem limit (see
# sirius_xmx_gb). Set BEFORE any Sirius() launch so the service JVM inherits it.
# NOTE: if SIRIUS's own launcher overrides -Xmx, verify with `ps -o args -C java`
# on the node and set SIRIUS's native memory config instead.
Sys.setenv("_JAVA_OPTIONS" = sprintf("-Xmx%dg", sirius_xmx_gb))

suppressPackageStartupMessages({
  library(Spectra); library(MsExperiment); library(xcms); library(RuSirius)
  library(dplyr); library(tidyr); library(readxl)
})
## standalone: no repo setup.R needed (readxl loaded above)

stopifnot(representative %in% c("primary", "richest"))
output_dir <- file.path(DATA_DIR, "results", representative)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## auto-discover mixtures from the class's seq files when not set explicitly
if (is.null(mixtures)) {
  .sq <- rbind(as.data.frame(read_xlsx(file.path(DATA_DIR, labs[1], study_group, "seq_pos.xlsx"))),
               as.data.frame(read_xlsx(file.path(DATA_DIR, labs[1], study_group, "seq_neg.xlsx"))))
  .mm <- gsub("[.]", "_", sub(".*_", "", sub("_MS2$", "", .sq[["Data File"]])))
  mixtures <- sort(unique(.mm[grepl("^[0-9]+_[0-9]+$", .mm)]))
  message("auto-discovered ", length(mixtures), " mixtures: ", paste(mixtures, collapse = ", "))
}

prec_thr_by_lab <- c(hmgu = hmgu_prec_thr, afekta = 0, icl = 0, cembio = 0)

extract_all_ms2 <- function(lab_name, study_group, mixture) {
  study_path <- file.path(DATA_DIR, lab_name, study_group)
  mzml_dir   <- file.path(study_path, "mzml")
  mix_label  <- mixture

  empty_result <- function(sp = Spectra(backend = MsBackendMemory())) {
    list(spectra_pos = sp, spectra_neg = sp,
         n_ms2_pos = 0L, n_ms2_neg = 0L,
         n_features_pos = 0L, n_features_neg = 0L)
  }

  seq_pos <- read_xlsx(file.path(study_path, "seq_pos.xlsx")) |> as.data.frame()
  seq_neg <- read_xlsx(file.path(study_path, "seq_neg.xlsx")) |> as.data.frame()
  seq_pos$polarity <- "pos"; seq_neg$polarity <- "neg"
  seq <- rbind(seq_pos, seq_neg)
  seq$filename <- paste0(seq$`Data File`, ".mzML")
  seq <- seq[seq$filename %in% list.files(mzml_dir), ]
  seq$mixture <- gsub("\\.", "_",
                      sub("target$", "",
                          sub(".*_", "",
                              sub("_MS2$", "", seq$`Data File`))))

  has_ms2_files <- any(grepl("MS2", seq$filename, fixed = TRUE))
  seq <- if (has_ms2_files) {
    seq[grepl("MS2", seq$filename, fixed = TRUE) & seq$mixture == mix_label, ]
  } else {
    seq[seq$mixture == mix_label, ]
  }
  if (nrow(seq) == 0L) return(empty_result())

  sp <- Spectra(file.path(mzml_dir, seq$filename), backend = MsBackendMzR())
  sp <- setBackend(sp, MsBackendMemory())

  # Waters (icl) PRE-STEP: the mzML records an unreliable / empty precursor m/z,
  # so re-derive it from the co-eluting MS1 apex. Without this every icl feature
  # gets a wrong precursor and recovery collapses to 0.
  if (lab_name == "icl") sp$precursorMz <- estimatePrecursorMz(sp)

  # afekta uses the RICHEST-PER-CE representative (see reduce fn below), validated
  # with a cleaner spectrum prep: 5%-of-base filter ONLY (no >=1000 floor), NO
  # top-60 cap, keep >=2 peaks. Every other lab keeps the production prep
  # (5%|1000 filter, cap60, >=3 peaks + single richest scan).
  # In the "richest" (union-baseline) pass EVERY lab uses the production prep + single-richest, so
  # afekta's cleaner prep is only active for the primary pass -> the richest side stays comparable
  # to the historical baseline.
  is_afekta <- lab_name == "afekta" && representative == "primary"

  # MS2 peak cleanup: keep peaks >=5% base peak (afekta) or >=5% base | >=1000 abs.
  sp <- filterIntensity(sp, intensity = function(x) {
    if (length(x) == 0L) return(logical(0))
    bp <- max(x, na.rm = TRUE)
    if (!is.finite(bp) || bp <= 0) return(rep(TRUE, length(x)))
    if (is_afekta) x > bp * 0.05 else x >= bp * 0.05 | x >= 1000
  }, msLevel. = 2L)
  sp <- applyProcessing(sp)
  sp <- sp[sp$msLevel != 2L | lengths(sp) >= (if (is_afekta) 2L else 3L)]

  # Cap each MS2 to its top-60 most intense peaks. SIRIUS tree/fingerprint cost
  # scales super-linearly with peak count, and neg-mode noise spectra carry up
  # to ~250 peaks (panel compounds have <=19), so a few dominate a whole batch.
  # Validated: every should-identify panel compound is below 60 peaks, so the
  # cap leaves their input unchanged. Also defines "richest" below by peak count.
  cap_top <- function(x, n = 60L) {
    if (nrow(x) <= n) return(x)
    x[sort(order(x[, "intensity"], decreasing = TRUE)[seq_len(n)]), , drop = FALSE]
  }

  # (1) GROUPING (all labs): precursor-tolerant merge, 20 ppm + 15 s single-link.
  group_precursor_rt <- function(file, rt, pmz, ppm = 20, rt_gap = 15) {
    g <- integer(length(rt)); ng <- 0L
    for (f in unique(file)) {
      idx <- which(file == f); o <- idx[order(pmz[idx])]; ps <- pmz[o]
      mzg <- cumsum(c(TRUE, (diff(ps) / ps[-length(ps)] * 1e6) > ppm))
      for (m in unique(mzg)) {
        ii <- o[mzg == m]; ii <- ii[order(rt[ii])]
        sub <- cumsum(c(TRUE, diff(rt[ii]) > rt_gap))
        for (s in unique(sub)) { ng <- ng + 1L; g[ii[sub == s]] <- ng }
      }
    }
    g
  }

  # (2a) RICHEST-SCAN representative (all labs except afekta): one MS2 per feature =
  # the scan with the most peaks; precursor := feature MEDIAN precursor m/z.
  reduce_richest <- function(ms2) {
    grp_f <- factor(ms2$group_idx)
    npk   <- lengths(ms2)
    g_med <- tapply(ms2$precursorMz, grp_f, median)
    sel   <- unlist(tapply(seq_along(ms2), grp_f,
                           function(ii) ii[which.max(npk[ii])]))
    keep  <- ms2[sel]
    keep$precursorMz <- as.numeric(g_med[as.character(keep$group_idx)])
    keep
  }

  # (2b) RICHEST-PER-CE representative (afekta): afekta acquires a collision-energy
  # ladder (10/20/40 eV). Keep the richest scan at EACH collision energy (~3 energy-
  # resolved spectra/feature) so SIRIUS builds an energy-resolved tree/fingerprint;
  # precursor := feature MEDIAN m/z. Full-panel head-to-head over 16 afekta mixtures:
  # net +6 vs single richest (+18 small-acid gains, -12), while an all-scans consensus
  # was net -5. If collisionEnergy is missing/single -> one richest scan (= 2a).
  reduce_richest_perCE <- function(ms2) {
    ce  <- suppressWarnings(round(as.numeric(ms2$collisionEnergy)))
    ce[is.na(ce)] <- -1L
    fce <- interaction(factor(ms2$group_idx), factor(ce), drop = TRUE)
    npk <- lengths(ms2)
    sel <- sort(unlist(tapply(seq_along(ms2), fce,
                              function(ii) ii[which.max(npk[ii])]), use.names = FALSE))
    keep  <- ms2[sel]
    g_med <- tapply(ms2$precursorMz, factor(ms2$group_idx), median)
    keep$precursorMz <- as.numeric(g_med[as.character(keep$group_idx)])
    keep
  }
  # (2c) ALL-SCANS representative (cembio): keep EVERY MS2 scan per feature (capped at the 40 richest so
  # scan-rich noise features don't choke SIRIUS), precursor := feature MEDIAN m/z; SIRIUS builds its own
  # consensus across them. On cembio's MS2-only CE-ladder data the single richest scan is a poor
  # representative on the group-2/3 mixtures (recovery collapses to 2-3/20); all-scans rebuilds it: over 6
  # mixtures top-1 recovery 63 vs richest 39 (+24), ~= semi-auto 67. richest-per-CE was WORSE (36).
  reduce_allscans <- function(ms2, cap = 40L) {
    grp_f <- factor(ms2$group_idx)
    g_med <- tapply(ms2$precursorMz, grp_f, median)
    npk   <- lengths(ms2)
    sel   <- unlist(tapply(seq_along(ms2), grp_f, function(ii)
                    if (length(ii) <= cap) ii else ii[order(npk[ii], decreasing = TRUE)[seq_len(cap)]]),
                    use.names = FALSE)
    keep  <- ms2[sort(sel)]
    keep$precursorMz <- as.numeric(g_med[as.character(keep$group_idx)])
    keep
  }
  # cembio + icl use ALL-SCANS (both suffer the noisy-richest problem: cembio's MS2-only weak mixtures,
  # icl's over-split Waters DDA). icl is a MODEST win (+1..+3/mix, a couple regressions -> union with richest
  # is the true ceiling), cembio a large one (+24). afekta uses richest-per-CE; hmgu single-scan -> richest.
  reduce_fun <- if (representative == "richest") reduce_richest                       # union baseline: all labs
                else if (lab_name %in% c("cembio", "icl")) reduce_allscans
                else if (is_afekta) reduce_richest_perCE
                else reduce_richest

  ms2 <- sp[sp$msLevel == 2L]
  ms2 <- filterEmptySpectra(ms2)
  if (length(ms2) == 0L) return(empty_result())
  if (!is_afekta) ms2 <- applyProcessing(addProcessing(ms2, cap_top, n = 60L))

  # Uniform precursor-m/z window (mz_range = 30-1050): keep only features within
  # the range of masses ever recorded in the semi-auto annotation; drops the
  # high-mass matrix tail and makes labs cross-comparable.
  ms2 <- ms2[!is.na(ms2$precursorMz) &
             ms2$precursorMz >= mz_range[1] & ms2$precursorMz <= mz_range[2]]
  if (length(ms2) == 0L) return(empty_result())

  ms2$group_idx <- group_precursor_rt(as.character(ms2$dataOrigin),
                                      ms2$rtime, ms2$precursorMz)

  # Stash the RAW per-scan precursor m/z BEFORE the reduce step below overwrites precursorMz with
  # the feature median. This is what the spectral library needs: one true experimental precursor
  # per contributing scan, not the median repeated. It rides along as a spectra variable through
  # the subset that reduce_fun keeps. (For icl this is the estimatePrecursorMz value, since the
  # Waters mzML precursor is empty; for the other labs it is the mzML precursor of that scan.)
  ms2$precursor_mz_raw <- ms2$precursorMz

  # MS1 for the isotope pattern (and hmgu's precursor-abundance filter).
  #
  # CEMBIO (adopted 2026-07-18): cembio was previously run with NO MS1, because the MS1 riding
  # inside its MS2 injection holds the feature precursor only ~74% of the time and attaching a
  # scan that lacks the precursor gives SIRIUS a meaningless isotope pattern. But cembio also
  # acquires a SEPARATE MS1 INJECTION which holds it ~100% of the time. Pair the two by polarity
  # + run number (NOT by filename: the 2_1 MS2 file is "...2.1target.mzML" while its MS1 file is
  # "...2.1.mzML"). They run back-to-back so RT drift is negligible (median ~1 s) and no alignment
  # model is needed. Measured over all 19 mixtures, top-1 panel recovery:
  #   richest 132 -> 138 (+6);  all-scans 185 -> 198 (+13);  semi-auto is 193, so all-scans+MS1
  #   beats the curated workflow for the first time. Union of the two MS1 passes = 206.
  ms1_injection <- lab_name == "cembio"
  ms1 <- sp[sp$msLevel == 1L]
  ms2_to_ms1_file <- NULL
  if (ms1_injection) {
    inj <- list.files(mzml_dir, full.names = TRUE)
    inj <- inj[!grepl("MS2", basename(inj), fixed = TRUE)]
    runno  <- function(b) { m <- regmatches(b, regexpr("_[0-9]{3}_", b)); if (length(m)) m else NA_character_ }
    pol_of <- function(b) if (grepl("_pos_", b, fixed = TRUE)) "_pos_" else "_neg_"
    ms2_files <- unique(as.character(ms2$dataOrigin))
    map <- setNames(rep(NA_character_, length(ms2_files)), ms2_files)
    for (f in ms2_files) {
      b <- basename(f); rn <- runno(b); if (is.na(rn)) next
      hit <- inj[grepl(rn, basename(inj), fixed = TRUE) & grepl(pol_of(b), basename(inj), fixed = TRUE)]
      if (length(hit)) map[[f]] <- hit[1]
    }
    got <- sum(!is.na(map))
    message("  MS1 injection: paired ", got, "/", length(ms2_files), " MS2 files")
    if (got > 0) {
      ms1 <- setBackend(Spectra(unique(na.omit(unname(map))), backend = MsBackendMzR()), MsBackendMemory())
      ms1 <- ms1[ms1$msLevel == 1L]
      ms2_to_ms1_file <- map
    } else ms1_injection <- FALSE
  }
  use_ms1 <- (lab_name != "cembio" || ms1_injection) && length(ms1) > 0L

  if (use_ms1) {
    grp_f  <- factor(ms2$group_idx)
    g_ids  <- as.integer(levels(grp_f))
    g_file <- tapply(as.character(ms2$dataOrigin), grp_f, `[`, 1L)
    g_rt   <- tapply(ms2$rtime, grp_f, median)
    g_mz   <- tapply(ms2$precursorMz, grp_f, median)
    if (ms1_injection) {
      # WINDOWS TRAP: Spectra normalises dataOrigin to BACKSLASHES while list.files() returns
      # forward slashes, so a full-path key silently matches nothing (0 features attached, which
      # reads like a negative result rather than a bug). Key by basename on both sides.
      by_ms1file <- split(seq_along(ms1), basename(as.character(ms1$dataOrigin)))
      ms1_by_file <- setNames(lapply(names(ms2_to_ms1_file), function(f) {
        t <- ms2_to_ms1_file[[f]]
        if (is.na(t)) return(integer(0))
        v <- by_ms1file[[basename(t)]]; if (is.null(v)) integer(0) else v }),
        names(ms2_to_ms1_file))
    } else {
      ms1_by_file <- split(seq_along(ms1), as.character(ms1$dataOrigin))
    }
    # Take the precursor APEX over the feature's OWN rt span, not the scan nearest its MEDIAN rt:
    # grouping produces features spanning 85-125 s, so the median is a meaningless midpoint and
    # picks a scan from the wrong part of the chromatogram. Then GATE on isotope quality: attach
    # only if M+0 clears the noise floor AND M+1 is present at a plausible ratio, else attach
    # nothing. A weak precursor with no M+1 tells SIRIUS "this molecule has no C-13", which
    # penalises carbon-rich formulas and DEMOTES the correct one. No MS1 beats a false MS1:
    # the gate recovered 5-oxo-L-proline and quinolinic acid while keeping every gain (198 vs 196).
    g_rtmin <- tapply(ms2$rtime, grp_f, min); g_rtmax <- tapply(ms2$rtime, grp_f, max)
    ms1_pd_all <- if (ms1_injection) peaksData(ms1) else NULL
    RT_PAD <- 10; PPM <- 20; MS1_MIN_INT <- 1000
    pick <- vapply(seq_along(g_ids), function(i) {
      cand <- ms1_by_file[[ g_file[i] ]]
      if (is.null(cand) || !length(cand)) return(NA_integer_)
      if (!ms1_injection) return(cand[which.min(abs(ms1$rtime[cand] - g_rt[i]))])
      rt <- ms1$rtime[cand]
      near <- cand[rt >= g_rtmin[i] - RT_PAD & rt <= g_rtmax[i] + RT_PAD]
      if (!length(near)) return(NA_integer_)
      mz <- g_mz[i]; tol <- mz * PPM / 1e6
      ints <- vapply(near, function(k) {
        p <- ms1_pd_all[[k]]; if (!nrow(p)) return(0)
        h <- abs(p[, "mz"] - mz) <= tol
        if (any(h)) max(p[h, "intensity"]) else 0
      }, numeric(1))
      if (all(ints <= 0)) return(NA_integer_)
      k <- near[which.max(ints)]
      p <- ms1_pd_all[[k]]
      h0 <- abs(p[, "mz"] - mz) <= tol
      if (!any(h0)) return(NA_integer_)
      i0 <- max(p[h0, "intensity"]); if (i0 < MS1_MIN_INT) return(NA_integer_)
      h1 <- abs(p[, "mz"] - (mz + 1.00336)) <= tol * 2
      if (!any(h1)) return(NA_integer_)
      r1 <- max(p[h1, "intensity"]) / i0
      if (r1 < 0.005 || r1 > 1) return(NA_integer_)
      k
    }, integer(1))
    if (ms1_injection)
      message("  MS1 attached to ", sum(!is.na(pick)), "/", length(g_ids), " features")
    keep_g  <- !is.na(pick)
    ms1_sel <- ms1[pick[keep_g]]
    ms1_sel$group_idx <- g_ids[keep_g]

    # Trim the attached MS1 to the precursor's ISOTOPE WINDOW. A full survey scan is ~600 peaks and
    # SIRIUS then fits the wrong envelope in crowded low-m/z regions. Use addProcessing on the
    # EXISTING Spectra: hand-building one breaks the RuSirius import ("object of type 'environment'
    # is not subsettable"), see reference_rusirius_import_spectra.
    if (ms1_injection && length(ms1_sel)) {
      ms1_sel$prec_mz <- as.numeric(g_mz[keep_g])
      # window: M-1.5 to M+7.5 Da. The +7.5 upper bound (widened from +5.5 on 2026-07-20) covers
      # the full isotope envelope SIRIUS can use for this panel, including the M+6 peak of the few
      # tri-chlorinated compounds; below +7.5 the panel is CHNOPS + <=1 Cl, all inside +5.5.
      iso_trim <- function(x, prec_mz, ...) {
        if (!nrow(x) || is.na(prec_mz[1])) return(x)
        x[x[, "mz"] >= prec_mz[1] - 1.5 & x[, "mz"] <= prec_mz[1] + 7.5, , drop = FALSE]
      }
      ms1_sel <- applyProcessing(addProcessing(ms1_sel, iso_trim, spectraVariables = "prec_mz"))
    }

    # Precursor-abundance filter (hmgu only): drop weak-precursor DDA-noise
    # features. Kept intensity = MAX precursor peak over the feature's SCANS
    # (each scan vs its own nearest MS1), NOT one median-RT scan. The single-scan
    # estimate sampled off-apex dips and wrongly dropped real multi-scan
    # compounds (e.g. O-phospho-L-tyrosine: median-RT read 253, group max 4719).
    prec_thr <- prec_thr_by_lab[[lab_name]]
    if (length(prec_thr) && !is.na(prec_thr) && prec_thr > 0) {
      ms1_pd <- peaksData(ms1)
      fo <- as.character(ms2$dataOrigin)
      scan_int <- vapply(seq_along(ms2), function(k) {
        cand <- ms1_by_file[[ fo[k] ]]
        if (is.null(cand) || !length(cand)) return(NA_real_)
        s <- cand[which.min(abs(ms1$rtime[cand] - ms2$rtime[k]))]
        p <- ms1_pd[[s]]
        w <- p[abs(p[, "mz"] - ms2$precursorMz[k]) <= 0.01, "intensity"]
        if (length(w)) max(w) else 0
      }, numeric(1))
      g_max <- tapply(scan_int, grp_f,
                      function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
      keep_ids <- as.integer(names(g_max))[is.na(g_max) | g_max >= prec_thr]
      ms1_sel  <- ms1_sel[ms1_sel$group_idx %in% keep_ids]
      ms2      <- ms2[ms2$group_idx %in% keep_ids]
      if (length(ms2) == 0L) return(empty_result())
    }

    ms2     <- reduce_fun(ms2)                           # richest scan (or richest-per-CE for afekta)
    ms1_sel <- ms1_sel[ms1_sel$group_idx %in% unique(ms2$group_idx)]
    sp      <- c(ms1_sel, ms2)
  } else {
    sp <- reduce_fun(ms2)                                # no MS1 available for this feature set
  }
  if (length(sp) == 0L) return(empty_result(sp))

  list(
    spectra_pos    = sp[sp$polarity == 1L],
    spectra_neg    = sp[sp$polarity == 0L],
    n_ms2_pos      = sum(sp$msLevel == 2L & sp$polarity == 1L),
    n_ms2_neg      = sum(sp$msLevel == 2L & sp$polarity == 0L),
    n_features_pos = length(unique(sp$group_idx[sp$polarity == 1L])),
    n_features_neg = length(unique(sp$group_idx[sp$polarity == 0L]))
  )
}

sirius_sum_path       <- function(lab, m) file.path(output_dir, sprintf("sirius_summary_%s_mix_%s.csv",   lab, m))
sirius_full_path      <- function(lab, m) file.path(output_dir, sprintf("sirius_full_%s_mix_%s.csv",      lab, m))
extraction_count_path <- function(lab, m) file.path(output_dir, sprintf("extraction_counts_%s_mix_%s.csv", lab, m))
provenance_path       <- function(lab, m) file.path(output_dir, sprintf("provenance_%s_mix_%s.csv",        lab, m))

# PROVENANCE: the exact MS2 scans handed to SIRIUS for each feature. sirius_full_* gives
# feature -> structure; this gives feature -> (file, scan, RT, precursor, CE). Join them on
# feature_id == xcms_fts to trace an annotation back to a real spectrum (spectral-library building).
# richest -> 1 row/feature; per-CE -> ~3 (one per collision energy); all-scans -> up to 40.
spectra_prov <- function(s) {
  if (length(s) == 0L) return(NULL)
  m <- s[s$msLevel == 2L]
  if (length(m) == 0L) return(NULL)
  sv  <- spectraVariables(m)
  get <- function(nm, fallback) if (nm %in% sv) m[[nm]] else fallback
  data.frame(
    feature_id       = as.integer(m$group_idx),
    file             = basename(as.character(m$dataOrigin)),
    scan_index       = as.integer(get("scanIndex", NA_integer_)),
    acquisition_num  = as.integer(get("acquisitionNum", NA_integer_)),
    rt_sec           = round(as.numeric(m$rtime), 2),
    precursor_mz_raw = as.numeric(get("precursor_mz_raw", m$precursorMz)),  # true per-scan precursor
    precursor_mz_med = as.numeric(m$precursorMz),                           # feature median (what SIRIUS saw)
    collision_energy = suppressWarnings(as.numeric(get("collisionEnergy", NA_real_))),
    n_peaks          = lengths(m),
    polarity         = ifelse(m$polarity == 1L, "pos", "neg"),
    stringsAsFactors = FALSE)
}

sirius_id_cols <- c("alignedFeatureId", "compoundId", "externalFeatureId",
                    "formulaId", "xcms_fts", "sirius_fts",
                    "molecularFormulaId", "structureCandidateInChIKey")

run_sirius_for_lab <- function(sp_pos, sp_neg,
                               username, password, projectId, path,
                               base_port, batch_size = 100L, max_try = 3L) {
  if (length(sp_pos) == 0L && length(sp_neg) == 0L)
    return(list(summary = NULL, full = NULL, complete = TRUE))

  groups_pos <- unique(sp_pos$group_idx)
  groups_neg <- unique(sp_neg$group_idx)
  pol_lookup <- c(
    setNames(rep("pos", length(groups_pos)), as.character(groups_pos)),
    setNames(rep("neg", length(groups_neg)), as.character(groups_neg))
  )

  # Split the mixture into interleaved batches: sort features by precursor mass,
  # then round-robin into `nb` batches so no batch is an all-high-mass block
  # (a solid high-mass block makes the SIRIUS daemon unresponsive). Each batch's
  # summary + full results are cached to CSV; a crash / internet blip resumes
  # from the first un-cached batch instead of redoing the whole lab-mixture.
  gm  <- tapply(c(sp_pos$precursorMz[sp_pos$msLevel == 2L], sp_neg$precursorMz[sp_neg$msLevel == 2L]),
                factor(c(sp_pos$group_idx[sp_pos$msLevel == 2L], sp_neg$group_idx[sp_neg$msLevel == 2L])),
                median)
  allg <- as.integer(names(gm)); o <- order(as.numeric(gm))
  nb  <- max(1L, ceiling(length(allg) / batch_size))
  bid <- integer(length(allg)); bid[o] <- ((seq_along(o) - 1L) %% nb) + 1L
  batches <- split(allg, bid)

  cache_dir <- file.path(path, "_batch_cache"); dir.create(cache_dir, showWarnings = FALSE)
  sum_cf  <- function(i) file.path(cache_dir, sprintf("%s_bsum_%03d.csv",  projectId, i))
  full_cf <- function(i) file.path(cache_dir, sprintf("%s_bfull_%03d.csv", projectId, i))

  fid <- formulaIdParam(
    numberOfCandidates = 5, instrument = "QTOF", numberOfCandidatesPerIonization = 2,
    massAccuracyMS2ppm = 10, filterByIsotopePattern = FALSE, isotopeMs2Settings = "SCORE",
    # de novo formula generation is combinatorial; off => fast bottom-up BIO search
    # (still finds the spike-in panel). ilpTimeout bounds the ILP tree step, and the
    # per-instance/decomposition seconds cap the formula enumeration on high-mass
    # features (hmgu) - the runaway that otherwise balloons memory + wall-time.
    performDeNovoBelowMz = 0, formulaSearchDBs = c("BIO"), ilpTimeout = TRUE,
    numberOfSecondsPerInstance = 60, numberOfSecondsPerDecomposition = 20)

  # One daemon, reused across batches (delete-replace); a fresh one (next port) is
  # spun up only after a failure. SIRIUS's workspace is single-lock, so we never
  # run two daemons for the same project at once.
  srs <- NULL; nxt <- base_port
  newsess <- function() { p <- nxt; nxt <<- nxt + 1L
    Sirius(username = username, password = password, projectId = projectId, path = path, port = p) }
  on.exit(try(shutdown(srs), silent = TRUE), add = TRUE)

  for (i in seq_along(batches)) {
    if (file.exists(sum_cf(i)) && file.exists(full_cf(i))) {
      message("    batch ", i, "/", length(batches), " cached"); next
    }
    bp <- sp_pos[sp_pos$group_idx %in% batches[[i]]]
    bn <- sp_neg[sp_neg$group_idx %in% batches[[i]]]
    np <- length(unique(bp$group_idx)); nn <- length(unique(bn$group_idx))
    exp_ids <- as.character(batches[[i]]); done <- FALSE
    for (att in 1:max_try) {
      ok <- tryCatch({
        if (is.null(srs)) srs <- newsess()
        if (np > 0) srs <- import(srs, spectra = bp, ms_column_name = "group_idx", deleteExistingFeatures = TRUE)
        if (nn > 0) srs <- import(srs, spectra = bn, ms_column_name = "group_idx", deleteExistingFeatures = np == 0)
        # re-import any silently dropped feature (SIRIUS occasionally skips one)
        miss <- setdiff(exp_ids, as.character(featuresId(srs, type = "xcms")))
        for (r in 1:2) {
          if (!length(miss)) break
          mp <- sp_pos[sp_pos$group_idx %in% miss]; mn <- sp_neg[sp_neg$group_idx %in% miss]
          if (length(mp)) srs <- import(srs, spectra = mp, ms_column_name = "group_idx", deleteExistingFeatures = FALSE)
          if (length(mn)) srs <- import(srs, spectra = mn, ms_column_name = "group_idx", deleteExistingFeatures = FALSE)
          pv <- length(miss); miss <- setdiff(exp_ids, as.character(featuresId(srs, type = "xcms")))
          if (length(miss) >= pv) break
        }
        run(srs, formulaIdParams = fid, predictParams = predictParam(),
            structureDbSearchParams = structureDbSearchParam(structureSearchDbs = c("BIO")),
            recompute = TRUE, wait = TRUE)
        b_sum <- summary(srs, result.type = "structure") |>
          mutate(across(any_of(c("confidenceExactMatch", "confidenceApproxMatch")), as.numeric),
                 across(any_of(sirius_id_cols), as.character),
                 polarity = pol_lookup[as.character(externalFeatureId)])
        b_full <- results(srs, result.type = "structureDb", topFormula = 8, topStructure = 20,
                          return.type = "data.frame") |>
          mutate(across(any_of(sirius_id_cols), as.character),
                 polarity = pol_lookup[as.character(xcms_fts)])
        write.csv(b_sum,  sum_cf(i),  row.names = FALSE)
        write.csv(b_full, full_cf(i), row.names = FALSE)
        message("    batch ", i, "/", length(batches), " done (", np, "+", nn, " feats, try ", att, ")")
        TRUE
      }, error = function(e) {
        message("    batch ", i, " try ", att, " FAILED: ", substr(conditionMessage(e), 1, 90))
        try(shutdown(srs), silent = TRUE); srs <<- NULL; FALSE
      })
      if (isTRUE(ok)) { done <- TRUE; break }
      Sys.sleep(15)  # let the killed daemon release before a fresh one starts
    }
    if (!done) warning("    batch ", i, "/", length(batches), " GAVE UP after ", max_try,
                       " tries - left un-cached; rerun to resume")
  }
  try(shutdown(srs), silent = TRUE)

  # Aggregate whatever batch caches exist. `complete` = every batch is cached, so
  # the caller only writes the (lab, mixture) "done" marker on a full result.
  complete   <- all(vapply(seq_along(batches),
                           function(i) file.exists(sum_cf(i)) && file.exists(full_cf(i)), logical(1)))
  sum_files  <- sort(list.files(cache_dir, sprintf("^%s_bsum_[0-9]+\\.csv$",  projectId), full.names = TRUE))
  full_files <- sort(list.files(cache_dir, sprintf("^%s_bfull_[0-9]+\\.csv$", projectId), full.names = TRUE))
  read_bind  <- function(fs) if (!length(fs)) NULL else bind_rows(lapply(fs, function(f)
    read.csv(f, stringsAsFactors = FALSE) |>
      mutate(across(any_of(c("confidenceExactMatch", "confidenceApproxMatch")), as.numeric),
             across(any_of(c(sirius_id_cols, "polarity", "adduct")), as.character))))
  list(summary = read_bind(sum_files), full = read_bind(full_files), complete = complete)
}

summary_results    <- setNames(vector("list", length(mixtures)), mixtures)
full_results       <- setNames(vector("list", length(mixtures)), mixtures)
extraction_summary <- list()

for (m in mixtures) {
  summary_results[[m]] <- list()
  full_results[[m]]    <- list()
  for (lab in labs) {
    sum_path  <- sirius_sum_path(lab, m)
    full_path <- sirius_full_path(lab, m)

    # provenance_only: skip anything whose provenance already exists -> the pass is resumable
    if (provenance_only && file.exists(provenance_path(lab, m))) {
      message("Provenance cached: ", lab, " / mix ", m); next
    }
    # provenance_only must otherwise re-extract even when the annotation is cached (that's the point)
    if (!provenance_only && file.exists(sum_path) && file.exists(full_path)) {
      message("Cached: ", lab, " / mix ", m)
      summary_results[[m]][[lab]] <- read.csv(sum_path, stringsAsFactors = FALSE) |>
        mutate(across(any_of(c("confidenceExactMatch", "confidenceApproxMatch")), as.numeric),
               across(any_of(c(sirius_id_cols, "lab", "polarity", "adduct")), as.character))
      full_results[[m]][[lab]] <- read.csv(full_path, stringsAsFactors = FALSE) |>
        mutate(across(any_of(c(sirius_id_cols, "lab", "polarity", "adduct")), as.character))
      if (file.exists(extraction_count_path(lab, m)))
        extraction_summary[[paste(m, lab)]] <-
          read.csv(extraction_count_path(lab, m), stringsAsFactors = FALSE,
                   colClasses = c(mixture = "character", lab = "character")) |> as_tibble()
      next
    }

    message("Extracting: ", lab, " / mix ", m)
    sp <- extract_all_ms2(lab, study_group, m)
    extraction_summary[[paste(m, lab)]] <- tibble(
      mixture = m, lab = lab,
      n_features_pos = sp$n_features_pos, n_features_neg = sp$n_features_neg,
      n_ms2_pos = sp$n_ms2_pos, n_ms2_neg = sp$n_ms2_neg)
    write.csv(extraction_summary[[paste(m, lab)]],
              extraction_count_path(lab, m), row.names = FALSE)

    # Which scans (file + scan no + RT + precursor + CE) back each feature -> spectral-library provenance.
    prov <- rbind(spectra_prov(sp$spectra_pos), spectra_prov(sp$spectra_neg))
    if (!is.null(prov)) write.csv(prov, provenance_path(lab, m), row.names = FALSE)
    if (isTRUE(provenance_only)) {
      message("  provenance only: ", NROW(prov), " scans / ",
              length(unique(prov$feature_id)), " features -> ", basename(provenance_path(lab, m)))
      rm(sp); gc(); next
    }

    if (length(sp$spectra_pos) == 0L && length(sp$spectra_neg) == 0L) {
      message("  no MS2 - writing empty cache")
      write.csv(tibble(lab = character(0)), sum_path,  row.names = FALSE)
      write.csv(tibble(lab = character(0)), full_path, row.names = FALSE)
      rm(sp); gc(); next
    }

    res <- run_sirius_for_lab(
      sp_pos       = sp$spectra_pos,
      sp_neg       = sp$spectra_neg,
      username     = sirius_user, password = sirius_pass,
      projectId    = paste0("ring_trial_allms2_", lab, "_mix_", m),
      path         = output_dir,
      base_port    = sirius_port,
      # cembio + icl all-scans import many spectra/feature -> smaller batches so a batch doesn't
      # choke SIRIUS (the stall that killed the earlier run). The richest pass is 1 spectrum/feature,
      # so it keeps the full batch size and runs fast.
      batch_size   = if (representative == "primary" && lab %in% c("cembio", "icl")) 30L else sirius_batch_size)

    summary_results[[m]][[lab]] <- if (is.null(res$summary)) tibble(lab = character(0)) else res$summary |> mutate(lab = lab)
    full_results[[m]][[lab]]    <- if (is.null(res$full))    tibble(lab = character(0)) else res$full    |> mutate(lab = lab)

    # Only write the (lab, mixture) "done" cache when EVERY batch completed. A
    # partial run (some batches gave up on a blip) is left un-cached so the next
    # run resumes the missing batches rather than trusting an incomplete result.
    if (isTRUE(res$complete)) {
      write.csv(summary_results[[m]][[lab]], sum_path,  row.names = FALSE)
      write.csv(full_results[[m]][[lab]],    full_path, row.names = FALSE)
      message("Done: ", lab, " / mix ", m)
    } else {
      message("INCOMPLETE: ", lab, " / mix ", m,
              " - some batches un-cached; rerun to resume")
    }

    rm(sp, res); gc()
  }
}

print(bind_rows(extraction_summary))
message("ALL DONE: mixtures ", paste(mixtures, collapse = ", "))
