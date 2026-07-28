# Build the spectral-library sheet: which PANEL COMPOUND was annotated, by which METHOD,
# from which SCAN of which FILE — one row per (compound x contributing scan).
#
# Joins, per lab x mixture x pass:
#   results_all_ms2/<pass>/sirius_full_<lab>_mix_<m>.csv   feature (xcms_fts) -> structure
#   results_all_ms2/<pass>/provenance_<lab>_mix_<m>.csv    feature -> (file, scan, RT, precursor, CE)
#
# provenance_*.csv is written by run_sirius_cluster.R. If it's missing for a result, regenerate it
# WITHOUT re-annotating: set  provenance_only <- TRUE  in the runner and re-run those lab/mixtures
# (extraction is deterministic, so the feature ids match the existing sirius_full).
#
# NOTE hmgu group-2 is vial-mislabeled (file 2_N holds mixture 2_(N+1); file 2_10 holds 2_2), so panel
# matching uses the TRUE mixture — otherwise every hmgu group-2 annotation would be dropped.
suppressWarnings(suppressMessages({ source(here::here("1_preprocessing", "setup.R")); library(writexl) }))
meta$mix <- chartr(".", "_", meta$Mixture); two <- function(x) substr(x, 1, 14)
RES <- here::here("2_annotation", "sirius_all_ms2", "results")   # all-MS2 annotations + provenance
OUT <- here::here("3_library_generation", "spectral_library_provenance.xlsx")
MIX <- c(paste0("1_", 1:6), paste0("2_", 1:10), paste0("3_", 1:3))
LABS <- c("afekta", "cembio", "icl", "hmgu")

## the finalized representative actually used for a (pass, lab)
method_of <- function(pass, lab) {
  if (pass == "richest") return("richest-scan")
  switch(lab, afekta = "richest-per-CE", cembio = "all-scans", icl = "all-scans", hmgu = "richest-scan")
}
## `pass` column label: the second run is stored under results/richest but is reported as "secondary"
pass_label <- function(p) if (p == "richest") "secondary" else p
## hmgu group-2 vial relabel: file label -> true mixture (cyclic, nothing lost)
true_mix <- function(lab, mx) {
  if (lab != "hmgu") return(mx)
  if (mx == "2_10") return("2_2")
  p <- as.integer(strsplit(mx, "_")[[1]])
  if (p[1] == 2L && p[2] >= 2L && p[2] <= 9L) return(sprintf("2_%d", p[2] + 1L))
  mx
}
## panel of a mixture: 2D key -> (name, full InChIKey, formula)
panel_tbl <- function(mx) {
  mm <- meta[meta$mix == mx & !is.na(meta$InChIKey), ]
  if (!nrow(mm)) return(NULL)
  k <- two(mm$InChIKey); keep <- !duplicated(k)
  data.frame(k2 = k[keep], compound = mm[["ChEBI name"]][keep],
             target_InChIKey = mm$InChIKey[keep],
             target_formula = if ("formula" %in% names(mm)) mm$formula[keep] else NA_character_,
             stringsAsFactors = FALSE)
}

rows <- list(); gaps <- list()
for (pass in c("primary", "richest")) for (lab in LABS) {
  if (pass == "richest" && lab == "hmgu") next          # hmgu: primary == richest, no second pass
  for (mx in MIX) {
    ffull <- file.path(RES, pass, sprintf("sirius_full_%s_mix_%s.csv", lab, mx))
    fprov <- file.path(RES, pass, sprintf("provenance_%s_mix_%s.csv", lab, mx))
    if (!file.exists(ffull)) next                        # that pass simply wasn't run -> not a gap
    if (!file.exists(fprov)) {                           # annotated but provenance missing -> flag it
      gaps[[length(gaps) + 1L]] <- data.frame(pass = pass_label(pass), lab = lab, mixture = mx,
        issue = "provenance missing - rerun runner with provenance_only <- TRUE", stringsAsFactors = FALSE); next
    }
    d <- read.csv(ffull, stringsAsFactors = FALSE, check.names = FALSE)
    d <- d[!is.na(d$inchiKey) & nzchar(d$inchiKey), ]; if (!nrow(d)) next
    d$csi <- suppressWarnings(as.numeric(d$csiScore))
    o <- d[order(d$xcms_fts, -d$csi), ]; top <- o[!duplicated(o$xcms_fts), ]   # SIRIUS's #1 per feature
    tm <- true_mix(lab, mx); pt <- panel_tbl(tm); if (is.null(pt)) next
    top$k2 <- two(top$inchiKey)
    hit <- merge(top, pt, by = "k2")                     # keep only PANEL compounds
    if (!nrow(hit)) next
    prov <- read.csv(fprov, stringsAsFactors = FALSE, check.names = FALSE)
    # take only what we need from the SIRIUS side: `polarity` (and others) exist in BOTH tables and
    # would be silently renamed to .x/.y by merge(), leaving the columns below empty.
    keep <- c("xcms_fts", "k2", "compound", "target_InChIKey", "structureName", "molecularFormula", "adduct", "csi")
    hit2 <- hit[, intersect(keep, names(hit)), drop = FALSE]
    j <- merge(hit2, prov, by.x = "xcms_fts", by.y = "feature_id")
    if (!nrow(j)) next
    rows[[length(rows) + 1L]] <- data.frame(
      lab = lab, file_mixture = mx, true_mixture = tm, mislabeled = (tm != mx),
      compound = j$compound, target_InChIKey = j$target_InChIKey, inchikey_2d = j$k2,
      sirius_structure = j$structureName, molecular_formula = j$molecularFormula, adduct = j$adduct,
      method = method_of(pass, lab), pass = pass_label(pass),
      feature_id = j$xcms_fts, file = j$file, scan_index = j$scan_index,
      acquisition_num = j$acquisition_num, rt_sec = j$rt_sec,
      # precursor_mz_raw = the true experimental precursor of THIS scan (what a library needs);
      # precursor_mz_med = the feature median that was actually fed to SIRIUS. Older provenance
      # CSVs only have `precursor_mz` (the median); fall back to it so a mixed set still builds.
      precursor_mz_raw = if ("precursor_mz_raw" %in% names(j)) j$precursor_mz_raw else j$precursor_mz,
      precursor_mz_med = if ("precursor_mz_med" %in% names(j)) j$precursor_mz_med else j$precursor_mz,
      collision_energy = j$collision_energy, n_peaks = j$n_peaks, polarity = j$polarity,
      csi_score = j$csi, stringsAsFactors = FALSE)
  }
}

ann <- if (length(rows)) do.call(rbind, rows) else
  data.frame(note = "no provenance yet - run the runner with provenance_only <- TRUE")
if (length(rows)) ann <- ann[order(ann$lab, ann$true_mixture, ann$compound, ann$method), ]

## one row per compound x method (scans collapsed) - the quick index
summ <- if (length(rows)) {
  s <- aggregate(cbind(n_scans = seq_len(nrow(ann))) ~ lab + true_mixture + compound + method, data = ann, FUN = length)
  s[order(s$lab, s$true_mixture, s$compound), ]
} else data.frame(note = "none")

sheets <- list(annotations = ann, by_compound_method = summ)
if (length(gaps)) sheets$provenance_gaps <- do.call(rbind, gaps)
write_xlsx(sheets, OUT)
cat("rows:", if (length(rows)) nrow(ann) else 0,
    "| compounds:", if (length(rows)) length(unique(paste(ann$lab, ann$true_mixture, ann$compound))) else 0,
    "| provenance gaps:", length(gaps), "\n")
cat("written:", OUT, "\n")
