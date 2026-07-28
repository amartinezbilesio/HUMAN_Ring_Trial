# Reproducibly build the three committed CSVs that let results_summary.qmd render on a fresh clone
# WITHOUT the (gitignored) results/ CSVs or the raw mzML. All three are written next to the qmd.
#   - recovery_bundle.csv : per (pass, lab, mix) and panel 2D-InChIKey, the best SIRIUS rank it
#                           reached (rank==1 is top-1; any presence is any-candidate).
#   - hmgu_grid.csv       : file x panel raw-MS2 mass-trigger counts for the hmgu vial-mislabel
#                           check (Evidence 2), which otherwise needs the hmgu mzML.
#   - hmgu_injorder.csv   : hmgu acquisition timestamps + the true mixture each file carries.
# Run this whenever the results change; commit the CSVs so the qmd is self-contained.
suppressPackageStartupMessages({library(here); library(readxl)})
RES <- here("2_annotation", "sirius_all_ms2", "results")
OUT <- here("2_annotation", "sirius_all_ms2")   # the CSVs sit directly next to the qmd
two <- function(x) substr(x, 1, 14)

## the qmd only ever intersects annotations with panel compounds, so keep only panel-relevant
## 2D-keys (union of every mixture's panel). Lossless for recovery; shrinks the bundle ~250x.
meta <- as.data.frame(read_excel(here("1_preprocessing", "standards.xlsx"),
                                 skip = 1, .name_repair = "minimal"), check.names = FALSE)
meta$mix <- chartr(".", "_", meta$Mixture)          # same key the qmd builds (panel_masses needs it)
allpanel <- unique(two(meta$InChIKey[!is.na(meta$InChIKey)]))

## ---- recovery bundle (slim; replaces read_pass / ranks / pre) ----
## Per (pass, lab, mix) and panel-relevant 2D-key, the BEST rank that key reached (position among a
## feature's structure candidates ordered by score, min across the features it appears in). rank==1
## is SIRIUS's #1 structure (top-1); any presence is any-candidate; the value drives the rank section.
rows <- list()
for (pass in c("primary", "richest", "primary_noms1")) {
  d <- file.path(RES, pass)
  if (!dir.exists(d)) next
  for (f in list.files(d, "^sirius_full_.*\\.csv$", full.names = TRUE)) {
    lab <- sub("^sirius_full_([a-z]+)_mix_.*", "\\1", basename(f))
    mx  <- gsub("^sirius_full_[a-z]+_mix_|\\.csv$", "", basename(f))
    x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    x <- x[!is.na(x$inchiKey) & nzchar(x$inchiKey), ]
    if (!nrow(x)) next
    o <- x[order(x$xcms_fts, -suppressWarnings(as.numeric(x$csiScore))), ]
    o$rk <- ave(seq_len(nrow(o)), o$xcms_fts, FUN = seq_along)   # rank within each feature's candidates
    o$k2 <- two(o$inchiKey)
    o <- o[o$k2 %in% allpanel, ]                                 # panel-relevant keys only
    if (!nrow(o)) next
    a <- aggregate(rk ~ k2, o, min)                             # best rank each panel key reached
    rows[[length(rows) + 1L]] <- data.frame(pass = pass, lab = lab, mix = mx,
      inchikey_2d = a$k2, rank = a$rk, stringsAsFactors = FALSE)
  }
}
rb <- do.call(rbind, rows)
write.csv(rb, file.path(OUT, "recovery_bundle.csv"), row.names = FALSE)
cat(sprintf("recovery_bundle.csv: %d rows, %d (pass,lab,mix) groups\n",
            nrow(rb), nrow(unique(rb[c("pass","lab","mix")]))))

## ---- hmgu vial-mislabel raw-MS2 evidence (replaces the two mzML-reading chunks) ----
## Slow (reads all hmgu mzML). Skipped when a good hmgu_grid.csv already exists, so reruns that only
## refresh the recovery bundle stay fast; delete hmgu_grid.csv to force a recompute.
.gp <- file.path(OUT, "hmgu_grid.csv")
.have_grid <- file.exists(.gp) && any(read.csv(.gp)$n > 0)
suppressPackageStartupMessages(library(Spectra))
sp_path <- here("1_preprocessing", "hmgu", "HE"); mzml_dir <- file.path(sp_path, "mzml")
if (.have_grid) {
  cat("hmgu_grid.csv already present and non-empty; skipped raw-evidence recompute\n")
} else if (dir.exists(mzml_dir) && length(list.files(mzml_dir, "mzML$"))) {
  seq <- rbind(as.data.frame(read_excel(file.path(sp_path, "seq_pos.xlsx"))),
               as.data.frame(read_excel(file.path(sp_path, "seq_neg.xlsx"))))
  seq$filename <- paste0(seq$`Data File`, ".mzML")
  seq <- seq[seq$filename %in% list.files(mzml_dir), ]
  seq$label <- chartr(".", "_", sub("target$", "", sub(".*_", "", sub("_MS2$", "", seq$`Data File`))))
  mixes <- sprintf("2_%d", 1:10); pr <- 1.007276
  adducts <- c(pr, 22.989218, 18.033823, 38.963158, -pr, 34.969402, 44.998201)
  panel_masses <- function(mx) { mm <- meta[meta$mix == mx & !is.na(meta$InChIKey), ]
    e2 <- substr(mm$InChIKey, 1, 14); dup <- e2[duplicated(e2) | duplicated(e2, fromLast = TRUE)]
    keep <- !(e2 %in% dup); as.numeric(mm$M[keep][!duplicated(e2[keep])]) }
  load_ms2 <- function(label) { fs <- seq$filename[seq$label == label]
    s <- setBackend(Spectra(file.path(mzml_dir, fs), backend = MsBackendMzR()), MsBackendMemory())
    m <- s[s$msLevel == 2L]; m <- m[!is.na(m$precursorMz)]
    list(prec = m$precursorMz, bpk = vapply(peaksData(m), function(pp) if (nrow(pp)) max(pp[, "intensity"]) else 0, numeric(1))) }
  strong <- function(d, M) any(vapply(M + adducts, function(a) {
    hit <- which(abs(d$prec - a) <= a * 15 / 1e6); length(hit) > 0 && max(d$bpk[hit]) > 1e4 }, logical(1)))
  panels <- setNames(lapply(mixes, panel_masses), mixes)
  grid <- expand.grid(file = mixes, panel = mixes, stringsAsFactors = FALSE); grid$n <- NA_integer_
  acqv <- setNames(rep(NA_character_, length(mixes)), mixes)
  for (fl in mixes) {
    d <- load_ms2(fl)
    for (pn in mixes) grid$n[grid$file == fl & grid$panel == pn] <- sum(vapply(panels[[pn]], function(M) strong(d, M), logical(1)))
    f1 <- file.path(mzml_dir, seq$filename[seq$label == fl & grepl("_pos_", seq$filename)][1])
    con <- file(f1, "r"); txt <- readLines(con, n = 500, warn = FALSE); close(con)
    ts <- regmatches(txt, regexpr('startTimeStamp="[^"]+"', txt)); ts <- ts[nzchar(ts)]
    acqv[fl] <- if (length(ts)) sub('startTimeStamp="([^"]+)"', "\\1", ts[1]) else NA_character_
  }
  write.csv(grid, file.path(OUT, "hmgu_grid.csv"), row.names = FALSE)
  best <- do.call(rbind, lapply(split(grid, grid$file), function(g) g[which.max(g$n), ]))
  write.csv(data.frame(file = mixes, acquired = acqv[mixes], true_mixture = best$panel[match(mixes, best$file)]),
            file.path(OUT, "hmgu_injorder.csv"), row.names = FALSE)
  cat("hmgu_grid.csv + hmgu_injorder.csv written\n")
} else cat("hmgu mzml not available; skipped raw-evidence precompute\n")
