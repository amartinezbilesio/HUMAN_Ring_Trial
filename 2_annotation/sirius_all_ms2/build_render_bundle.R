# Reproducibly build recovery_bundle.csv, the committed CSV that lets results_summary.qmd render its
# recovery tables on a fresh clone WITHOUT the (gitignored) results/ CSVs. Written next to the qmd.
#   - recovery_bundle.csv : per (pass, lab, mix) and panel 2D-InChIKey, the best SIRIUS rank it
#                           reached (rank==1 is top-1; any presence is any-candidate).
# The hmgu vial-mislabel evidence is NOT precomputed here; results_summary.qmd recomputes it directly
# from the raw hmgu mzML so it stays reproducible. Run this whenever the annotation results change.
suppressPackageStartupMessages({library(here); library(readxl)})
RES <- here("2_annotation", "sirius_all_ms2", "results")
OUT <- here("2_annotation", "sirius_all_ms2")   # the CSV sits directly next to the qmd
two <- function(x) substr(x, 1, 14)

## the qmd only ever intersects annotations with panel compounds, so keep only panel-relevant
## 2D-keys (union of every mixture's panel). Lossless for recovery; shrinks the bundle ~250x.
meta <- as.data.frame(read_excel(here("1_preprocessing", "standards.xlsx"),
                                 skip = 1, .name_repair = "minimal"), check.names = FALSE)
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
