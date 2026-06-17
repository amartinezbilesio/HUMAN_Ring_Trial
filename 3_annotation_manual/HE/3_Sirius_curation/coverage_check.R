suppressPackageStartupMessages({ library(readxl); library(here) })
a <- commandArgs(trailingOnly = TRUE); lab <- a[1]; mix <- a[2]
two <- function(x) toupper(substr(as.character(x), 1, 14))
cat(sprintf("\n############ %s mix %s ############\n", lab, mix))

## EXPECTED panel: HE, this mixture
raw <- as.data.frame(read_xlsx(here("1_preprocessing","standards.xlsx"), col_names = FALSE))
d <- raw[-1, ]
sel <- !is.na(d[[11]]) & d[[11]] == "Human Endosome" & gsub("\\.","_", as.character(d[[12]])) == mix & !is.na(d[[8]])
ep <- data.frame(k2 = two(d[[8]][sel]), nm = d[[1]][sel], stringsAsFactors = FALSE)
ep <- ep[!duplicated(ep$k2) & nzchar(ep$k2), ]
cat("EXPECTED panel (HE mix", mix, "):", nrow(ep), "compounds\n")

## NEW all-MS2
D  <- here("3_annotation_manual","HE","3_Sirius_curation","results_all_ms2_pg")
sm <- read.csv(file.path(D, sprintf("sirius_summary_%s_mix_%s.csv", lab, mix)), stringsAsFactors = FALSE)
sm <- sm[!is.na(sm$inchiKey) & nzchar(sm$inchiKey), ]; sm$k2 <- two(sm$inchiKey)
conf <- tapply(sm$confidenceApproxMatch, sm$k2, function(v) suppressWarnings(max(v, na.rm = TRUE)))
newK <- unique(sm$k2)
hiK  <- function(th) names(conf)[is.finite(conf) & conf >= th]
cat(sprintf("NEW all-MS2: %d features, %d distinct compounds (any) | %d @>=0.5 | %d @>=0.7\n",
            length(unique(sm$externalFeatureId)), length(newK), length(hiK(.5)), length(hiK(.7))))

## MANUAL CURATION (ground truth, + MS2 status)
cf <- here("3_annotation_manual","HE","3_Sirius_curation", paste0(lab, "_lab_annotations_sirius_curated.xlsx"))
cx <- as.data.frame(read_xlsx(cf)); cx$mix <- gsub("\\.","_", as.character(cx$mixture)); cx$k2 <- two(cx$target_InChIKey)
m <- cx[cx$mix == mix & !is.na(cx$k2) & nzchar(cx$k2), ]
ms2cnt <- tapply(suppressWarnings(as.numeric(m$ms2_count)), m$k2, function(v) max(v, na.rm = TRUE))
manK <- names(ms2cnt); man_ms2 <- names(ms2cnt)[ms2cnt > 0]; man_ms1only <- names(ms2cnt)[!(ms2cnt > 0)]
cat(sprintf("MANUAL curation: %d compounds (%d with MS2, %d MS1-only)\n", length(manK), length(man_ms2), length(man_ms1only)))

## COVERAGE OF EXPECTED
ek <- ep$k2
cov_any <- intersect(ek, newK); cov_man <- intersect(ek, manK)
cat("\n--- coverage of EXPECTED panel ---\n")
cat(sprintf("  all-MS2 : %d / %d (any) | %d @>=0.5 | %d @>=0.7\n", length(cov_any), nrow(ep), length(intersect(ek, hiK(.5))), length(intersect(ek, hiK(.7)))))
cat(sprintf("  manual  : %d / %d\n", length(cov_man), nrow(ep)))
cat(sprintf("  expected found by all-MS2 that manual MISSED: %d\n", length(setdiff(cov_any, cov_man))))
cat(sprintf("  expected found by manual that all-MS2 MISSED: %d\n", length(setdiff(cov_man, cov_any))))

## the misses + MS2 status
miss <- setdiff(ek, newK)
cat("\n--- expected compounds MISSED by all-MS2 (and whether they have MS2) ---\n")
for (k in miss) {
  st <- if (k %in% man_ms2) "HAS MS2 (real gap)" else if (k %in% man_ms1only) "MS1-only (unreachable)" else "not in manual"
  cat(sprintf("  %-15s %-32s %s\n", k, ep$nm[match(k, ep$k2)], st))
}
