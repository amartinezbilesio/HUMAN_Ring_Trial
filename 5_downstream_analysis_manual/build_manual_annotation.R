## Build feat_rts.csv + annotation_identity.csv (canonical downstream schema)
## from the MANUAL CURATION xlsx (one per lab). These replace the SIRIUS
## all-MS2 producer's outputs so the downstream pipeline runs on the
## manually-curated identifications instead.
##
## Column mapping (curated -> canonical):
##   target_ChEBI.name        -> compound
##   target_InChIKey[1:14]    -> inchikey_2d
##   mixture (underscore)     -> mixture
##   adduct (pipe-sep, unnest)-> adduct ; polarity from adduct sign (+/-)
##   mz_<adduct>              -> ionMass (per-adduct m/z column)
##   rtmed                    -> rt_sec   (chromatographic apex, not MS2-trigger)
##   best_ConfidenceScoreExact-> cosmic / confidenceExactMatch
##   best_structurePerIdRank  -> best_structure_rank
##   target_formula           -> molecularFormula / expected_formula
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr); library(here)
})
labs    <- c("afekta", "cembio", "hmgu", "icl")
base    <- here("3_annotation_manual", "HE", "3_Sirius_curation")
out_dir <- here("5_downstream_analysis_manual", "object")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

two    <- function(x) toupper(substr(as.character(x), 1, 14))
pol_of <- function(add) ifelse(grepl("\\+$", add), "pos", "neg")

cur <- bind_rows(lapply(labs, function(l) {
  d <- as.data.frame(read_xlsx(file.path(base, paste0(l, "_lab_annotations_sirius_curated.xlsx"))),
                     check.names = FALSE)
  d$lab <- l; d
}))
cur$mixture <- gsub("\\.", "_", as.character(cur$mixture))

## one row per (curated compound x adduct)
long <- cur |>
  mutate(.adducts = str_split(adduct, fixed("|"))) |>
  unnest(.adducts) |>
  mutate(adduct = str_trim(.adducts)) |>
  filter(nzchar(adduct))

## ionMass = value of the matching mz_<adduct> column for that row
long$ionMass <- vapply(seq_len(nrow(long)), function(i) {
  col <- paste0("mz_", long$adduct[i])
  if (col %in% names(long)) suppressWarnings(as.numeric(long[[col]][i])) else NA_real_
}, numeric(1))

unmatched <- setdiff(unique(long$adduct), sub("^mz_", "", grep("^mz_", names(long), value = TRUE)))
if (length(unmatched)) message("  !! adducts with no mz_ column: ", paste(unmatched, collapse = ", "))

long <- long |>
  mutate(polarity    = pol_of(adduct),
         compound    = `target_ChEBI.name`,
         inchikey_2d = two(`target_InChIKey`),
         rt_sec      = suppressWarnings(as.numeric(rtmed)),
         cosmic      = suppressWarnings(as.numeric(`best_ConfidenceScoreExact`)),
         strank      = suppressWarnings(as.numeric(`best_structurePerIdRank`)))

## ---- feat_rts.csv (per feature; pipeline aggregates per adduct) ----
feat_rts <- long |>
  filter(!is.na(ionMass), !is.na(rt_sec), !is.na(compound)) |>
  transmute(mixture, lab, polarity,
            xcms_fts = paste(lab, mixture, polarity, compound, adduct, sep = "::"),
            compound, inchikey_2d, rt_sec, ionMass, adduct, cosmic)
write.csv(feat_rts, file.path(out_dir, "feat_rts.csv"), row.names = FALSE)

## ---- annotation_identity.csv (per mixture,lab,compound,polarity) ----
fin_min <- function(x) { x <- x[is.finite(x)]; if (length(x)) min(x) else NA_real_ }
fin_max <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x) else NA_real_ }
annotation_identity <- long |>
  filter(!is.na(compound)) |>
  group_by(mixture, lab, compound, inchikey_2d, polarity) |>
  summarise(best_formula_rank   = NA_real_,
            best_structure_rank = fin_min(strank),
            adducts             = paste(unique(adduct), collapse = ";"),
            molecularFormula    = dplyr::first(`target_formula`),
            expected_formula    = dplyr::first(`target_formula`),
            expected_inchikey   = dplyr::first(`target_InChIKey`),
            mz_observed         = suppressWarnings(median(ionMass, na.rm = TRUE)),
            confidenceExactMatch  = fin_max(cosmic),
            confidenceApproxMatch = NA_real_,
            confidence_tier       = "unfiltered",
            .groups = "drop")
write.csv(annotation_identity, file.path(out_dir, "annotation_identity.csv"), row.names = FALSE)

cat(sprintf("feat_rts.csv            : %d rows (%d labs, %d mixtures)\n",
            nrow(feat_rts), dplyr::n_distinct(feat_rts$lab), dplyr::n_distinct(feat_rts$mixture)))
cat(sprintf("annotation_identity.csv : %d rows (%d distinct compounds)\n",
            nrow(annotation_identity), dplyr::n_distinct(annotation_identity$compound)))
cat("\nfeat_rts per lab:\n"); print(table(feat_rts$lab, feat_rts$polarity))
cat("\nafekta mix 1_5 feat_rts (subset-test target):\n")
print(feat_rts |> filter(lab == "afekta", mixture == "1_5") |>
        select(compound, polarity, adduct, ionMass, rt_sec, cosmic) |> head(12))
