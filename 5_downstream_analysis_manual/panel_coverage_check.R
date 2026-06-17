suppressPackageStartupMessages({ library(readxl); library(here) })
two <- function(x) toupper(substr(as.character(x), 1, 14))

## FULL HE (Human Endosome) panel from standards.xlsx
raw <- as.data.frame(read_xlsx(here("1_preprocessing", "standards.xlsx"), col_names = FALSE))
d   <- raw[-1, ]
he  <- d[!is.na(d[[11]]) & d[[11]] == "Human Endosome", ]
he_k   <- two(he[[8]]); he_mix <- gsub("\\.", "_", as.character(he[[12]]))
ok     <- nzchar(he_k) & !is.na(he_k)
cat("HE panel — unique compounds (2D InChIKey):", length(unique(he_k[ok])), "\n")
cat("HE panel — unique (mixture, compound) spike-ins:",
    nrow(unique(data.frame(m = he_mix[ok], k = he_k[ok]))), "\n")
cat("HE panel — rows in standards.xlsx:", nrow(he), " | mixtures:", length(unique(he_mix[ok])), "\n")

## MANUAL curation coverage
ai <- read.csv(here("5_downstream_analysis_manual", "object", "annotation_identity.csv"),
               stringsAsFactors = FALSE, colClasses = c(mixture = "character"))
cat("\nmanual curation — unique compounds:", length(unique(ai$compound)),
    "| unique inchikey_2d:", length(unique(ai$inchikey_2d)), "\n")
cat("manual curation — unique (mixture, compound):", nrow(unique(ai[, c("mixture", "compound")])), "\n")

## panel coverage (mixture-aware)
he_mk <- unique(paste(he_mix[ok], he_k[ok]))
ai_mk <- unique(paste(ai$mixture, toupper(ai$inchikey_2d)))
cat("\n(mixture,compound) panel spike-ins curated by >=1 lab:", sum(he_mk %in% ai_mk), "of", length(he_mk), "\n")
cat("panel spike-ins NOT in manual curation:", sum(!(he_mk %in% ai_mk)), "\n")
cat("manual-curation (mixture,compound) NOT in panel (off-panel):", sum(!(ai_mk %in% he_mk)), "\n")
