## Extract tables and printed cat() output from every rendered HTML so
## we can see what landed without manually opening each file.
files <- c(
  prod  = "3_annotation_manual/HE/3_Sirius_curation/sirius_annotation_all_ms2.html",
  naps  = "5_downstream_analysis/3_naps_extraction.html",
  rtaln = "5_downstream_analysis/4_rt_alignment.html",
  norm  = "5_downstream_analysis/5_normalization.html",
  tier  = "5_downstream_analysis/6_confidence_tier.html",
  icc   = "5_downstream_analysis/7_icc_reproducibility.html",
  level = "5_downstream_analysis/8_data_level_progression.html",
  pcdig = "5_downstream_analysis/9_per_compound_diagnostics.html"
)

extract_one <- function(nm, path) {
  if (!file.exists(path)) { cat("MISSING:", path, "\n"); return(invisible()) }
  cat("\n\n############# ", nm, " :: ", basename(path), " #############\n")
  txt <- paste(readLines(path, warn = FALSE), collapse = " ")

  # Headers
  hdrs <- regmatches(txt, gregexpr("<h[1-3][^>]*>[^<]+</h[1-3]>", txt))[[1]]
  hdrs <- gsub("<[^>]+>", "", hdrs)
  cat("\n=== Section headers ===\n")
  for (h in hdrs) cat("  -", h, "\n")

  # Tables
  tabs <- regmatches(txt, gregexpr("<table[\\s\\S]*?</table>", txt, perl = TRUE))[[1]]
  cat("\n=== Tables (", length(tabs), ") ===\n", sep = "")
  shown <- 0L
  for (i in seq_along(tabs)) {
    t <- tabs[i]
    cap <- regmatches(t, regexpr("<caption[^>]*>[^<]+</caption>", t))
    cap_txt <- if (length(cap)) gsub("<[^>]+>", "", cap) else "(no caption)"
    rows <- regmatches(t, gregexpr("<tr[^>]*>[\\s\\S]*?</tr>", t, perl = TRUE))[[1]]
    if (length(rows) > 25L || length(rows) < 1L) {
      cat("    [", length(rows), " rows, skipped] ", cap_txt, "\n", sep = "")
      next
    }
    shown <- shown + 1L
    cat("\n  TABLE ", i, ": ", cap_txt, "\n", sep = "")
    for (r in rows) {
      cells <- regmatches(r, gregexpr("<t[hd][^>]*>[\\s\\S]*?</t[hd]>", r, perl = TRUE))[[1]]
      cells <- gsub("<[^>]+>", "", cells)
      cells <- gsub("[\r\n]+", " ", cells)
      cells <- gsub("&nbsp;", " ", cells)
      cells <- gsub(" +", " ", cells)
      cells <- gsub("&amp;", "&", cells)
      cells <- gsub("&lt;", "<", cells)
      cells <- gsub("&gt;", ">", cells)
      cells <- trimws(cells)
      if (length(cells)) cat("    ", paste(cells, collapse = " | "), "\n")
    }
  }

  # Printed outputs (cat / print blocks in <pre>)
  pres <- regmatches(txt, gregexpr("<pre[^>]*>[\\s\\S]*?</pre>", txt, perl = TRUE))[[1]]
  cat("\n=== Printed outputs (", length(pres), " <pre> blocks) ===\n", sep = "")
  pshown <- 0L
  for (p in pres) {
    p <- gsub("<[^>]+>", "", p)
    p <- gsub("&gt;", ">", p); p <- gsub("&lt;", "<", p)
    p <- gsub("&amp;", "&", p); p <- gsub("&quot;", '"', p)
    p <- gsub("&#39;", "'", p)
    p <- trimws(p)
    if (!nchar(p)) next
    first <- strsplit(p, "\n")[[1]][1]
    # Skip code echoes
    if (grepl("^(library\\(|setwd|source|theme_|geom_|aes\\(|scale_|labs\\(|ggplot|knitr::)",
              first)) next
    if (grepl("^(#\\||<-|=|%>%|\\|>)\\s", first)) next
    if (!grepl("[0-9]", p)) next
    pshown <- pshown + 1L
    if (pshown > 20L) { cat("    [", length(pres) - 20, " more outputs skipped]\n"); break }
    if (nchar(p) > 1500L) p <- paste0(substr(p, 1, 1500), " [...]")
    cat("\n  OUTPUT ", pshown, ":\n", p, "\n", sep = "")
  }
}

for (nm in names(files)) extract_one(nm, files[[nm]])
