#!/usr/bin/env Rscript
## build_consensus_no_bin.R
##
## Compute consensus peaks across labs WITHOUT m/z binning.
## Uses pairwise mclosest mappings (all combinations) and merges via graph
## connected components. Designed to run on HPC nodes; supports parallel chunking.
##
## Usage:
## Rscript build_consensus_no_bin.R --input_dir ../object --out_prefix ../object/consensus_all --cores 8 --chunk_size 5000
##
## Notes:
## - The script expects CSV files named: detected_peaks_<lab>_HE.csv for labs: afekta, hmgu, icl, cembio
## - Output: <out_prefix>.(rds|csv) with consensus table and an edges file
## - Recommended: run per-pair on separate jobs or run once on a multicore node.
##   For typical data sizes here, start with 8-16 cores and 4-8 GB per core.

suppressPackageStartupMessages({
  library(data.table)
  library(MetaboCoreUtils)
  library(BiocParallel)
  library(igraph)
})

args <- commandArgs(trailingOnly = TRUE)
message(
  "Command-line args:",
  if (length(args)) paste(args, collapse = " ") else " <none>"
)

parse_arg <- function(name, default = NULL) {
  # Accept both "--name=value" and "--name value" styles
  eq_pat <- paste0("--", name, "=")
  space_pat <- paste0("--", name)

  # --name=value
  m <- args[startsWith(args, eq_pat)]
  if (length(m) > 0) {
    return(sub(eq_pat, "", m[1]))
  }

  # --name value
  i <- which(args == space_pat)
  if (length(i) > 0 && length(args) >= i[1] + 1) {
    return(args[i[1] + 1])
  }

  default
}

# Defaults assume script is run from repository root
input_dir <- parse_arg(
  "input_dir",
  file.path("5_downstream_analysis", "object")
)
out_prefix <- parse_arg(
  "out_prefix",
  file.path("5_downstream_analysis", "object", "consensus_all")
)
cores <- as.integer(parse_arg(
  "cores",
  Sys.getenv("SLURM_CPUS_ON_NODE", unset = 8)
))
chunk_size <- as.integer(parse_arg("chunk_size", 5000))
ppm_tol <- as.numeric(parse_arg("ppm_tol", 10))
rt_tol <- as.numeric(parse_arg("rt_tol", 10))

input_dir <- normalizePath(input_dir, mustWork = FALSE)
out_prefix <- normalizePath(out_prefix, mustWork = FALSE)
message("Resolved input dir: ", input_dir)
message("Out prefix: ", out_prefix)

# Verify expected input files exist; if not, try to discover them recursively
labs <- c("afekta", "hmgu", "icl", "cembio")
expected_files <- paste0("detected_peaks_", labs, "_HE.csv")
missing <- expected_files[!file.exists(file.path(input_dir, expected_files))]
if (length(missing) > 0) {
  message(
    "Some expected files are missing in the provided input_dir: ",
    input_dir
  )
  message("Missing: ", paste(missing, collapse = ", "))

  # Search recursively under current working directory for candidate directories
  candidates <- list.files(
    path = ".",
    pattern = "detected_peaks_.*_HE\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(candidates) > 0) {
    cand_dirs <- unique(dirname(candidates))
    found <- NULL
    for (d in cand_dirs) {
      if (all(file.exists(file.path(d, expected_files)))) {
        found <- d
        break
      }
    }
    if (!is.null(found)) {
      input_dir <- normalizePath(found)
      message("Auto-detected input_dir: ", input_dir)
      missing <- expected_files[
        !file.exists(file.path(input_dir, expected_files))
      ]
    } else {
      message(
        "Did not find a single directory containing all expected files. Found candidates in: ",
        paste(cand_dirs, collapse = ", ")
      )
    }
  } else {
    message(
      "No candidate detected_peaks files found under current directory (recursive search)."
    )
  }

  if (length(missing) > 0) {
    stop(
      "Missing required detected peaks CSV files. Either place them in '",
      input_dir,
      "' or provide --input_dir with the correct path. Missing: ",
      paste(missing, collapse = ", ")
    )
  }
}
message("Cores: ", cores, "  Chunk size: ", chunk_size)

labs <- c("afekta", "hmgu", "icl", "cembio")
files <- file.path(input_dir, paste0("detected_peaks_", labs, "_HE.csv"))
names(files) <- labs

for (f in files) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}


# --- NEW: Per-sample consensus matching ---
# Assumes all input tables have a 'sample' column (numeric or character)
all_samples <- Reduce(union, lapply(dt_list, function(dt) unique(dt$sample)))
all_samples <- sort(all_samples)
message("Found ", length(all_samples), " unique samples.")

# BiocParallel parameter
BPPARAM <- MulticoreParam(workers = cores)

for (sample_val in all_samples) {
  message("Processing sample: ", sample_val)
  # Subset each lab's table to this sample
  dt_sample <- lapply(dt_list, function(dt) dt[dt$sample == sample_val, ])
  names(dt_sample) <- labs

  # If all are empty, skip
  if (all(sapply(dt_sample, nrow) == 0)) {
    next
  }

  # helper: midpoints for mz and rt
  get_mid <- function(dt) {
    mz <- (dt$mzmin + dt$mzmax) / 2
    rt <- (dt$rtmin + dt$rtmax) / 2
    list(mz = mz, rt = rt)
  }
  peaks_mid <- lapply(dt_sample, get_mid)

  compute_mclosest_chunks <- function(
    from_lab,
    to_lab,
    chunk_size,
    peak_cols = 1:4
  ) {
    df_from <- dt_sample[[from_lab]]
    df_to <- dt_sample[[to_lab]]
    n_from <- nrow(df_from)
    if (n_from == 0 || nrow(df_to) == 0) {
      return(data.table(from = integer(), to = integer()))
    }

    chunks <- split(seq_len(n_from), ceiling(seq_len(n_from) / chunk_size))

    res_chunks <- bplapply(
      chunks,
      function(idx_chunk) {
        x <- as.matrix(df_from[idx_chunk, peak_cols, with = FALSE])
        y <- as.matrix(df_to[, peak_cols, with = FALSE])
        m <- tryCatch(
          mclosest(x, y, tolerance = c(10, 10, 0, 0), ppm = c(0, 0, 10, 10)),
          error = function(e) rep(NA_integer_, nrow(x))
        )
        data.table(from = idx_chunk, to = ifelse(is.na(m), NA_integer_, m))
      },
      BPPARAM = BPPARAM
    )

    data.table::rbindlist(res_chunks)
  }

  # Build edge list by computing all pairwise mclosest (both directions)
  pairs <- t(combn(labs, 2))
  edges <- data.table(from = character(), to = character())

  for (k in seq_len(nrow(pairs))) {
    a <- pairs[k, 1]
    b <- pairs[k, 2]
    message("Computing mclosest: ", a, " -> ", b)
    m_ab <- compute_mclosest_chunks(a, b, chunk_size)
    if (nrow(m_ab) > 0) {
      m_ab <- m_ab[!is.na(to)]
      if (nrow(m_ab) > 0) {
        edges <- rbind(
          edges,
          data.table(
            from = paste0(a, ":", m_ab$from),
            to = paste0(b, ":", m_ab$to)
          )
        )
      }
    }

    message("Computing mclosest: ", b, " -> ", a)
    m_ba <- compute_mclosest_chunks(b, a, chunk_size)
    if (nrow(m_ba) > 0) {
      m_ba <- m_ba[!is.na(to)]
      if (nrow(m_ba) > 0) {
        edges <- rbind(
          edges,
          data.table(
            from = paste0(b, ":", m_ba$from),
            to = paste0(a, ":", m_ba$to)
          )
        )
      }
    }
  }

  message("Total raw edges: ", nrow(edges))

  # Add isolated nodes (peaks never matched) to vertex set
  all_nodes <- unlist(lapply(labs, function(lb) {
    paste0(lb, ":", seq_len(nrow(dt_sample[[lb]])))
  }))
  vertices <- data.table(name = unique(c(all_nodes, edges$from, edges$to)))

  # Build graph and extract components
  g <- graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)
  comp <- components(g)
  membership <- comp$membership

  node_df <- data.table(name = names(membership), comp = as.integer(membership))
  node_df[, c("lab", "idx") := tstrsplit(name, ":", fixed = TRUE)]
  node_df[, idx := as.integer(idx)]

  # Aggregate per component
  consensus_list <- node_df[,
    .(
      labs_present = list(unique(lab)),
      idxs = list(.SD[, .(lab, idx)]),
      n_labs = uniqueN(lab)
    ),
    by = comp
  ]

  # For each component create a single-row mapping
  consensus_rows <- lapply(seq_len(nrow(consensus_list)), function(i) {
    row <- consensus_list[i]
    members <- row$idxs[[1]]
    idx_afekta <- if ("afekta" %in% members$lab) {
      members[idx == min(idx) & lab == "afekta"]$idx[1]
    } else {
      NA_integer_
    }
    idx_hmgu <- if ("hmgu" %in% members$lab) {
      members[idx == min(idx) & lab == "hmgu"]$idx[1]
    } else {
      NA_integer_
    }
    idx_icl <- if ("icl" %in% members$lab) {
      members[idx == min(idx) & lab == "icl"]$idx[1]
    } else {
      NA_integer_
    }
    idx_cembio <- if ("cembio" %in% members$lab) {
      members[idx == min(idx) & lab == "cembio"]$idx[1]
    } else {
      NA_integer_
    }

    mz_vals <- c()
    rt_vals <- c()
    if (!is.na(idx_afekta)) {
      mz_vals <- c(mz_vals, peaks_mid$afekta$mz[idx_afekta])
      rt_vals <- c(rt_vals, peaks_mid$afekta$rt[idx_afekta])
    }
    if (!is.na(idx_hmgu)) {
      mz_vals <- c(mz_vals, peaks_mid$hmgu$mz[idx_hmgu])
      rt_vals <- c(rt_vals, peaks_mid$hmgu$rt[idx_hmgu])
    }
    if (!is.na(idx_icl)) {
      mz_vals <- c(mz_vals, peaks_mid$icl$mz[idx_icl])
      rt_vals <- c(rt_vals, peaks_mid$icl$rt[idx_icl])
    }
    if (!is.na(idx_cembio)) {
      mz_vals <- c(mz_vals, peaks_mid$cembio$mz[idx_cembio])
      rt_vals <- c(rt_vals, peaks_mid$cembio$rt[idx_cembio])
    }

    consensus_mz <- if (length(mz_vals) > 0) {
      mean(mz_vals, na.rm = TRUE)
    } else {
      NA_real_
    }
    consensus_rt <- if (length(rt_vals) > 0) {
      mean(rt_vals, na.rm = TRUE)
    } else {
      NA_real_
    }

    data.table(
      comp = row$comp,
      idx_afekta = idx_afekta,
      idx_hmgu = idx_hmgu,
      idx_icl = idx_icl,
      idx_cembio = idx_cembio,
      n_labs_matched = row$n_labs,
      consensus_mz = consensus_mz,
      consensus_rt = consensus_rt
    )
  })

  consensus_dt <- data.table::rbindlist(consensus_rows)

  # Map chrom_peak_id values back
  consensus_dt[,
    chrom_peak_id_afekta := ifelse(
      !is.na(idx_afekta),
      dt_sample$afekta$chrom_peak_id[idx_afekta],
      NA_character_
    )
  ]
  consensus_dt[,
    chrom_peak_id_hmgu := ifelse(
      !is.na(idx_hmgu),
      dt_sample$hmgu$chrom_peak_id[idx_hmgu],
      NA_character_
    )
  ]
  consensus_dt[,
    chrom_peak_id_icl := ifelse(
      !is.na(idx_icl),
      dt_sample$icl$chrom_peak_id[idx_icl],
      NA_character_
    )
  ]
  consensus_dt[,
    chrom_peak_id_cembio := ifelse(
      !is.na(idx_cembio),
      dt_sample$cembio$chrom_peak_id[idx_cembio],
      NA_character_
    )
  ]

  # Add consensus id
  consensus_dt[, consensus_id := .I]
  consensus_dt[, sample := sample_val]

  # Reorder columns
  setcolorder(
    consensus_dt,
    c(
      "sample",
      "consensus_id",
      "consensus_mz",
      "consensus_rt",
      "n_labs_matched",
      "idx_afekta",
      "chrom_peak_id_afekta",
      "idx_hmgu",
      "chrom_peak_id_hmgu",
      "idx_icl",
      "chrom_peak_id_icl",
      "idx_cembio",
      "chrom_peak_id_cembio",
      "comp"
    )
  )

  # Save per-sample results immediately
  out_rds <- paste0(out_prefix, "_sample", sample_val, ".rds")
  out_csv <- paste0(out_prefix, "_sample", sample_val, ".csv")
  fwrite(consensus_dt, out_csv)
  fwrite(edges, paste0(out_prefix, "_sample", sample_val, "_edges.csv"))
  saveRDS(list(consensus = consensus_dt, edges = edges), file = out_rds)

  message("Done. Sample ", sample_val, " consensus rows: ", nrow(consensus_dt))
  message("Saved: ", out_rds, " and ", out_csv)
}
