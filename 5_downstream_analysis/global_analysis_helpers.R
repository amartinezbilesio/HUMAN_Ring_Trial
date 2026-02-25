#' Helper functions for global analysis
#' 
#' This script contains reusable functions for the inter-laboratory
#' LC-MS reproducibility analysis.

# ==== Polarity Iteration Helpers ====

#' Run function for both polarities and combine results
run_for_polarities <- function(fn, ...) {
  bind_rows(
    fn(1L, "Positive", ...),
    fn(0L, "Negative", ...)
  )
}

#' Run function for both polarities (single argument version)
run_for_polarities_simple <- function(fn) {
  bind_rows(fn(1L), fn(0L))
}

# ==== Lab Pairs Definition ====
LAB_PAIRS <- list(
  c("afekta", "cembio"),
  c("afekta", "hmgu"),
  c("afekta", "icl"),
  c("cembio", "hmgu"),
  c("cembio", "icl"),
  c("hmgu", "icl")
)

# ==== Variance Decomposition ====

#' Compute variance decomposition from linear model
#' @param data Data frame with lab, mixture columns
#' @param metric_name Column name of metric to analyze
#' @return List with model, r_squared, variance_table, coefficients
check_metric_quality <- function(data, metric_name) {
  if (metric_name %in% c("tic", "bpc")) {
    model <- lm(log2(get(metric_name)) ~ lab + mixture, data = data)
  } else {
    model <- lm(get(metric_name) ~ lab + mixture, data = data)
  }
  
  r_squared <- summary(model)$r.squared
  an <- anova(model)
  ss_total <- sum(an[, "Sum Sq"])
  
  variance_data <- data.frame(
    Source = c("Lab Variance", "Mixture Variance", "Residual Variance"),
    Percentage = c(
      an["lab", "Sum Sq"] / ss_total * 100,
      an["mixture", "Sum Sq"] / ss_total * 100,
      an["Residuals", "Sum Sq"] / ss_total * 100
    )
  )
  
  coefs <- coef(model)
  coef_df <- data.frame(
    coefficient = names(coefs),
    estimate = coefs,
    type = case_when(
      grepl("^lab", names(coefs)) ~ "Lab Effect",
      grepl("^mixture", names(coefs)) ~ "Mixture Effect",
      TRUE ~ "Intercept"
    )
  ) |> filter(type != "Intercept")
  
  list(
    model = model,
    r_squared = r_squared,
    variance_table = variance_data,
    coefficients = coef_df
  )
}

#' Get variance percentages from linear model
get_variance <- function(data, value_col) {
  f <- as.formula(paste0("log2(", value_col, ") ~ lab + mixture"))
  an <- anova(lm(f, data = data))
  ss_total <- sum(an[, "Sum Sq"])
  data.frame(
    Mixture_Var = an["mixture", "Sum Sq"] / ss_total * 100,
    Lab_Var = an["lab", "Sum Sq"] / ss_total * 100
  )
}

# ==== ICC Computation ====

#' Compute intraclass correlation coefficient
compute_icc <- function(data, metric_name) {
  if (grepl("ratio", metric_name, ignore.case = TRUE)) {
    data$value <- data[[metric_name]]
  } else {
    data$value <- log2(data[[metric_name]])
  }
  
  between_var <- data |>
    group_by(mixture) |>
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") |>
    pull(mean_val) |>
    var(na.rm = TRUE)
  
  within_var <- data |>
    group_by(mixture) |>
    summarise(var_val = var(value, na.rm = TRUE), .groups = "drop") |>
    pull(var_val) |>
    mean(na.rm = TRUE)
  
  data.frame(
    Metric = metric_name,
    ICC = between_var / (between_var + within_var),
    Between_mixture_var = between_var,
    Within_mixture_var = within_var
  )
}

# ==== RT-Binned Correlation Analysis ====

#' Compute RT-binned correlations for TIC data across labs
#' @param tic_data Chromatograms object
#' @param bin_size RT bin size in seconds (default 5)
#' @return List with correlation matrices per polarity-mixture
compute_rt_binned_correlations <- function(tic_data, bin_size = 5) {
  rt_range <- range(rtime(tic_data[1L])[[1]])
  breaks <- seq(rt_range[1], rt_range[2], by = bin_size)
  
  polarities <- unique(tic_data$polarity)
  mixtures <- unique(tic_data$mixture)
  
  results <- list()
  
  for (pol in polarities) {
    for (mix in mixtures) {
      subset_tic <- tic_data[which(tic_data$polarity == pol & tic_data$mixture == mix)]
      if (length(subset_tic) == 0) next
      
      int_matrix <- data.frame(
        afekta = numeric(length(breaks) - 1),
        icl = numeric(length(breaks) - 1),
        hmgu = numeric(length(breaks) - 1),
        cembio = numeric(length(breaks) - 1)
      )
      
      for (i in seq_along(subset_tic)) {
        lab_name <- subset_tic$lab[i]
        peaks <- peaksData(subset_tic[i])[[1]]
        int_matrix[[lab_name]] <- tapply(
          peaks$intensity, 
          cut(peaks$rtime, breaks = breaks), 
          sum, na.rm = TRUE
        )
      }
      
      cor_mat <- cor(int_matrix, method = "pearson", use = "pairwise.complete.obs")
      polarity_label <- ifelse(pol == 1, "pos", "neg")
      results[[paste0(polarity_label, "_", mix)]] <- cor_mat
    }
  }
  results
}

#' Extract pairwise correlations from correlation results
extract_pairwise_cors <- function(correlation_results) {
  pairwise_cors <- data.frame()
  
  for (result_name in names(correlation_results)) {
    cor_mat <- correlation_results[[result_name]]
    for (pair in LAB_PAIRS) {
      pairwise_cors <- rbind(pairwise_cors, data.frame(
        mixture = result_name,
        pair = paste(pair[1], pair[2], sep = ":"),
        correlation = cor_mat[pair[1], pair[2]]
      ))
    }
  }
  
  pairwise_cors |>
    separate(mixture, into = c("polarity", "mix"), sep = "_", extra = "merge") |>
    mutate(polarity = factor(polarity, levels = c("neg", "pos")))
}

# ==== Plotting Helpers ====

#' Create residual plot from model with lab and mixture styling
#' @param model Linear model object
#' @param title Plot title
#' @param data Original data used for the model (optional, for lab/mixture styling)
plot_residuals <- function(model, title = "Residual Plot", data = NULL) {
  plot_df <- data.frame(Fitted = fitted(model), Residuals = residuals(model))
  
  if (!is.null(data) && "lab" %in% names(data) && "mixture" %in% names(data)) {
    plot_df$lab <- data$lab
    plot_df$mixture <- as.factor(data$mixture)
    p <- ggplot(plot_df, aes(x = Fitted, y = Residuals, color = mixture, shape = lab)) +
      geom_point(alpha = 0.7, size = 2.5) +
      scale_color_viridis_d(option = "turbo") +
      scale_shape_manual(values = c(16, 17, 15, 18)) +
      labs(color = "Mixture", shape = "Laboratory")
  } else {
    p <- ggplot(plot_df, aes(x = Fitted, y = Residuals)) +
      geom_point(alpha = 0.6, size = 2)
  }
  
  p + geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = title, x = "Fitted Values", y = "Residuals") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

#' Plot violin with jitter for intensity metrics
plot_violin <- function(data, metric) {
  if (metric %in% c("bpc", "tic")) {
    data <- data |> mutate(across(all_of(metric), log2))
    ylab <- paste(metric, "(log2)")
    title <- paste("Violin plot of", metric, "(log2 scale)")
  } else {
    ylab <- metric
    title <- paste("Violin plot of", metric)
  }
  
  ggplot(data, aes(x = lab, y = .data[[metric]], fill = polarity)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_jitter(aes(color = as.factor(mixture)), width = 0.2, size = 1.5, alpha = 0.7) +
    scale_color_viridis_d(option = "turbo") +
    facet_wrap(~ data_level) +
    theme_minimal() +
    labs(title = title, x = "Laboratory", y = ylab, color = "Mixture")
}

#' Plot correlation heatmap
plot_correlation_heatmap <- function(cor_data, title) {
  ggplot(cor_data, aes(x = pair, y = mix, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = round(correlation, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0.5, limits = c(0, 1)) +
    facet_wrap(~polarity, ncol = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = "Lab Pair", y = "Mixture", fill = "Pearson r")
}

#' Plot correlation distribution (beeswarm)
plot_correlation_distribution <- function(cor_data, title) {
  ggplot(cor_data, aes(x = pair, y = correlation, color = pair)) +
    geom_violin(alpha = 0.3, show.legend = FALSE) +
    ggbeeswarm::geom_beeswarm(size = 2, alpha = 0.7, show.legend = FALSE) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, 
                 color = "black", show.legend = FALSE) +
    facet_wrap(~polarity, ncol = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, subtitle = "Each point represents one mixture",
         x = "Lab Pair", y = "Pearson Correlation")
}

#' Plot quartile variance distribution
plot_quartile_variance <- function(data, metric_type = "TIC") {
  quartile_map <- c(
    "Q1" = "Q1 (Polar)", "Q2" = "Q2 (Semi-polar)",
    "Q3" = "Q3 (Lipids)", "Q4" = "Q4 (Wash)"
  )
  
  plot_data <- data |>
    mutate(Residual_Var = 100 - (Mixture_Var + Lab_Var)) |>
    filter(grepl(metric_type, Metric)) |>
    mutate(Quartile = case_when(
      grepl("Q1", Metric) ~ quartile_map["Q1"],
      grepl("Q2", Metric) ~ quartile_map["Q2"],
      grepl("Q3", Metric) ~ quartile_map["Q3"],
      grepl("Q4", Metric) ~ quartile_map["Q4"]
    )) |>
    pivot_longer(cols = c(Mixture_Var, Lab_Var, Residual_Var),
                 names_to = "Source", values_to = "Variance") |>
    mutate(Source = factor(Source,
                           levels = c("Residual_Var", "Lab_Var", "Mixture_Var"),
                           labels = c("Unexplained (Noise)", "Laboratory Bias", "Mixture (Biology)")))
  
  ggplot(plot_data, aes(x = Quartile, y = Variance, fill = Source)) +
    geom_col(width = 0.7, color = "white") +
    facet_wrap(~Polarity) +
    scale_fill_manual(values = c("gray85", "#E7B800", "#00AFBB")) +
    labs(title = paste(metric_type, "Variance Distribution Across RT Quartiles"),
         y = "Variance Contribution (%)", x = NULL) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top", panel.grid.major.x = element_blank())
}

# ==== Data Loading Helpers ====

#' Load peak tables for all labs
load_peak_tables <- function(base_path) {
  list(
    afekta = read.csv(file.path(base_path, "object/detected_peaks_afekta_HE.csv")),
    hmgu = read.csv(file.path(base_path, "object/detected_peaks_hmgu_HE.csv")),
    icl = read.csv(file.path(base_path, "object/detected_peaks_icl_HE.csv")),
    cembio = read.csv(file.path(base_path, "object/detected_peaks_cembio_HE.csv"))
  )
}

#' Join intensity data from all labs to consensus table
join_intensity_data <- function(ct_filtered, peak_tables) {
  ct_filtered |>
    select(consensus_id, sample, idx_afekta, 
           chrom_peak_id_afekta, chrom_peak_id_icl, 
           chrom_peak_id_hmgu, chrom_peak_id_cembio) |>
    left_join(peak_tables$afekta |> select(chrom_peak_id, into, rt),
              by = c("chrom_peak_id_afekta" = "chrom_peak_id")) |>
    rename(intensity_afekta = into, rt_afekta = rt) |>
    left_join(peak_tables$icl |> select(chrom_peak_id, into, rt),
              by = c("chrom_peak_id_icl" = "chrom_peak_id")) |>
    rename(intensity_icl = into, rt_icl = rt) |>
    left_join(peak_tables$hmgu |> select(chrom_peak_id, into, rt),
              by = c("chrom_peak_id_hmgu" = "chrom_peak_id")) |>
    rename(intensity_hmgu = into, rt_hmgu = rt) |>
    left_join(peak_tables$cembio |> select(chrom_peak_id, into, rt),
              by = c("chrom_peak_id_cembio" = "chrom_peak_id")) |>
    rename(intensity_cembio = into, rt_cembio = rt)
}

#' Assign RT quartiles to peaks
assign_rt_quartiles <- function(intensity_table, global_breaks) {
  intensity_table |>
    mutate(
      mean_rt = rowMeans(select(pick(everything()), starts_with("rt_")), na.rm = TRUE),
      quartile = case_when(
        mean_rt >= global_breaks[1] & mean_rt < global_breaks[2] ~ "Q1",
        mean_rt >= global_breaks[2] & mean_rt < global_breaks[3] ~ "Q2",
        mean_rt >= global_breaks[3] & mean_rt < global_breaks[4] ~ "Q3",
        mean_rt >= global_breaks[4] & mean_rt <= global_breaks[5] ~ "Q4",
        TRUE ~ NA_character_
      )
    )
}
