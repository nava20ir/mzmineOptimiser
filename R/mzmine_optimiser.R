
#' Read MZmine optimiser TOML configuration
#'
#' @param config_file Path to a TOML config file.
#' @return A nested list from `RcppTOML::parseTOML()` with normalised empty values.
#' @export
read_mzmine_config <- function(config_file = "config.toml") {
  cfg <- RcppTOML::parseTOML(config_file)

  required_sections <- c("paths", "qc", "optimization", "parameters")
  missing_sections <- setdiff(required_sections, names(cfg))
  if (length(missing_sections) > 0) {
    stop("Missing TOML section(s): ", paste(missing_sections, collapse = ", "), call. = FALSE)
  }

  required_paths <- c("base_batch", "mzmine_exe", "optim_root", "final_optimised_batch")
  missing_paths <- required_paths[!nzchar(vapply(required_paths, function(x) cfg$paths[[x]] %||% "", character(1)))]
  if (length(missing_paths) > 0) {
    stop("Missing paths value(s): ", paste(missing_paths, collapse = ", "), call. = FALSE)
  }

  cfg$paths$mzuser <- cfg$paths$mzuser %||% ""
  if (identical(cfg$paths$mzuser, "")) cfg$paths$mzuser <- NA_character_

  cfg$qc$qc_all <- unlist(cfg$qc$qc_all, use.names = FALSE)
  cfg$qc$final_input_files <- unlist(cfg$qc$final_input_files %||% character(0), use.names = FALSE)
  if (length(cfg$qc$final_input_files) == 0) cfg$qc$final_input_files <- NULL

  cfg$optimization$n_phase1_runs <- as.integer(cfg$optimization$n_phase1_runs)
  cfg$optimization$top_n_phase2 <- as.integer(cfg$optimization$top_n_phase2)
  cfg$optimization$feature_csv_name <- cfg$optimization$feature_csv_name %||% "features.csv"
  cfg$optimization$random_seed <- as.integer(cfg$optimization$random_seed %||% 123)

  cfg
}

#' Write a template config.toml file
#'
#' @param output_file Path where the template should be written.
#' @param overwrite Overwrite existing file?
#' @return Invisibly returns `output_file`.
#' @export
write_template_config <- function(output_file = "config.toml", overwrite = FALSE) {
  if (file.exists(output_file) && !overwrite) {
    stop("File already exists: ", output_file, ". Use overwrite = TRUE.", call. = FALSE)
  }
  template <- system.file("extdata", "config.toml", package = "mzmineOptimiser")
  if (!nzchar(template)) stop("Template config not found inside package.", call. = FALSE)
  file.copy(template, output_file, overwrite = overwrite)
  invisible(output_file)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

read_batch <- function(batch_file) {
  xml2::read_xml(batch_file)
}

write_batch <- function(doc, output_file) {
  xml2::write_xml(doc, output_file, options = "format")
  output_file
}

get_batchsteps <- function(doc) {
  xml2::xml_find_all(doc, ".//batchstep")
}

get_batchstep_by_module <- function(doc, module_name) {
  steps <- get_batchsteps(doc)
  steps[xml2::xml_attr(steps, "module_name") == module_name]
}

get_param_direct <- function(node, param_name) {
  xml2::xml_find_first(node, sprintf("./parameter[@name=\"%s\"]", param_name))
}

get_param_anywhere <- function(node, param_name) {
  xml2::xml_find_first(node, sprintf(".//parameter[@name=\"%s\"]", param_name))
}

set_param_text <- function(node, param_name, value, anywhere = FALSE) {
  p <- if (anywhere) get_param_anywhere(node, param_name) else get_param_direct(node, param_name)
  if (inherits(p, "xml_missing")) stop(glue::glue("Parameter not found: {param_name}"), call. = FALSE)
  xml2::xml_text(p) <- as.character(value)
  invisible(node)
}

set_range_param <- function(node, param_name, min_value = NULL, max_value = NULL) {
  p <- get_param_anywhere(node, param_name)
  if (inherits(p, "xml_missing")) stop(glue::glue("Range parameter not found: {param_name}"), call. = FALSE)
  if (!is.null(min_value)) {
    min_node <- xml2::xml_find_first(p, "./min")
    if (inherits(min_node, "xml_missing")) stop(glue::glue("No <min> found in parameter: {param_name}"), call. = FALSE)
    xml2::xml_text(min_node) <- as.character(min_value)
  }
  if (!is.null(max_value)) {
    max_node <- xml2::xml_find_first(p, "./max")
    if (inherits(max_node, "xml_missing")) stop(glue::glue("No <max> found in parameter: {param_name}"), call. = FALSE)
    xml2::xml_text(max_node) <- as.character(max_value)
  }
  invisible(node)
}

set_current_file <- function(param_node, new_file) {
  cf <- xml2::xml_find_first(param_node, "./current_file")
  if (inherits(cf, "xml_missing")) {
    xml2::xml_add_child(param_node, "current_file", new_file)
  } else {
    xml2::xml_text(cf) <- new_file
  }
  invisible(param_node)
}

replace_import_files <- function(doc, files) {
  import_step <- get_batchstep_by_module(doc, "Import MS data")
  if (length(import_step) != 1) stop("Expected exactly one 'Import MS data' step.", call. = FALSE)
  file_param <- get_param_direct(import_step[[1]], "File names")
  if (inherits(file_param, "xml_missing")) stop("Could not find parameter 'File names'.", call. = FALSE)
  old_file_nodes <- xml2::xml_find_all(file_param, "./file")
  xml2::xml_remove(old_file_nodes)
  for (f in files) xml2::xml_add_child(file_param, "file", f)
  invisible(doc)
}

get_mass_detection_step <- function(doc, target = c("MS1", "MSn")) {
  target <- match.arg(target)
  steps <- get_batchstep_by_module(doc, "Mass detection")
  if (length(steps) == 0) stop("No Mass detection steps found.", call. = FALSE)
  for (st in steps) {
    ms_filter <- xml2::xml_find_first(st, ".//parameter[@name='MS level filter']")
    if (!inherits(ms_filter, "xml_missing")) {
      selected <- xml2::xml_attr(ms_filter, "selected") %||% ""
      txt <- stringr::str_trim(xml2::xml_text(ms_filter))
      if (target == "MS1" && (stringr::str_detect(selected, "MS1") || txt == "1")) return(st)
      if (target == "MSn" && (stringr::str_detect(selected, "MSn") || txt %in% c("2", "3"))) return(st)
    }
  }
  stop(glue::glue("Could not identify {target} Mass detection step."), call. = FALSE)
}

set_ms1_noise <- function(doc, noise_level) {
  ms1_step <- get_mass_detection_step(doc, "MS1")
  auto_noise <- xml2::xml_find_first(ms1_step, ".//parameter[@name='Mass detector']/module[@name='Auto']/parameter[@name='Noise level']")
  if (inherits(auto_noise, "xml_missing")) stop("Could not find MS1 Auto Noise level.", call. = FALSE)
  xml2::xml_text(auto_noise) <- as.character(noise_level)
  invisible(doc)
}

set_chrom_builder_params <- function(doc, min_consecutive_scans, min_consecutive_intensity, min_absolute_height) {
  step <- get_batchstep_by_module(doc, "Chromatogram builder")
  if (length(step) != 1) stop("Expected exactly one Chromatogram builder step.", call. = FALSE)
  st <- step[[1]]
  set_param_text(st, "Minimum consecutive scans", min_consecutive_scans)
  set_param_text(st, "Minimum intensity for consecutive scans", min_consecutive_intensity)
  set_param_text(st, "Minimum absolute height", min_absolute_height)
  invisible(doc)
}

set_resolver_params <- function(doc, chromatographic_threshold, resolver_min_height, peak_top_edge_ratio, resolver_min_scans, peak_duration_min, peak_duration_max) {
  step <- get_batchstep_by_module(doc, "Local minimum feature resolver")
  if (length(step) != 1) stop("Expected exactly one Local minimum feature resolver step.", call. = FALSE)
  st <- step[[1]]
  set_param_text(st, "Chromatographic threshold", chromatographic_threshold)
  set_param_text(st, "Minimum absolute height", resolver_min_height)
  set_param_text(st, "Min ratio of peak top/edge", peak_top_edge_ratio)
  set_param_text(st, "Minimum scans (data points)", resolver_min_scans)
  set_range_param(st, "Peak duration range (min/mobility)", min_value = peak_duration_min, max_value = peak_duration_max)
  invisible(doc)
}

apply_parameter_set <- function(doc, pars) {
  set_ms1_noise(doc, noise_level = pars$ms1_noise)
  set_chrom_builder_params(doc,
    min_consecutive_scans = pars$min_consecutive_scans,
    min_consecutive_intensity = pars$min_consecutive_intensity,
    min_absolute_height = pars$chrom_min_absolute_height)
  set_resolver_params(doc,
    chromatographic_threshold = pars$chromatographic_threshold,
    resolver_min_height = pars$resolver_min_height,
    peak_top_edge_ratio = pars$peak_top_edge_ratio,
    resolver_min_scans = pars$resolver_min_scans,
    peak_duration_min = pars$peak_duration_min,
    peak_duration_max = pars$peak_duration_max)
  invisible(doc)
}

set_csv_export_path <- function(doc, csv_file) {
  steps <- get_batchstep_by_module(doc, "Export to CSV file (modular)")
  if (length(steps) < 1) stop("Could not find CSV export step.", call. = FALSE)
  filename_param <- get_param_direct(steps[[1]], "Filename")
  if (inherits(filename_param, "xml_missing")) stop("Could not find CSV export Filename parameter.", call. = FALSE)
  set_current_file(filename_param, csv_file)
  invisible(doc)
}

set_autosave_path <- function(doc, batch_save_file) {
  steps <- get_batchstep_by_module(doc, "Save current batch")
  if (length(steps) < 1) return(invisible(doc))
  p <- get_param_direct(steps[[1]], "Save batch to")
  if (!inherits(p, "xml_missing")) set_current_file(p, batch_save_file)
  invisible(doc)
}

make_phase1_detection_doc <- function(doc) {
  steps <- get_batchsteps(doc)
  keep_modules <- c("Import MS data", "Mass detection", "Chromatogram builder", "Smoothing", "Local minimum feature resolver", "Export to CSV file (modular)", "Save current batch")
  for (st in steps) if (!xml2::xml_attr(st, "module_name") %in% keep_modules) xml2::xml_remove(st)
  invisible(doc)
}

make_phase2_qc_validation_doc <- function(doc) {
  steps <- get_batchsteps(doc)
  keep_modules <- c("Import MS data", "Mass detection", "Chromatogram builder", "Smoothing", "Local minimum feature resolver", "13C isotope filter (formerly: isotope grouper)", "Isotopic peaks finder", "Join aligner", "Feature list rows filter", "Export to CSV file (modular)", "Save current batch")
  for (st in steps) if (!xml2::xml_attr(st, "module_name") %in% keep_modules) xml2::xml_remove(st)
  invisible(doc)
}

run_mzmine_batch <- function(batch_file, mzmine_exe, mzuser = NA_character_) {
  args <- c()
  if (!is.na(mzuser) && nzchar(mzuser)) args <- c(args, "-user", mzuser)
  args <- c(args, "-batch", batch_file)
  message("Running MZmine:")
  message(mzmine_exe, " ", paste(args, collapse = " "))
  system2(command = mzmine_exe, args = args, stdout = TRUE, stderr = TRUE)
}

pick_col <- function(df, patterns) {
  nms <- names(df); nms_low <- stringr::str_to_lower(nms)
  for (p in patterns) {
    hits <- which(stringr::str_detect(nms_low, stringr::str_to_lower(p)))
    if (length(hits) > 0) return(nms[hits[1]])
  }
  NA_character_
}

read_mzmine_feature_csv <- function(csv_file) {
  if (!fs::file_exists(csv_file)) {
    warning("CSV file does not exist: ", csv_file)
    return(NULL)
  }
  df <- tryCatch(readr::read_csv(csv_file, show_col_types = FALSE, guess_max = 10000), error = function(e) {
    warning("Could not read CSV: ", csv_file); warning(e$message); NULL
  })
  if (is.null(df) || nrow(df) == 0) return(NULL)
  mz_col <- pick_col(df, c("^row m/z$", "row mz", "m/z", "mz"))
  rt_col <- pick_col(df, c("^row retention time$", "retention time", "rt"))
  height_col <- pick_col(df, c("height", "peak height", "maximum height", "max height"))
  area_col <- pick_col(df, c("area", "peak area"))
  ms2_col <- pick_col(df, c("ms/ms", "ms2", "fragment", "spectra", "scan"))
  df2 <- df
  df2$mz_value <- if (!is.na(mz_col)) suppressWarnings(as.numeric(df2[[mz_col]])) else NA_real_
  df2$rt_value <- if (!is.na(rt_col)) suppressWarnings(as.numeric(df2[[rt_col]])) else NA_real_
  df2$height_value <- if (!is.na(height_col)) suppressWarnings(as.numeric(df2[[height_col]])) else NA_real_
  df2$area_value <- if (!is.na(area_col)) suppressWarnings(as.numeric(df2[[area_col]])) else NA_real_
  df2$has_ms2 <- if (!is.na(ms2_col)) !is.na(df2[[ms2_col]]) & as.character(df2[[ms2_col]]) != "" & as.character(df2[[ms2_col]]) != "NA" else FALSE
  df2
}

identify_abundance_cols <- function(df) {
  numeric_cols <- names(df)[purrr::map_lgl(df, is.numeric)]
  exclude_patterns <- c("row", "mz", "m/z", "retention", "rt", "id", "charge", "score", "similarity", "cos", "network", "annotation", "formula", "mass", "mobility", "ccs", "isotope", "height_value", "area_value", "mz_value", "rt_value")
  numeric_cols[!stringr::str_detect(stringr::str_to_lower(numeric_cols), paste(exclude_patterns, collapse = "|"))]
}

score_phase1 <- function(features) {
  if (is.null(features) || nrow(features) == 0) return(tibble::tibble(n_features = 0, n_ms2_features = 0, ms2_fraction = 0, median_height = NA_real_, phase1_score = -Inf))
  n_features <- nrow(features)
  n_ms2_features <- sum(features$has_ms2, na.rm = TRUE)
  ms2_fraction <- ifelse(n_features > 0, n_ms2_features / n_features, 0)
  median_height <- stats::median(features$height_value, na.rm = TRUE)
  phase1_score <- log10(n_features + 1) + 10 * ms2_fraction + ifelse(is.finite(median_height), 0.5 * log10(median_height + 1), 0) - 0.0005 * n_features
  tibble::tibble(n_features = n_features, n_ms2_features = n_ms2_features, ms2_fraction = ms2_fraction, median_height = median_height, phase1_score = phase1_score)
}

score_phase2 <- function(features) {
  if (is.null(features) || nrow(features) == 0) return(tibble::tibble(n_features = 0, n_ms2_features = 0, ms2_fraction = 0, median_height = NA_real_, median_detection_frequency = 0, median_qc_rsd = NA_real_, phase2_score = -Inf))
  n_features <- nrow(features)
  n_ms2_features <- sum(features$has_ms2, na.rm = TRUE)
  ms2_fraction <- ifelse(n_features > 0, n_ms2_features / n_features, 0)
  median_height <- stats::median(features$height_value, na.rm = TRUE)
  abundance_cols <- identify_abundance_cols(features)
  if (length(abundance_cols) == 0) {
    warning("No abundance columns detected. Phase 2 score will use fallback metrics.")
    phase2_score <- log10(n_features + 1) + 10 * ms2_fraction + ifelse(is.finite(median_height), 0.5 * log10(median_height + 1), 0) - 0.0005 * n_features
    return(tibble::tibble(n_features = n_features, n_ms2_features = n_ms2_features, ms2_fraction = ms2_fraction, median_height = median_height, median_detection_frequency = NA_real_, median_qc_rsd = NA_real_, phase2_score = phase2_score))
  }
  abundance_matrix <- features |>
    dplyr::select(dplyr::all_of(abundance_cols)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(as.numeric(.x), 0))) |>
    as.matrix()
  detection_frequency <- rowMeans(abundance_matrix > 0, na.rm = TRUE)
  row_mean <- rowMeans(abundance_matrix, na.rm = TRUE)
  row_sd <- apply(abundance_matrix, 1, stats::sd, na.rm = TRUE)
  row_rsd <- ifelse(row_mean > 0, row_sd / row_mean, NA_real_)
  median_detection_frequency <- stats::median(detection_frequency, na.rm = TRUE)
  median_qc_rsd <- stats::median(row_rsd, na.rm = TRUE)
  phase2_score <- 50 * median_detection_frequency + 20 * ms2_fraction - 10 * median_qc_rsd + ifelse(is.finite(median_height), 0.5 * log10(median_height + 1), 0) - 0.0005 * n_features
  tibble::tibble(n_features = n_features, n_ms2_features = n_ms2_features, ms2_fraction = ms2_fraction, median_height = median_height, median_detection_frequency = median_detection_frequency, median_qc_rsd = median_qc_rsd, phase2_score = phase2_score)
}

make_param_grid <- function(cfg) {
  p <- cfg$parameters
  tidyr::crossing(
    ms1_noise = unlist(p$ms1_noise),
    min_consecutive_scans = unlist(p$min_consecutive_scans),
    min_consecutive_intensity = unlist(p$min_consecutive_intensity),
    chrom_min_absolute_height = unlist(p$chrom_min_absolute_height),
    chromatographic_threshold = unlist(p$chromatographic_threshold),
    resolver_min_height = unlist(p$resolver_min_height),
    peak_top_edge_ratio = unlist(p$peak_top_edge_ratio),
    resolver_min_scans = unlist(p$resolver_min_scans),
    peak_duration_min = unlist(p$peak_duration_min),
    peak_duration_max = unlist(p$peak_duration_max)
  )
}

make_phase1_batch <- function(base_batch, output_batch, qc_file, output_csv, pars) {
  doc <- read_batch(base_batch)
  replace_import_files(doc, files = qc_file)
  make_phase1_detection_doc(doc)
  apply_parameter_set(doc, pars)
  set_csv_export_path(doc, output_csv)
  set_autosave_path(doc, stringr::str_replace(output_batch, "\\.mzbatch$", "_autosaved.mzbatch"))
  write_batch(doc, output_batch)
}

make_phase2_batch <- function(base_batch, output_batch, qc_files, output_csv, pars) {
  doc <- read_batch(base_batch)
  replace_import_files(doc, files = qc_files)
  make_phase2_qc_validation_doc(doc)
  apply_parameter_set(doc, pars)
  set_csv_export_path(doc, output_csv)
  set_autosave_path(doc, stringr::str_replace(output_batch, "\\.mzbatch$", "_autosaved.mzbatch"))
  write_batch(doc, output_batch)
}

make_final_optimised_full_batch <- function(base_batch, output_batch, best_row, final_input_files = NULL) {
  doc <- read_batch(base_batch)
  if (!is.null(final_input_files)) replace_import_files(doc, files = final_input_files)
  apply_parameter_set(doc, best_row)
  final_csv <- stringr::str_replace(output_batch, "\\.mzbatch$", "_feature_table.csv")
  set_csv_export_path(doc, final_csv)
  set_autosave_path(doc, stringr::str_replace(output_batch, "\\.mzbatch$", "_autosaved.mzbatch"))
  write_batch(doc, output_batch)
}

#' Run the complete MZmine optimisation workflow from config.toml
#'
#' @param config_file Path to TOML configuration file.
#' @return A list containing final best parameters and output file paths.
#' @export
run_mzmine_optimisation <- function(config_file = "config.toml") {
  cfg <- read_mzmine_config(config_file)
  base_batch <- cfg$paths$base_batch
  mzmine_exe <- cfg$paths$mzmine_exe
  mzuser <- cfg$paths$mzuser
  optim_root <- cfg$paths$optim_root
  final_optimised_batch <- cfg$paths$final_optimised_batch
  qc_single <- cfg$qc$qc_single
  qc_all <- cfg$qc$qc_all
  final_input_files <- cfg$qc$final_input_files
  n_phase1_runs <- cfg$optimization$n_phase1_runs
  top_n_phase2 <- cfg$optimization$top_n_phase2
  feature_csv_name <- cfg$optimization$feature_csv_name

  fs::dir_create(optim_root)

  param_grid_full <- make_param_grid(cfg)
  set.seed(cfg$optimization$random_seed)
  n_to_sample <- min(n_phase1_runs, nrow(param_grid_full))
  param_grid_phase1 <- dplyr::slice_sample(param_grid_full, n = n_to_sample)
  readr::write_csv(param_grid_phase1, fs::path(optim_root, "phase1_parameter_grid.csv"))
  message("Phase 1 runs: ", nrow(param_grid_phase1))

  phase1_results <- vector("list", nrow(param_grid_phase1))
  for (i in seq_len(nrow(param_grid_phase1))) {
    pars <- param_grid_phase1[i, ]
    run_id <- sprintf("phase1_run_%03d", i)
    run_dir <- fs::path(optim_root, "phase1", run_id)
    fs::dir_create(run_dir)
    batch_i <- fs::path(run_dir, paste0(run_id, ".mzbatch"))
    csv_i <- fs::path(run_dir, feature_csv_name)
    log_i <- fs::path(run_dir, "mzmine_log.txt")
    message("\n================================================")
    message("Phase 1: ", run_id, " of ", nrow(param_grid_phase1))
    message("================================================")
    make_phase1_batch(base_batch, batch_i, qc_single, csv_i, pars)
    mzmine_output <- tryCatch(run_mzmine_batch(batch_i, mzmine_exe, mzuser), error = function(e) { warning("MZmine failed for ", run_id, ": ", e$message); e$message })
    writeLines(as.character(mzmine_output), log_i)
    features <- read_mzmine_feature_csv(csv_i)
    score_tbl <- score_phase1(features)
    phase1_results[[i]] <- dplyr::bind_cols(tibble::tibble(run_id = run_id, batch_file = batch_i, csv_file = csv_i), pars, score_tbl)
    readr::write_csv(dplyr::bind_rows(phase1_results), fs::path(optim_root, "phase1_results_live.csv"))
    print(phase1_results[[i]])
  }

  phase1_tbl <- dplyr::bind_rows(phase1_results)
  readr::write_csv(phase1_tbl, fs::path(optim_root, "phase1_results_final.csv"))

  top_phase1 <- phase1_tbl |>
    dplyr::filter(is.finite(.data$phase1_score)) |>
    dplyr::arrange(dplyr::desc(.data$phase1_score)) |>
    dplyr::slice_head(n = top_n_phase2)
  if (nrow(top_phase1) == 0) stop("No successful Phase 1 run found.", call. = FALSE)
  readr::write_csv(top_phase1, fs::path(optim_root, "top_phase1_for_phase2.csv"))
  message("\nTop Phase 1 parameter sets selected for Phase 2:")
  print(top_phase1)

  phase2_results <- vector("list", nrow(top_phase1))
  for (i in seq_len(nrow(top_phase1))) {
    pars <- top_phase1[i, ]
    run_id <- sprintf("phase2_run_%03d", i)
    run_dir <- fs::path(optim_root, "phase2", run_id)
    fs::dir_create(run_dir)
    batch_i <- fs::path(run_dir, paste0(run_id, ".mzbatch"))
    csv_i <- fs::path(run_dir, feature_csv_name)
    log_i <- fs::path(run_dir, "mzmine_log.txt")
    message("\n================================================")
    message("Phase 2: ", run_id, " of ", nrow(top_phase1))
    message("================================================")
    make_phase2_batch(base_batch, batch_i, qc_all, csv_i, pars)
    mzmine_output <- tryCatch(run_mzmine_batch(batch_i, mzmine_exe, mzuser), error = function(e) { warning("MZmine failed for ", run_id, ": ", e$message); e$message })
    writeLines(as.character(mzmine_output), log_i)
    features <- read_mzmine_feature_csv(csv_i)
    score_tbl <- score_phase2(features)
    phase2_results[[i]] <- dplyr::bind_cols(
      tibble::tibble(run_id = run_id, batch_file = batch_i, csv_file = csv_i),
      pars |> dplyr::select(.data$ms1_noise, .data$min_consecutive_scans, .data$min_consecutive_intensity, .data$chrom_min_absolute_height, .data$chromatographic_threshold, .data$resolver_min_height, .data$peak_top_edge_ratio, .data$resolver_min_scans, .data$peak_duration_min, .data$peak_duration_max),
      score_tbl)
    readr::write_csv(dplyr::bind_rows(phase2_results), fs::path(optim_root, "phase2_results_live.csv"))
    print(phase2_results[[i]])
  }

  phase2_tbl <- dplyr::bind_rows(phase2_results)
  readr::write_csv(phase2_tbl, fs::path(optim_root, "phase2_results_final.csv"))

  best <- phase2_tbl |>
    dplyr::filter(is.finite(.data$phase2_score)) |>
    dplyr::arrange(dplyr::desc(.data$phase2_score)) |>
    dplyr::slice(1)
  if (nrow(best) == 0) stop("No successful Phase 2 run found.", call. = FALSE)
  readr::write_csv(best, fs::path(optim_root, "best_parameters_final.csv"))
  message("\nFinal best parameters:")
  print(best)

  make_final_optimised_full_batch(base_batch, final_optimised_batch, best, final_input_files)
  message("\nFinal optimised full batch written to:")
  message(final_optimised_batch)

  final_qc_only_batch <- stringr::str_replace(final_optimised_batch, "\\.mzbatch$", "_QC_all_only.mzbatch")
  make_final_optimised_full_batch(base_batch, final_qc_only_batch, best, qc_all)
  message("\nFinal optimised QC-validation full batch written to:")
  message(final_qc_only_batch)

  plot_dir <- fs::path(optim_root, "plots")
  fs::dir_create(plot_dir)

  p1 <- ggplot2::ggplot(phase1_tbl, ggplot2::aes(x = .data$n_features, y = .data$phase1_score)) + ggplot2::geom_point(size = 2) + ggplot2::theme_bw() + ggplot2::labs(title = "Phase 1: score vs feature count", x = "Number of features", y = "Phase 1 score")
  ggplot2::ggsave(fs::path(plot_dir, "phase1_score_vs_feature_count.png"), p1, width = 7, height = 5, dpi = 150)
  p2 <- ggplot2::ggplot(phase1_tbl, ggplot2::aes(x = factor(.data$ms1_noise), y = .data$phase1_score)) + ggplot2::geom_boxplot() + ggplot2::theme_bw() + ggplot2::labs(title = "Phase 1: effect of MS1 noise level", x = "MS1 noise level", y = "Phase 1 score")
  ggplot2::ggsave(fs::path(plot_dir, "phase1_score_by_ms1_noise.png"), p2, width = 7, height = 5, dpi = 150)
  p3 <- ggplot2::ggplot(phase2_tbl, ggplot2::aes(x = .data$median_detection_frequency, y = .data$phase2_score)) + ggplot2::geom_point(size = 3) + ggplot2::theme_bw() + ggplot2::labs(title = "Phase 2: score vs QC detection frequency", x = "Median QC detection frequency", y = "Phase 2 score")
  ggplot2::ggsave(fs::path(plot_dir, "phase2_score_vs_detection_frequency.png"), p3, width = 7, height = 5, dpi = 150)
  p4 <- ggplot2::ggplot(phase2_tbl, ggplot2::aes(x = .data$median_qc_rsd, y = .data$phase2_score)) + ggplot2::geom_point(size = 3) + ggplot2::theme_bw() + ggplot2::labs(title = "Phase 2: score vs QC RSD", x = "Median QC RSD", y = "Phase 2 score")
  ggplot2::ggsave(fs::path(plot_dir, "phase2_score_vs_qc_rsd.png"), p4, width = 7, height = 5, dpi = 150)

  summary_file <- fs::path(optim_root, "optimisation_summary.txt")
  summary_lines <- c(
    "MZmine optimisation without benchmark file", "==========================================", "",
    paste("Base batch:", base_batch), paste("Single QC used in Phase 1:", qc_single), paste("Number of QC files used in Phase 2:", length(qc_all)),
    paste("Phase 1 runs:", nrow(param_grid_phase1)), paste("Phase 2 runs:", nrow(top_phase1)), "",
    "Best parameters:", paste("MS1 noise:", best$ms1_noise), paste("Minimum consecutive scans:", best$min_consecutive_scans), paste("Minimum consecutive intensity:", best$min_consecutive_intensity), paste("Chromatogram builder minimum absolute height:", best$chrom_min_absolute_height), paste("Resolver chromatographic threshold:", best$chromatographic_threshold), paste("Resolver minimum height:", best$resolver_min_height), paste("Peak top/edge ratio:", best$peak_top_edge_ratio), paste("Resolver minimum scans:", best$resolver_min_scans), paste("Peak duration min:", best$peak_duration_min), paste("Peak duration max:", best$peak_duration_max), "",
    "Best Phase 2 metrics:", paste("Number of features:", best$n_features), paste("Number of MS2 features:", best$n_ms2_features), paste("MS2 fraction:", best$ms2_fraction), paste("Median height:", best$median_height), paste("Median QC detection frequency:", best$median_detection_frequency), paste("Median QC RSD:", best$median_qc_rsd), paste("Phase 2 score:", best$phase2_score), "",
    paste("Final optimised full batch:", final_optimised_batch), paste("Final QC-only validation batch:", final_qc_only_batch))
  writeLines(summary_lines, summary_file)
  message("\nSummary written to:")
  message(summary_file)
  message("\nDone.")

  invisible(list(best = best, phase1 = phase1_tbl, phase2 = phase2_tbl, final_optimised_batch = final_optimised_batch, final_qc_only_batch = final_qc_only_batch, summary_file = summary_file, plot_dir = plot_dir))
}
