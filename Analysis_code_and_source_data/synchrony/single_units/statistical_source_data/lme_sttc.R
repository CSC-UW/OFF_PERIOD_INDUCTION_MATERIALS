############################################
# Dyadic LME (frequentist) with df corrections
# - outcome: dz = atanh(STTC_rebound) - atanh(STTC_late)
# - fixed effect: probe (optrode vs contra)
# - random effects: (1 + probe | mouse)  +  multi-membership (unit_i + unit_j)
# - inference: Satterthwaite df and Kenward-Roger df
############################################

# -------- USER: set your data directory and CSV files ----------
data_dir <- "C:\\Users\\driessen2\\Documents\\sttc_lme_stats"
csv_files <- c(
  "dfsom.csv",
  "dfsom_fix.csv",
  "dfacr.csv",
  "dfacr_fix.csv",
  "dfhalo.csv",
  "dfhalo_fix.csv"
)

# -------- USER: set column names if needed -------
mouse_col       <- "mouse"
probe_col       <- "probe"
unit_i_col      <- "unit_i"
unit_j_col      <- "unit_j"
sttc_late_col   <- "sttc_late"
sttc_reb_col    <- "sttc_rebound"

contra_label    <- "contra"
optrode_label   <- "optrode"

# -------- packages (install if needed) ----------
pkgs_cran <- c("data.table", "dplyr", "Matrix", "lme4", "lmerTest", "pbkrtest")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

# lmerMultiMember is on GitHub; install once if you don't have it
if (!requireNamespace("lmerMultiMember", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("jvparidon/lmerMultiMember")
}

library(data.table)
library(dplyr)
library(Matrix)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(lmerMultiMember)

# -------- helper: fisher z with clipping ----------
fisher_z <- function(r, eps = 1e-6) {
  r <- pmin(pmax(r, -1 + eps), 1 - eps)
  atanh(r)
}

# -------- main analysis function ----------
run_lme_analysis <- function(csv_path) {

  cat("\n\n########################################################\n")
  cat("# Starting analysis for:", basename(csv_path), "\n")
  cat("########################################################\n\n")

  # -------- load data ----------
  df <- data.table::fread(csv_path) %>% as.data.frame()

  # -------- basic validation ----------
  required_cols <- c(mouse_col, probe_col, unit_i_col, unit_j_col, sttc_late_col, sttc_reb_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # -------- compute dz (paired difference) ----------
  df <- df %>%
    mutate(
      sttc_late = as.numeric(.data[[sttc_late_col]]),
      sttc_reb  = as.numeric(.data[[sttc_reb_col]]),
      dz = fisher_z(sttc_reb) - fisher_z(sttc_late)
    )

  # -------- encode probe indicator (0=contra, 1=optrode) ----------
  df <- df %>%
    mutate(
      mouse = factor(.data[[mouse_col]]),
      probe = factor(.data[[probe_col]])
    )

  bad_probe <- setdiff(levels(df$probe), c(contra_label, optrode_label))
  if (length(bad_probe) > 0) {
    stop(paste("Unexpected probe labels present:", paste(bad_probe, collapse = ", ")))
  }

  df <- df %>%
    mutate(
      is_optrode = ifelse(as.character(probe) == optrode_label, 1, 0)
    )

  # -------- make globally-unique unit IDs ----------
  df <- df %>%
    mutate(
      unit_i_id = paste(as.character(mouse), as.character(probe), as.character(.data[[unit_i_col]]), sep = "||"),
      unit_j_id = paste(as.character(mouse), as.character(probe), as.character(.data[[unit_j_col]]), sep = "||")
    )

  # Ensure no missing values in the variables used by the model
  df_cc <- df %>% dplyr::filter(complete.cases(dz, is_optrode, mouse, unit_i_id, unit_j_id))

  # Rebuild the membership matrix using df_cc
  unit_levels <- sort(unique(c(df_cc$unit_i_id, df_cc$unit_j_id)))
  n_obs  <- nrow(df_cc)
  n_unit <- length(unit_levels)

  i_idx <- match(df_cc$unit_i_id, unit_levels)
  j_idx <- match(df_cc$unit_j_id, unit_levels)

  row_index <- rep(seq_len(n_obs), times = 2)
  col_index <- c(i_idx, j_idx)
  x_vals    <- rep(1.0, times = 2 * n_obs)

  W_unit_obs_by_unit <- Matrix::sparseMatrix(
    i = row_index, j = col_index, x = x_vals,
    dims = c(n_obs, n_unit),
    dimnames = list(NULL, unit_levels)
  )

  # lmerMultiMember expects: (units x observations)
  W_unit <- Matrix::t(W_unit_obs_by_unit)

  # Sanity checks
  stopifnot(ncol(W_unit) == nrow(df_cc))
  stopifnot(all(abs(Matrix::colSums(W_unit) - 2) < 1e-8))

  # -------- fit the dyadic mixed model ----------
  mod_full <- lmerMultiMember::lmer(
    dz ~ is_optrode + (1 + is_optrode | mouse) + (1 | unit),
    data = df_cc,
    memberships = list(unit = W_unit),
    REML = TRUE
  )

  cat("\n=== Model fit summary ===\n")
  print(summary(mod_full))

  cat("\nSingular fit? (variance components on boundary)\n")
  print(isSingular(mod_full, tol = 1e-5))

  # -------- inference: Satterthwaite ----------
  anova_satt <- anova(mod_full, ddf = "Satterthwaite")
  cat("\n=== Fixed effect test (Satterthwaite) ===\n")
  print(anova_satt)

  # -------- inference: Kenward-Roger ----------
  anova_kr <- tryCatch(
    anova(mod_full, ddf = "Kenward-Roger"),
    error = function(e) {
      cat("\n  KR test failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(anova_kr)) {
    cat("\n=== Fixed effect test (Kenward-Roger) ===\n")
    print(anova_kr)
  }

  # -------- extract key numbers for reporting ----------
  # beta, SE, t from summary coefficients
  coef_tab <- summary(mod_full)$coefficients
  beta_hat <- coef_tab["is_optrode", "Estimate"]
  se_hat   <- coef_tab["is_optrode", "Std. Error"]
  tval     <- coef_tab["is_optrode", "t value"]

  # Satterthwaite df and p from anova (lmerMultiMember doesn't produce
  # the augmented coefficient table with df and Pr(>|t|) columns)
  satt_df <- anova_satt["is_optrode", "DenDF"]
  satt_F  <- anova_satt["is_optrode", "F value"]
  satt_p  <- anova_satt["is_optrode", "Pr(>F)"]

  # Kenward-Roger df and p from anova (if available)
  if (!is.null(anova_kr)) {
    kr_df <- anova_kr["is_optrode", "DenDF"]
    kr_F  <- anova_kr["is_optrode", "F value"]
    kr_p  <- anova_kr["is_optrode", "Pr(>F)"]
  } else {
    kr_df <- NA
    kr_F  <- NA
    kr_p  <- NA
  }

  # Effect size: Cohen's f^2 = F / DenDF (NumDF = 1)
  # Benchmarks: 0.02 = small, 0.15 = medium, 0.35 = large (Cohen, 1988)
  kr_f2 <- if (!is.na(kr_F) && !is.na(kr_df)) kr_F / kr_df else NA

  cat("\n=== PRIMARY RESULT (Satterthwaite, two-tailed) ===\n")
  cat(sprintf("  beta = %.6f, SE = %.6f, t(%.1f) = %.4f, p = %.6f\n",
              beta_hat, se_hat, satt_df, tval, satt_p))

  if (!is.null(anova_kr)) {
    cat(sprintf("  Kenward-Roger: F(1, %.1f) = %.4f, p = %.6f, Cohen's f2 = %.4f\n",
                kr_df, kr_F, kr_p, kr_f2))
  }

  # -------- data summary ----------
  cat("\n=== Data summary ===\n")
  cat(sprintf("  N observations (unit pairs): %d\n", nrow(df_cc)))
  cat(sprintf("  N mice: %d\n", length(unique(df_cc$mouse))))
  cat(sprintf("  N unique units: %d\n", n_unit))

  probe_summary <- df_cc %>%
    group_by(probe) %>%
    summarise(
      n_pairs = n(),
      n_units = length(unique(c(unit_i_id, unit_j_id))),
      mean_dz = mean(dz, na.rm = TRUE),
      sd_dz   = sd(dz, na.rm = TRUE),
      .groups = "drop"
    )
  cat("\n  Per-probe summary:\n")
  print(as.data.frame(probe_summary))

  mouse_probe_summary <- df_cc %>%
    group_by(mouse, probe) %>%
    summarise(n_pairs = n(), .groups = "drop")
  cat("\n  Pairs per mouse x probe:\n")
  print(as.data.frame(mouse_probe_summary))

  # -------- save all outputs ----------
  csv_basename <- tools::file_path_sans_ext(basename(csv_path))
  out_dir      <- dirname(csv_path)
  out_txt      <- file.path(out_dir, paste0(csv_basename, "_lme_results.txt"))
  out_rdata    <- file.path(out_dir, paste0(csv_basename, "_lme_results.RData"))

  sink(out_txt)
  cat("============================================\n")
  cat("Dyadic LME results for:", csv_path, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("============================================\n")

  cat("\n=== Data summary ===\n")
  cat(sprintf("  N observations (unit pairs): %d\n", nrow(df_cc)))
  cat(sprintf("  N mice: %d\n", length(unique(df_cc$mouse))))
  cat(sprintf("  N unique units: %d\n", n_unit))
  cat("\n  Per-probe summary:\n")
  print(as.data.frame(probe_summary))
  cat("\n  Pairs per mouse x probe:\n")
  print(as.data.frame(mouse_probe_summary))

  cat("\n=== Model fit summary ===\n")
  print(summary(mod_full))

  cat("\nSingular fit?\n")
  print(isSingular(mod_full, tol = 1e-5))

  cat("\n=== Variance components ===\n")
  print(VarCorr(mod_full))

  cat("\n=== Fixed effects (from summary) ===\n")
  print(coef_tab)

  cat("\n=== PRIMARY RESULT (Satterthwaite, two-tailed) ===\n")
  cat(sprintf("  beta = %.6f, SE = %.6f, t(%.1f) = %.4f, p = %.6f\n",
              beta_hat, se_hat, satt_df, tval, satt_p))

  cat("\n=== Fixed effect test (Satterthwaite ANOVA) ===\n")
  print(anova_satt)

  cat("\n=== Fixed effect test (Kenward-Roger ANOVA) ===\n")
  if (!is.null(anova_kr)) {
    print(anova_kr)
    cat(sprintf("\n  Kenward-Roger: F(1, %.1f) = %.4f, p = %.6f, Cohen's f2 = %.4f\n",
                kr_df, kr_F, kr_p, kr_f2))
  } else {
    cat("  KR test failed (see console output for details)\n")
  }

  cat("\n=== Random effects (BLUPs by mouse) ===\n")
  print(ranef(mod_full)$mouse)

  cat("\n=== Session info ===\n")
  print(sessionInfo())
  sink()

  cat(sprintf("\nResults saved to: %s\n", out_txt))

  # Save structured R objects for programmatic access
  results <- list(
    csv_path         = csv_path,
    n_obs            = nrow(df_cc),
    n_mice           = length(unique(df_cc$mouse)),
    n_units          = n_unit,
    probe_summary    = probe_summary,
    model_full       = mod_full,
    coef_table       = coef_tab,
    anova_satt       = anova_satt,
    anova_kr         = anova_kr,
    beta_hat         = beta_hat,
    se_hat           = se_hat,
    t_value          = tval,
    satt_df          = satt_df,
    satt_F           = satt_F,
    satt_p           = satt_p,
    kr_df            = kr_df,
    kr_F             = kr_F,
    kr_p             = kr_p,
    kr_f2            = kr_f2,
    is_singular      = isSingular(mod_full, tol = 1e-5),
    var_corr         = VarCorr(mod_full),
    session_info     = sessionInfo()
  )
  save(results, file = out_rdata)
  cat(sprintf("R objects saved to: %s\n", out_rdata))

  return(results)
}

# -------- run analysis for all CSV files ----------
all_results <- list()
for (csv_file in csv_files) {
  csv_path <- file.path(data_dir, csv_file)
  if (!file.exists(csv_path)) {
    cat(sprintf("\nWARNING: %s not found, skipping.\n", csv_path))
    next
  }
  all_results[[csv_file]] <- tryCatch(
    run_lme_analysis(csv_path),
    error = function(e) {
      cat(sprintf("\nERROR processing %s: %s\n", csv_file, conditionMessage(e)))
      NULL
    }
  )
}

# -------- print summary across all datasets ----------
cat("\n\n########################################################\n")
cat("# SUMMARY ACROSS ALL DATASETS (Satterthwaite)\n")
cat("########################################################\n\n")
cat(sprintf("%-20s  %10s  %10s  %10s  %10s  %s\n",
            "Dataset", "beta", "t", "df", "p", "Notes"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (name in names(all_results)) {
  r <- all_results[[name]]
  if (!is.null(r)) {
    notes <- ifelse(r$is_singular, "[SINGULAR]", "")
    cat(sprintf("%-20s  %+10.6f  %+10.4f  %10.1f  %10.6f  %s\n",
                name, r$beta_hat, r$t_value, r$satt_df, r$satt_p, notes))
  } else {
    cat(sprintf("%-20s  FAILED\n", name))
  }
}

cat("\n")
cat(sprintf("%-20s  %10s  %10s  %10s  %10s  %10s  %s\n",
            "Dataset", "F", "NumDF", "DenDF", "p(KR)", "f2", "Notes"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (name in names(all_results)) {
  r <- all_results[[name]]
  if (!is.null(r)) {
    notes <- ifelse(r$is_singular, "[SINGULAR]", "")
    if (!is.na(r$kr_p)) {
      cat(sprintf("%-20s  %10.4f  %10d  %10.1f  %10.6f  %10.4f  %s\n",
                  name, r$kr_F, 1L, r$kr_df, r$kr_p, r$kr_f2, notes))
    } else {
      cat(sprintf("%-20s  KR failed  %s\n", name, notes))
    }
  } else {
    cat(sprintf("%-20s  FAILED\n", name))
  }
}

# -------- save summary CSV across all datasets ----------
summary_rows <- list()
for (name in names(all_results)) {
  r <- all_results[[name]]
  if (!is.null(r)) {
    # Parse mouse line and pair type from filename
    # e.g. "dfsom.csv" -> SOM, homeostatic; "dfsom_fix.csv" -> SOM, fixed
    bn <- tools::file_path_sans_ext(name)
    is_fix <- grepl("_fix$", bn)
    line_raw <- sub("^df", "", sub("_fix$", "", bn))
    mouse_line <- toupper(line_raw)
    pair_type <- ifelse(is_fix, "fixed", "homeostatic")

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Dataset        = name,
      Mouse_Line     = mouse_line,
      Pair_Type      = pair_type,
      N_Pairs        = r$n_obs,
      N_Mice         = r$n_mice,
      N_Units        = r$n_units,
      Beta           = r$beta_hat,
      SE             = r$se_hat,
      Satt_t         = r$t_value,
      Satt_df        = r$satt_df,
      Satt_p         = r$satt_p,
      KR_F           = r$kr_F,
      KR_df          = r$kr_df,
      KR_p           = r$kr_p,
      Cohens_f2      = r$kr_f2,
      Singular_Fit   = r$is_singular,
      stringsAsFactors = FALSE
    )
  }
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  summary_csv_path <- file.path(data_dir, "sttc_lme_summary_table.csv")
  write.csv(summary_df, summary_csv_path, row.names = FALSE)
  cat(sprintf("\n=== Summary table saved to: %s ===\n", summary_csv_path))
} else {
  cat("\nNo successful results to save.\n")
}

cat("\n=== Session info ===\n")
print(sessionInfo())
