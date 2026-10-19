#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
required_packages(c("data.table", "TwoSampleMR", "ieugwasr", "ggplot2"))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 3L)) {
  stop("Usage: 06_two_sample_mr.R config.yml [exposure_id outcome_id]", call. = FALSE)
}
cfg <- load_project_config(args[[1]])
pairs <- if (length(args) == 3L) list(args[2:3]) else list(c("ra", "mg"), c("mg", "ra"))

run_pair <- function(exposure_id, outcome_id) {
  if (!all(c(exposure_id, outcome_id) %in% names(cfg$traits))) stop("Unknown trait ID.", call. = FALSE)
  pars <- cfg$parameters$mr
  exp_cfg <- cfg$traits[[exposure_id]]
  out_cfg <- cfg$traits[[outcome_id]]
  exposure <- read_standard_gwas(project_path(cfg, exp_cfg$file, TRUE))[p < pars$p_threshold]
  exposure[, F_statistic := (beta / se)^2]
  exposure <- exposure[F_statistic >= pars$f_min]
  if (nrow(exposure) == 0L) stop("No strong exposure instruments for ", exposure_id, call. = FALSE)

  clump_input <- data.frame(rsid = exposure$SNP, pval = exposure$p, id = exposure_id)
  clumped <- ieugwasr::ld_clump(
    dat = clump_input,
    clump_kb = pars$clump_kb,
    clump_r2 = pars$clump_r2,
    clump_p = pars$p_threshold,
    bfile = project_path(cfg, cfg$references$plink_eur_prefix),
    plink_bin = tool_path(cfg, cfg$tools$plink)
  )
  exposure <- exposure[SNP %in% clumped$rsid]
  outcome <- read_standard_gwas(project_path(cfg, out_cfg$file, TRUE))[SNP %in% exposure$SNP]

  exp_dat <- TwoSampleMR::format_data(
    as.data.frame(exposure), type = "exposure", phenotype_col = "phenotype",
    snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p",
    samplesize_col = "n"
  )
  exp_dat$exposure <- exp_cfg$label
  out_dat <- TwoSampleMR::format_data(
    as.data.frame(outcome), type = "outcome", phenotype_col = "phenotype",
    snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p",
    samplesize_col = "n"
  )
  out_dat$outcome <- out_cfg$label
  dat <- TwoSampleMR::harmonise_data(exp_dat, out_dat, action = pars$harmonise_action)
  dat <- dat[dat$mr_keep, , drop = FALSE]
  if (nrow(dat) == 0L) stop("No instruments remain after harmonization.", call. = FALSE)

  steiger <- tryCatch(TwoSampleMR::steiger_filtering(dat), error = function(e) NULL)
  if (!is.null(steiger)) dat <- steiger[is.na(steiger$steiger_dir) | steiger$steiger_dir, , drop = FALSE]
  if (nrow(dat) == 0L) stop("No instruments remain after directionality filtering.", call. = FALSE)

  fit <- TwoSampleMR::mr(dat)
  heterogeneity <- if (nrow(dat) > 1L) TwoSampleMR::mr_heterogeneity(dat) else data.frame()
  pleiotropy <- if (nrow(dat) > 2L) TwoSampleMR::mr_pleiotropy_test(dat) else data.frame()
  single_snp <- TwoSampleMR::mr_singlesnp(dat)
  leave_one_out <- if (nrow(dat) > 2L) TwoSampleMR::mr_leaveoneout(dat) else data.frame()

  presso <- NULL
  outlier_indices <- integer()
  if (nrow(dat) >= 4L && requireNamespace("MRPRESSO", quietly = TRUE)) {
    set.seed(cfg$project$seed)
    presso <- tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = dat, NbDistribution = pars$presso_distributions,
        SignifThreshold = 0.05
      ),
      error = function(e) structure(list(message = conditionMessage(e)), class = "mr_presso_error")
    )
    if (!inherits(presso, "mr_presso_error")) {
      candidate <- presso[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]
      candidate <- suppressWarnings(as.integer(candidate))
      candidate <- candidate[!is.na(candidate)]
      outlier_indices <- unique(candidate[candidate >= 1 & candidate <= nrow(dat)])
    }
  }
  corrected_dat <- if (length(outlier_indices)) dat[-outlier_indices, , drop = FALSE] else dat
  corrected_fit <- if (nrow(corrected_dat) > 0L) TwoSampleMR::mr(corrected_dat) else data.frame()

  out_dir <- output_path(cfg, "mr", paste0(exposure_id, "_to_", outcome_id))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(dat, file.path(out_dir, "harmonized_instruments.tsv"), sep = "\t")
  data.table::fwrite(data.table::as.data.table(TwoSampleMR::generate_odds_ratios(fit)), file.path(out_dir, "mr_results.tsv"), sep = "\t")
  data.table::fwrite(data.table::as.data.table(TwoSampleMR::generate_odds_ratios(corrected_fit)), file.path(out_dir, "mr_results_presso_corrected.tsv"), sep = "\t")
  data.table::fwrite(data.table::as.data.table(heterogeneity), file.path(out_dir, "heterogeneity.tsv"), sep = "\t")
  data.table::fwrite(data.table::as.data.table(pleiotropy), file.path(out_dir, "egger_intercept.tsv"), sep = "\t")
  data.table::fwrite(data.table::as.data.table(single_snp), file.path(out_dir, "single_snp.tsv"), sep = "\t")
  if (nrow(leave_one_out)) data.table::fwrite(leave_one_out, file.path(out_dir, "leave_one_out.tsv"), sep = "\t")
  saveRDS(list(raw = presso, outlier_indices = outlier_indices), file.path(out_dir, "mr_presso.rds"))

  plot_and_save <- function(plot_list, filename) {
    if (length(plot_list)) try(ggplot2::ggsave(file.path(out_dir, filename), plot_list[[1]], width = 7, height = 7), silent = TRUE)
  }
  plot_and_save(TwoSampleMR::mr_scatter_plot(fit, dat), "scatter.pdf")
  plot_and_save(TwoSampleMR::mr_forest_plot(single_snp), "forest.pdf")
  plot_and_save(TwoSampleMR::mr_funnel_plot(single_snp), "funnel.pdf")
  if (nrow(leave_one_out)) plot_and_save(TwoSampleMR::mr_leaveoneout_plot(leave_one_out), "leave_one_out.pdf")
}

for (pair in pairs) run_pair(pair[[1]], pair[[2]])
