#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
required_packages(c("cause", "data.table", "ieugwasr", "ggplot2"))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 3L)) stop("Usage: 07_cause.R config.yml [exposure_id outcome_id]", call. = FALSE)
cfg <- load_project_config(args[[1]])
pairs <- if (length(args) == 3L) list(args[2:3]) else list(c("ra", "mg"), c("mg", "ra"))

run_cause <- function(exposure_id, outcome_id) {
  exposure <- read_standard_gwas(project_path(cfg, cfg$traits[[exposure_id]]$file, TRUE))
  outcome <- read_standard_gwas(project_path(cfg, cfg$traits[[outcome_id]]$file, TRUE))
  X <- cause::gwas_merge(
    exposure, outcome,
    snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("beta", "beta"),
    se_cols = c("se", "se"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2")
  )
  p_lookup <- exposure[, .(snp = SNP, p_exposure = p)]
  X <- merge(X, p_lookup, by = "snp", all.x = TRUE)
  X <- X[is.finite(p_exposure)]
  if (nrow(X) == 0L) stop("No variants remain for CAUSE.", call. = FALSE)

  set.seed(cfg$project$seed)
  nuisance_snps <- if (nrow(X) <= 1e6) X$snp else sample(X$snp, 1e6, replace = FALSE)
  params <- cause::est_cause_params(X, nuisance_snps)

  pars <- cfg$parameters$cause
  clump_input <- data.frame(rsid = X$snp, pval = X$p_exposure, id = exposure_id)
  clumped <- ieugwasr::ld_clump(
    dat = clump_input, clump_kb = pars$clump_kb, clump_r2 = pars$clump_r2,
    clump_p = pars$p_threshold, bfile = project_path(cfg, cfg$references$plink_eur_prefix),
    plink_bin = tool_path(cfg, cfg$tools$plink)
  )
  variants <- intersect(clumped$rsid, X$snp)
  if (length(variants) < 2L) stop("Too few clumped variants for CAUSE.", call. = FALSE)
  result <- cause::cause(X = X, variants = variants, param_ests = params, force = isTRUE(pars$force))

  out_dir <- output_path(cfg, "cause", paste0(exposure_id, "_to_", outcome_id))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(X, file.path(out_dir, "merged_data.rds"))
  saveRDS(params, file.path(out_dir, "nuisance_parameters.rds"))
  saveRDS(result, file.path(out_dir, "cause_result.rds"))
  writeLines(capture.output(summary(result, ci_size = 0.95)), file.path(out_dir, "summary.txt"))
  grDevices::pdf(file.path(out_dir, "diagnostic_plots.pdf"), width = 8, height = 6)
  print(plot(result$sharing)); print(plot(result$causal)); print(plot(result)); print(plot(result, type = "data"))
  grDevices::dev.off()
}

for (pair in pairs) run_cause(pair[[1]], pair[[2]])
