#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
required_packages(c("MAGMA.Celltyping", "MungeSumstats", "data.table"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: 09_magma_celltyping.R config.yml", call. = FALSE)
cfg <- load_project_config(args[[1]])
ctd <- readRDS(project_path(cfg, cfg$expression$ctd_rds, TRUE))
pars <- cfg$parameters$magma

for (trait_id in c("ra", "mg")) {
  trait <- cfg$traits[[trait_id]]
  out_dir <- output_path(cfg, "magma_celltyping", trait_id)
  storage <- file.path(out_dir, "cache")
  dir.create(storage, recursive = TRUE, showWarnings = FALSE)
  formatted <- file.path(out_dir, paste0(trait_id, ".formatted.tsv.gz"))
  MungeSumstats::format_sumstats(
    path = project_path(cfg, trait$file, TRUE), save_path = formatted,
    ref_genome = "GRCh37", return_data = FALSE
  )
  genes_out <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = formatted, genome_build = "GRCh37",
    upstream_kb = pars$upstream_kb, downstream_kb = pars$downstream_kb,
    N = trait$n_case + trait$n_control, population = "eur", storage_dir = storage
  )
  result <- MAGMA.Celltyping::celltype_associations_pipeline(
    magma_dirs = dirname(genes_out), ctd = ctd, ctd_species = "human",
    ctd_name = "study_ctd", run_linear = isTRUE(pars$run_linear),
    run_top10 = isTRUE(pars$run_top10), save_dir = out_dir
  )
  merged <- MAGMA.Celltyping::merge_results(MAGMA_results = result, save_dir = out_dir)
  data.table::fwrite(data.table::as.data.table(merged), file.path(out_dir, "merged_results.tsv"), sep = "\t")
  saveRDS(result, file.path(out_dir, "magma_celltyping_result.rds"))
}
