#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(dirname(normalizePath(script_file)))), "R", "common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: prepare_mixer.R config.yml trait_id output.tsv.gz", call. = FALSE)
cfg <- load_project_config(args[[1]])
trait_id <- args[[2]]
if (!trait_id %in% names(cfg$traits)) stop("Unknown trait: ", trait_id, call. = FALSE)
trait <- cfg$traits[[trait_id]]
dat <- read_standard_gwas(project_path(cfg, trait$file, TRUE))
effective_n <- 4 / (1 / trait$n_case + 1 / trait$n_control)
out <- dat[, .(SNP, CHR, BP, A1, A2, N = effective_n, Z)]
dir.create(dirname(args[[3]]), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(out, args[[3]], sep = "\t", na = "nan")
