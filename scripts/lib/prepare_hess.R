#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(dirname(normalizePath(script_file)))), "R", "common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: prepare_hess.R config.yml trait_id output.tsv.gz", call. = FALSE)
cfg <- load_project_config(args[[1]])
dat <- read_standard_gwas(project_path(cfg, cfg$traits[[args[[2]]]]$file, TRUE))
out <- dat[, .(SNP, CHR, BP, A1, A2, Z, N = n)]
dir.create(dirname(args[[3]]), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(out, args[[3]], sep = "\t")
