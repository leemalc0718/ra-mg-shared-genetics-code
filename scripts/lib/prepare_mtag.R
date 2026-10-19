#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(dirname(normalizePath(script_file)))), "R", "common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: prepare_mtag.R config.yml trait_id output.tsv.gz", call. = FALSE)
cfg <- load_project_config(args[[1]])
trait <- cfg$traits[[args[[2]]]]
dat <- read_standard_gwas(project_path(cfg, trait$file, TRUE))
fraction <- cfg$parameters$mtag$n_min_fraction_of_p90
n_threshold <- fraction * as.numeric(stats::quantile(dat$n, 0.90, na.rm = TRUE))
maf_min <- cfg$parameters$mtag$maf_min
dat <- dat[n >= n_threshold & pmin(eaf, 1 - eaf) >= maf_min]
dir.create(dirname(args[[3]]), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(dat, args[[3]], sep = "\t")
