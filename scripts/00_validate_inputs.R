#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: 00_validate_inputs.R config/config.yml [--strict-references]", call. = FALSE)

cfg <- load_project_config(args[[1]])
for (trait_id in names(cfg$traits)) {
  trait <- cfg$traits[[trait_id]]
  path <- project_path(cfg, trait$file, must_work = TRUE)
  dat <- read_standard_gwas(path)
  validate_standard_gwas(dat, trait_id)
  expected_n <- trait$n_case + trait$n_control
  message(trait_id, ": ", nrow(dat), " variants; configured total N = ", expected_n)
}

if ("--strict-references" %in% args) {
  refs <- unlist(cfg$references, use.names = TRUE)
  plain_refs <- refs[!grepl("@", refs, fixed = TRUE)]
  missing <- names(plain_refs)[!file.exists(vapply(plain_refs, function(x) project_path(cfg, x), character(1)))]
  if (length(missing) > 0L) stop("Missing reference paths: ", paste(missing, collapse = ", "), call. = FALSE)
}
message("Input validation passed.")
