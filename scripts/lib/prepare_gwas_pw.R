#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(dirname(normalizePath(script_file)))), "R", "common.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare_gwas_pw.R config.yml output.txt.gz", call. = FALSE)
cfg <- load_project_config(args[[1]])
ra <- read_standard_gwas(project_path(cfg, cfg$traits$ra$file, TRUE))
mg <- read_standard_gwas(project_path(cfg, cfg$traits$mg$file, TRUE))
dat <- harmonize_pair(ra, mg, c(".RA", ".MG"), drop_palindromic = TRUE)
dat <- dat[CHR.RA == CHR.MG & BP.RA == BP.MG]
dat[allele_flip, `:=`(beta.MG = -beta.MG, Z.MG = -Z.MG, eaf.MG = 1 - eaf.MG)]
out <- dat[, .(SNPID = SNP, CHR = CHR.RA, POS = BP.RA, Z_RA = Z.RA,
               V_RA = se.RA^2, Z_MG = Z.MG, V_MG = se.MG^2)]
data.table::setorder(out, CHR, POS)
dir.create(dirname(args[[2]]), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(out, args[[2]], sep = "\t")
