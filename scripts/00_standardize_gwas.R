#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
required_packages(c("data.table", "optparse"))

option_list <- list(
  optparse::make_option("--input", type = "character"),
  optparse::make_option("--output", type = "character"),
  optparse::make_option("--snp", type = "character", default = "SNP"),
  optparse::make_option("--chr", type = "character", default = "CHR"),
  optparse::make_option("--bp", type = "character", default = "BP"),
  optparse::make_option("--a1", type = "character", default = "A1"),
  optparse::make_option("--a2", type = "character", default = "A2"),
  optparse::make_option("--beta", type = "character", default = "beta"),
  optparse::make_option("--se", type = "character", default = "se"),
  optparse::make_option("--p", type = "character", default = "p"),
  optparse::make_option("--eaf", type = "character", default = "eaf"),
  optparse::make_option("--n", type = "character", default = "n"),
  optparse::make_option("--info", type = "character", default = NA_character_),
  optparse::make_option("--maf-min", type = "double", default = 0.01),
  optparse::make_option("--info-min", type = "double", default = 0.90),
  optparse::make_option("--deduplicate", action = "store_true", default = FALSE)
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
if (is.null(opt$input) || is.null(opt$output)) stop("--input and --output are required.", call. = FALSE)

raw <- data.table::fread(opt$input)
mapping <- c(SNP = opt$snp, CHR = opt$chr, BP = opt$bp, A1 = opt$a1, A2 = opt$a2,
             beta = opt$beta, se = opt$se, p = opt$p, eaf = opt$eaf, n = opt$n)
missing <- setdiff(unname(mapping), names(raw))
if (length(missing) > 0L) stop("Source columns not found: ", paste(missing, collapse = ", "), call. = FALSE)

dat <- raw[, ..mapping]
data.table::setnames(dat, names(dat), names(mapping))
numeric_cols <- c("CHR", "BP", "beta", "se", "p", "eaf", "n")
dat[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]
dat[, `:=`(SNP = as.character(SNP), A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)))]
if (!is.na(opt$info)) {
  if (!opt$info %in% names(raw)) stop("INFO source column not found: ", opt$info, call. = FALSE)
  dat[, info := as.numeric(raw[[opt$info]])]
} else {
  dat[, info := NA_real_]
}
dat[, `:=`(Z = beta / se, MAF = pmin(eaf, 1 - eaf))]
dat <- dat[is.finite(Z) & se > 0 & p > 0 & p <= 1 & eaf > 0 & eaf < 1 & MAF >= opt$`maf-min`]
if (!all(is.na(dat$info))) dat <- dat[is.finite(info) & info >= opt$`info-min`]
dat <- dat[CHR %in% 1:22 & BP >= 1 & A1 %in% c("A", "C", "G", "T") & A2 %in% c("A", "C", "G", "T") & A1 != A2]

if (anyDuplicated(dat$SNP)) {
  if (!opt$deduplicate) stop("Duplicated SNP identifiers remain; resolve them or use --deduplicate.", call. = FALSE)
  data.table::setorder(dat, SNP, p)
  dat <- unique(dat, by = "SNP")
}
validate_standard_gwas(dat)
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(dat[, .(SNP, CHR, BP, A1, A2, beta, se, p, eaf, n, info, Z)], opt$output, sep = "\t")
message("Wrote ", nrow(dat), " variants to ", opt$output)
