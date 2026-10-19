#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
required_packages(c("data.table", "ieugwasr"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 10b_mtag_postprocess.R config.yml trait1_mtag.txt trait2_mtag.txt", call. = FALSE)
cfg <- load_project_config(args[[1]])

read_mtag <- function(path, label) {
  x <- data.table::fread(path)
  required <- c("snp_name", "mtag_beta", "mtag_pval")
  if (!all(required %in% names(x))) stop("MTAG output missing columns: ", paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
  x <- x[, .(SNP = as.character(snp_name), mtag_beta = as.numeric(mtag_beta), mtag_p = as.numeric(mtag_pval))]
  data.table::setnames(x, c("mtag_beta", "mtag_p"), paste0(c("mtag_beta_", "mtag_p_"), label))
  x
}

m1 <- read_mtag(args[[2]], "ra")
m2 <- read_mtag(args[[3]], "mg")
mtag <- merge(m1, m2, by = "SNP")
mtag[, joint_mtag_p := pmin(mtag_p_ra, mtag_p_mg)]
mtag <- mtag[joint_mtag_p < 5e-8]

ra <- read_standard_gwas(project_path(cfg, cfg$traits$ra$file, TRUE))
mg <- read_standard_gwas(project_path(cfg, cfg$traits$mg$file, TRUE))
original <- harmonize_pair(ra, mg, c(".RA", ".MG"), drop_palindromic = TRUE)
original <- original[CHR.RA == CHR.MG & BP.RA == BP.MG]
original[allele_flip, `:=`(beta.MG = -beta.MG, Z.MG = -Z.MG, eaf.MG = 1 - eaf.MG)]
dat <- merge(mtag, original[, .(SNP, CHR = CHR.RA, BP = BP.RA, beta_RA = beta.RA,
  beta_MG = beta.MG, p_RA = p.RA, p_MG = p.MG)], by = "SNP")
dat <- dat[p_RA < 0.05 & p_MG < 0.05 & sign(beta_RA) == sign(beta_MG)]

if (nrow(dat)) {
  clump_input <- data.frame(rsid = dat$SNP, pval = dat$joint_mtag_p, id = "RA_MG")
  clumped <- ieugwasr::ld_clump(clump_input, clump_kb = 500, clump_r2 = 0.2,
    clump_p = 5e-8, bfile = project_path(cfg, cfg$references$plink_eur_prefix),
    plink_bin = tool_path(cfg, cfg$tools$plink))
  dat <- dat[SNP %in% clumped$rsid]
}
out <- output_path(cfg, "mtag", "significant_clumped_loci.tsv")
data.table::fwrite(dat, out, sep = "\t")
