#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
required_packages(c("data.table", "coloc"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 12_coloc.R config.yml loci.tsv", call. = FALSE)
cfg <- load_project_config(args[[1]])
loci <- data.table::fread(args[[2]])
needed <- c("lead_snp", "CHR", "BP")
if (!all(needed %in% names(loci))) stop("Locus file requires: ", paste(needed, collapse = ", "), call. = FALSE)

ra <- read_standard_gwas(project_path(cfg, cfg$traits$ra$file, TRUE))
mg <- read_standard_gwas(project_path(cfg, cfg$traits$mg$file, TRUE))
pars <- cfg$parameters$coloc
out <- output_path(cfg, "coloc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
summaries <- list()

for (i in seq_len(nrow(loci))) {
  locus <- loci[i]
  ra_locus <- ra[CHR == locus$CHR & BP >= locus$BP - pars$window_bp & BP <= locus$BP + pars$window_bp]
  mg_locus <- mg[CHR == locus$CHR & BP >= locus$BP - pars$window_bp & BP <= locus$BP + pars$window_bp]
  dat <- harmonize_pair(ra_locus, mg_locus, c(".RA", ".MG"), drop_palindromic = TRUE)
  dat <- dat[CHR.RA == CHR.MG & BP.RA == BP.MG]
  dat[allele_flip, `:=`(beta.MG = -beta.MG, Z.MG = -Z.MG, eaf.MG = 1 - eaf.MG)]
  dat <- dat[is.finite(beta.RA) & is.finite(se.RA) & se.RA > 0 & is.finite(eaf.RA) & eaf.RA > 0 & eaf.RA < 1 &
             is.finite(beta.MG) & is.finite(se.MG) & se.MG > 0 & is.finite(eaf.MG) & eaf.MG > 0 & eaf.MG < 1]
  if (nrow(dat) < 50L) warning("Only ", nrow(dat), " harmonized variants at ", locus$lead_snp)
  if (nrow(dat) == 0L) next

  d_ra <- list(snp = dat$SNP, beta = dat$beta.RA, varbeta = dat$se.RA^2,
    MAF = pmin(dat$eaf.RA, 1 - dat$eaf.RA), N = cfg$traits$ra$n_case + cfg$traits$ra$n_control,
    s = case_fraction(cfg$traits$ra), type = "cc", position = dat$BP.RA)
  d_mg <- list(snp = dat$SNP, beta = dat$beta.MG, varbeta = dat$se.MG^2,
    MAF = pmin(dat$eaf.MG, 1 - dat$eaf.MG), N = cfg$traits$mg$n_case + cfg$traits$mg$n_control,
    s = case_fraction(cfg$traits$mg), type = "cc", position = dat$BP.MG)
  coloc::check_dataset(d_ra); coloc::check_dataset(d_mg)
  result <- coloc::coloc.abf(d_ra, d_mg, p1 = pars$p1, p2 = pars$p2, p12 = pars$p12)
  prefix <- file.path(out, locus$lead_snp)
  data.table::fwrite(dat, paste0(prefix, "_harmonized.tsv.gz"), sep = "\t")
  per_snp <- data.table::as.data.table(result$results)
  per_snp[, lead_snp := locus$lead_snp]
  data.table::fwrite(per_snp, paste0(prefix, "_snp_results.tsv.gz"), sep = "\t")
  summary <- data.table::as.data.table(as.list(result$summary))
  summary[, `:=`(lead_snp = locus$lead_snp, CHR = locus$CHR, BP = locus$BP,
                  window_bp_each_side = pars$window_bp, n_overlap_snps = nrow(dat),
                  coloc_supported = get("PP.H4.abf") > pars$pp_h4_threshold)]
  summaries[[length(summaries) + 1L]] <- summary
}
data.table::fwrite(data.table::rbindlist(summaries, fill = TRUE), file.path(out, "coloc_summary.tsv"), sep = "\t")
