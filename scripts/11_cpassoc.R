#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
required_packages(c("data.table", "ieugwasr"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 11_cpassoc.R config.yml", call. = FALSE)
cfg <- load_project_config(args[[1]])
source(project_path(cfg, cfg$tools$cpassoc_functions, TRUE))

ra <- read_standard_gwas(project_path(cfg, cfg$traits$ra$file, TRUE))
mg <- read_standard_gwas(project_path(cfg, cfg$traits$mg$file, TRUE))
dat <- harmonize_pair(ra, mg, c(".RA", ".MG"), drop_palindromic = TRUE)
dat <- dat[CHR.RA == CHR.MG & BP.RA == BP.MG]
dat[allele_flip, `:=`(beta.MG = -beta.MG, Z.MG = -Z.MG, eaf.MG = 1 - eaf.MG)]
sumstat <- as.matrix(dat[, .(Z.RA, Z.MG)])

pruned <- data.table::fread(project_path(cfg, cfg$references$cpassoc_pruned_snps, TRUE), header = FALSE)[[1]]
corr_dat <- dat[SNP %in% pruned & abs(Z.RA) <= 1.96 & abs(Z.MG) <= 1.96]
if (nrow(corr_dat) < 20L) stop("Too few null, LD-pruned SNPs to estimate the CPASSOC correlation matrix.", call. = FALSE)
cor_matrix <- diag(2)
cor_matrix[1, 2] <- cor_matrix[2, 1] <- stats::cor(corr_dat$Z.RA, corr_dat$Z.MG, use = "complete.obs")
sample_size <- c(cfg$traits$ra$n_case + cfg$traits$ra$n_control,
                 cfg$traits$mg$n_case + cfg$traits$mg$n_control)

test_shom <- SHom(sumstat, sample_size, cor_matrix)
p_shom <- stats::pchisq(test_shom, df = 1, lower.tail = FALSE)
pars <- cfg$parameters$cpassoc
folds <- cut(seq_len(nrow(sumstat)), breaks = pars$folds, labels = FALSE)
gamma_par <- EstimateGamma(N = pars$gamma_simulations, sample_size, cor_matrix)
p_shet <- rep(NA_real_, nrow(sumstat))
for (i in seq_len(pars$folds)) {
  idx <- which(folds == i)
  test_shet <- SHet(sumstat[idx, , drop = FALSE], sample_size, cor_matrix)
  p_shet[idx] <- stats::pgamma(test_shet - gamma_par[3], shape = gamma_par[1], scale = gamma_par[2], lower.tail = FALSE)
}
results <- dat[, .(SNP, CHR = CHR.RA, BP = BP.RA, p_RA = p.RA, p_MG = p.MG)]
results[, `:=`(p_SHom = p_shom, p_SHet = p_shet)]
significant <- results[p_SHet < pars$significance & p_RA < pars$single_trait_p & p_MG < pars$single_trait_p]

if (nrow(significant)) {
  clump_input <- data.frame(rsid = significant$SNP, pval = significant$p_SHet, id = "RA_MG")
  clumped <- ieugwasr::ld_clump(clump_input, clump_kb = 500, clump_r2 = 0.2,
    clump_p = pars$significance, bfile = project_path(cfg, cfg$references$plink_eur_prefix),
    plink_bin = tool_path(cfg, cfg$tools$plink))
  significant <- significant[SNP %in% clumped$rsid]
}
out <- output_path(cfg, "cpassoc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(results, file.path(out, "all_results.tsv.gz"), sep = "\t")
data.table::fwrite(significant, file.path(out, "significant_clumped_loci.tsv"), sep = "\t")
data.table::fwrite(data.table::as.data.table(cor_matrix, keep.rownames = "trait"), file.path(out, "correlation_matrix.tsv"), sep = "\t")
