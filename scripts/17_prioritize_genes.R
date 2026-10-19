#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file))), "R", "common.R"))
required_packages("data.table")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 17_prioritize_genes.R config.yml", call. = FALSE)
cfg <- load_project_config(args[[1]])
root <- output_path(cfg)

fusion_manifest <- data.table::fread(project_path(cfg, cfg$expression$fusion_manifest, TRUE))
fusion_all <- list()
for (trait in c("ra", "mg")) for (tissue in fusion_manifest$tissue) {
  files <- list.files(file.path(root, "fusion_twas", trait, tissue), pattern = "^chr[0-9]+\\.dat$", full.names = TRUE)
  if (!length(files)) next
  x <- data.table::rbindlist(lapply(files, data.table::fread), fill = TRUE)
  if (!all(c("ID", "TWAS.P") %in% names(x))) stop("Unexpected FUSION output columns for ", trait, "/", tissue, call. = FALSE)
  x[, `:=`(trait = trait, tissue = tissue, FDR = stats::p.adjust(TWAS.P, method = "BH"))]
  fusion_all[[length(fusion_all) + 1L]] <- x
}
if (!length(fusion_all)) stop("No FUSION result files were found.", call. = FALSE)
fusion <- data.table::rbindlist(fusion_all, fill = TRUE)
fusion_sig <- fusion[FDR < 0.05]
fusion_shared <- fusion_sig[, if (all(c("ra", "mg") %in% trait)) .SD, by = .(tissue, ID)]

smr_manifest <- data.table::fread(project_path(cfg, cfg$expression$smr_manifest, TRUE))
smr_all <- list()
for (trait in c("ra", "mg")) for (tissue in smr_manifest$tissue) {
  path <- file.path(root, "smr", paste0(trait, "_", tissue, ".smr"))
  if (!file.exists(path)) next
  x <- data.table::fread(path)
  gene_col <- intersect(c("Gene", "Probe"), names(x))[1]
  if (is.na(gene_col) || !all(c("p_SMR", "p_HEIDI") %in% names(x))) stop("Unexpected SMR output columns: ", path, call. = FALSE)
  x[, gene_id := as.character(get(gene_col))]
  x[, `:=`(trait = trait, tissue = tissue, FDR = stats::p.adjust(p_SMR, method = "BH"))]
  smr_all[[length(smr_all) + 1L]] <- x
}
smr <- data.table::rbindlist(smr_all, fill = TRUE)
heidi_n_col <- intersect(c("nsnp_HEIDI", "n_HEIDI"), names(smr))[1]
if (nrow(smr)) {
  smr[, heidi_n := if (!is.na(heidi_n_col)) as.numeric(get(heidi_n_col)) else NA_real_]
  smr_sig <- smr[FDR < 0.05 & p_HEIDI > 0.05 & (is.na(heidi_n) | heidi_n > 10)]
  smr_shared <- smr_sig[, if (all(c("ra", "mg") %in% trait)) .SD, by = .(tissue, gene_id)]
} else {
  smr_sig <- smr_shared <- data.table::data.table()
}

out <- file.path(root, "gene_prioritization"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(fusion, file.path(out, "fusion_all.tsv.gz"), sep = "\t")
data.table::fwrite(fusion_shared, file.path(out, "fusion_shared_fdr05.tsv"), sep = "\t")
data.table::fwrite(smr, file.path(out, "smr_all.tsv.gz"), sep = "\t")
data.table::fwrite(smr_shared, file.path(out, "smr_shared_fdr05_heidi.tsv"), sep = "\t")
