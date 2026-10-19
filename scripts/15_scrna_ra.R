#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
root <- dirname(dirname(normalizePath(script_file)))
source(file.path(root, "R", "common.R")); source(file.path(root, "R", "scrna_common.R"))
required_packages(c("Seurat", "SeuratObject", "data.table", "ggplot2", "patchwork"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 15_scrna_ra.R config.yml", call. = FALSE)
cfg <- load_project_config(args[[1]])
sc <- cfg$single_cell$ra
input_dir <- project_path(cfg, sc$input_dir, TRUE)

objects <- lapply(sc$sample_names, function(sample_name) {
  counts <- Seurat::Read10X(file.path(input_dir, sample_name), gene.column = sc$gene_column, unique.features = TRUE)
  object <- Seurat::CreateSeuratObject(counts, project = sample_name, min.cells = 3, min.features = 200)
  object$sample <- sample_name
  object$group <- "RA"
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
  subset(object, subset = nFeature_RNA > sc$min_features & nFeature_RNA < sc$max_features & percent.mt < sc$max_percent_mt)
})
names(objects) <- sc$sample_names
object <- if (length(objects) == 1L) objects[[1]] else merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects))
out <- output_path(cfg, "scrna", "ra"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
object <- run_scrna_workflow(object, sc$dims, cfg$parameters$scrna$variable_features,
  cfg$parameters$scrna$clustering_resolution, cfg$project$seed, out)
object <- annotate_clusters(object, project_path(cfg, sc$annotation_file, TRUE))
save_scrna_outputs(object, unlist(cfg$parameters$scrna$candidate_genes), out, "ra")
