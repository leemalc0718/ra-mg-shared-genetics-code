#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
root <- dirname(dirname(normalizePath(script_file)))
source(file.path(root, "R", "common.R")); source(file.path(root, "R", "scrna_common.R"))
required_packages(c("Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "patchwork"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 16_scrna_mg.R config.yml", call. = FALSE)
cfg <- load_project_config(args[[1]])
sc <- cfg$single_cell$mg
files <- list.files(project_path(cfg, sc$input_dir, TRUE), pattern = sc$file_pattern, full.names = TRUE)
if (!length(files)) stop("No MG single-cell matrices matched the configured pattern.", call. = FALSE)

objects <- lapply(files, function(path) {
  sample_name <- sub("\\.txt(\\.gz)?$", "", basename(path))
  raw <- data.table::fread(path, data.table = FALSE)
  genes <- make.unique(as.character(raw[[1]]))
  matrix_data <- as.matrix(raw[, -1, drop = FALSE]); storage.mode(matrix_data) <- "double"
  rownames(matrix_data) <- genes
  counts <- methods::as(matrix_data, "dgCMatrix")
  object <- Seurat::CreateSeuratObject(counts, project = sample_name, min.cells = 3, min.features = 200)
  object$sample <- sample_name
  object$group <- if (grepl(sc$disease_pattern, sample_name)) "MG" else if (grepl(sc$control_pattern, sample_name)) "Control" else "Unknown"
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
  subset(object, subset = nFeature_RNA > sc$min_features & nFeature_RNA < sc$max_features & percent.mt < sc$max_percent_mt)
})
names(objects) <- vapply(objects, function(x) as.character(unique(x$sample)), character(1))
object <- if (length(objects) == 1L) objects[[1]] else merge(objects[[1]], y = objects[-1], add.cell.ids = names(objects))
out <- output_path(cfg, "scrna", "mg"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
object <- run_scrna_workflow(object, sc$dims, cfg$parameters$scrna$variable_features,
  cfg$parameters$scrna$clustering_resolution, cfg$project$seed, out)
object <- annotate_clusters(object, project_path(cfg, sc$annotation_file, TRUE))
save_scrna_outputs(object, unlist(cfg$parameters$scrna$candidate_genes), out, "mg")
