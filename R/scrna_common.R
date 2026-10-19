run_scrna_workflow <- function(object, dims, variable_features, resolution, seed, output_dir) {
  object[["RNA"]] <- SeuratObject::JoinLayers(object[["RNA"]])
  Seurat::DefaultAssay(object) <- "RNA"
  set.seed(seed)
  object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = variable_features)
  object <- Seurat::ScaleData(object, features = rownames(object))
  object <- Seurat::RunPCA(object, features = Seurat::VariableFeatures(object), npcs = max(50, dims))
  elbow <- Seurat::ElbowPlot(object, ndims = max(50, dims))
  ggplot2::ggsave(file.path(output_dir, "elbow_plot.pdf"), elbow, width = 6, height = 5)
  object <- Seurat::FindNeighbors(object, dims = seq_len(dims))
  object <- Seurat::FindClusters(object, resolution = resolution)
  object <- Seurat::RunUMAP(object, dims = seq_len(dims))
  object <- Seurat::RunTSNE(object, dims = seq_len(dims))
  object
}

annotate_clusters <- function(object, mapping_path) {
  mapping <- data.table::fread(mapping_path, colClasses = c(cluster = "character"))
  required <- c("cluster", "detailed_cell_type", "broad_cell_type")
  if (!all(required %in% names(mapping))) stop("Annotation table requires: ", paste(required, collapse = ", "), call. = FALSE)
  observed <- sort(unique(as.character(Seurat::Idents(object))))
  supplied <- sort(unique(mapping$cluster))
  if (!identical(observed, supplied)) {
    stop("Observed clusters do not match the annotation table. Observed: ", paste(observed, collapse = ", "),
         "; supplied: ", paste(supplied, collapse = ", "), call. = FALSE)
  }
  idx <- match(as.character(Seurat::Idents(object)), mapping$cluster)
  object$detailed_cell_type <- mapping$detailed_cell_type[idx]
  object$broad_cell_type <- mapping$broad_cell_type[idx]
  Seurat::Idents(object) <- "broad_cell_type"
  object
}

save_scrna_outputs <- function(object, genes, output_dir, prefix) {
  cluster_plot <- Seurat::DimPlot(object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggplot2::theme_classic()
  celltype_plot <- Seurat::DimPlot(object, reduction = "umap", group.by = "broad_cell_type", label = TRUE, repel = TRUE) + ggplot2::theme_classic()
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_clusters_umap.pdf")), cluster_plot, width = 8, height = 7)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_celltypes_umap.pdf")), celltype_plot, width = 9, height = 7)

  object[["RNA"]] <- SeuratObject::JoinLayers(object[["RNA"]])
  Seurat::Idents(object) <- "seurat_clusters"
  markers <- Seurat::FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  data.table::fwrite(markers, file.path(output_dir, paste0(prefix, "_cluster_markers.tsv.gz")), sep = "\t")
  Seurat::Idents(object) <- "broad_cell_type"
  present <- intersect(genes, rownames(object))
  missing <- setdiff(genes, present)
  if (length(missing)) warning("Candidate genes absent from the expression matrix: ", paste(missing, collapse = ", "))
  if (length(present)) {
    plots <- lapply(present, function(gene) {
      Seurat::VlnPlot(object, features = gene, group.by = "broad_cell_type", pt.size = 0) +
        ggplot2::ggtitle(gene) + ggplot2::labs(x = "Cell type", y = "Normalized expression") +
        ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
    })
    combined <- patchwork::wrap_plots(plots, ncol = 2)
    ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_candidate_gene_expression.pdf")), combined, width = 10, height = 10)
  }
  saveRDS(object, file.path(output_dir, paste0(prefix, "_annotated.rds")))
  object
}
