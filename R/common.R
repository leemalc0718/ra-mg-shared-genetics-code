required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

load_project_config <- function(path) {
  required_packages("yaml")
  path <- normalizePath(path, mustWork = TRUE)
  cfg <- yaml::read_yaml(path)
  attr(cfg, "config_path") <- path
  attr(cfg, "repo_root") <- dirname(dirname(path))
  cfg
}

project_path <- function(cfg, path, must_work = FALSE) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("A non-empty path is required.", call. = FALSE)
  }
  expanded <- path.expand(path)
  if (!grepl("^/", expanded)) {
    expanded <- file.path(attr(cfg, "repo_root"), expanded)
  }
  normalizePath(expanded, mustWork = must_work)
}

tool_path <- function(cfg, value) {
  if (grepl("/", value, fixed = TRUE)) return(project_path(cfg, value, must_work = TRUE))
  resolved <- Sys.which(value)
  if (!nzchar(resolved)) stop("Tool not found on PATH: ", value, call. = FALSE)
  unname(resolved)
}

output_path <- function(cfg, ...) {
  root <- project_path(cfg, cfg$project$output_dir, must_work = FALSE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  file.path(root, ...)
}

read_standard_gwas <- function(path) {
  required_packages("data.table")
  dat <- data.table::fread(path)
  required <- c("SNP", "CHR", "BP", "A1", "A2", "beta", "se", "p", "eaf", "n", "Z")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("Missing standardized GWAS columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  numeric_cols <- c("CHR", "BP", "beta", "se", "p", "eaf", "n", "Z")
  dat[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]
  dat[, `:=`(SNP = as.character(SNP), A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)))]
  dat
}

validate_standard_gwas <- function(dat, label = "GWAS") {
  problems <- character()
  if (anyDuplicated(dat$SNP)) problems <- c(problems, "duplicated SNP identifiers")
  if (any(!dat$CHR %in% 1:22, na.rm = TRUE)) problems <- c(problems, "non-autosomal CHR values")
  if (any(dat$BP < 1, na.rm = TRUE)) problems <- c(problems, "invalid BP values")
  valid_alleles <- c("A", "C", "G", "T")
  if (any(!dat$A1 %in% valid_alleles | !dat$A2 %in% valid_alleles)) problems <- c(problems, "non-SNV alleles")
  if (any(dat$A1 == dat$A2)) problems <- c(problems, "identical A1 and A2")
  if (any(!is.finite(dat$beta)) || any(!is.finite(dat$se)) || any(dat$se <= 0)) problems <- c(problems, "invalid beta or se")
  if (any(!is.finite(dat$p)) || any(dat$p <= 0 | dat$p > 1)) problems <- c(problems, "invalid p values")
  if (any(!is.finite(dat$eaf)) || any(dat$eaf <= 0 | dat$eaf >= 1)) problems <- c(problems, "invalid eaf values")
  if (any(!is.finite(dat$n)) || any(dat$n <= 0)) problems <- c(problems, "invalid sample sizes")
  if (length(problems) > 0L) {
    stop(label, " failed validation: ", paste(unique(problems), collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

case_fraction <- function(trait_cfg) {
  trait_cfg$n_case / (trait_cfg$n_case + trait_cfg$n_control)
}

harmonize_pair <- function(x, y, suffixes = c(".x", ".y"), drop_palindromic = TRUE) {
  required_packages("data.table")
  dat <- merge(x, y, by = "SNP", suffixes = suffixes)
  x1 <- paste0("A1", suffixes[1]); x2 <- paste0("A2", suffixes[1])
  y1 <- paste0("A1", suffixes[2]); y2 <- paste0("A2", suffixes[2])
  complement <- c(A = "T", T = "A", C = "G", G = "C")
  y1_comp <- unname(complement[dat[[y1]]])
  y2_comp <- unname(complement[dat[[y2]]])
  dat[, allele_match := (get(x1) == get(y1) & get(x2) == get(y2)) |
                          (get(x1) == y1_comp & get(x2) == y2_comp)]
  dat[, allele_flip := (get(x1) == get(y2) & get(x2) == get(y1)) |
                         (get(x1) == y2_comp & get(x2) == y1_comp)]
  if (drop_palindromic) {
    dat[, palindromic := paste0(get(x1), get(x2)) %in% c("AT", "TA", "CG", "GC")]
  } else {
    dat[, palindromic := FALSE]
  }
  dat <- dat[(allele_match | allele_flip) & !palindromic]
  dat
}

script_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg) != 1L) stop("Unable to locate the running script.", call. = FALSE)
  dirname(dirname(normalizePath(file_arg)))
}
