# BIOF3001 Project 7 — Shared helpers for run-based statistical analyses.
# Objective: provide standard matrix paths, province mapping, binary cluster
# detection loading, group-size screening and TSV saving for downstream R tests
# (PERMANOVA, PERMDISP, pairwise and other analyses).
# Inputs: results/statistical_analysis/run_based_analysis/ matrices and mapping.
# Outputs: none by itself; helper functions are sourced by downstream R scripts.
# Requirements: base R. Source from repository root with:
#   source("scripts/statistical_analysis/common.R")
# Settings: N_PERM=9999, MIN_GROUP_N=5, SEED=42 (overridable via environment).
# Responsible contributor: [confirm before submission]

# Shared helpers for the 2026 run-based viral-community analyses.
# Run all R scripts from the repository root.
BASE <- "results/statistical_analysis/run_based_analysis"
PERM <- as.integer(Sys.getenv("N_PERM", "9999"))
MIN_N <- as.integer(Sys.getenv("MIN_GROUP_N", "5"))
SEED <- as.integer(Sys.getenv("SEED", "42"))
CONDITIONS <- expand.grid(identity = c(90L, 70L),
                          filter = c("no_singletons", "min5", "min6"),
                          stringsAsFactors = FALSE)

matrix_path <- function(identity, filter) {
  condition <- paste0(identity, "_aS60_aL60")
  if (filter == "no_singletons")
    return(file.path(BASE, condition, "informative_runs_only.tsv"))
  file.path(BASE, paste0(filter, "_observations"),
            paste0(condition, "_informative_runs_", filter, ".tsv"))
}

province_lookup <- function() {
  file <- file.path(BASE, "china_province_audit", "china_srr_provinces.tsv")
  if (!file.exists(file)) stop("Missing province mapping: ", file)
  p <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("sample", "province") %in% names(p))) stop("Province mapping needs sample and province")
  if (anyDuplicated(p$sample)) stop("Duplicate sample IDs in province metadata")
  setNames(trimws(as.character(p$province)), as.character(p$sample))
}

read_run <- function(identity, filter) {
  path <- matrix_path(identity, filter)
  if (!file.exists(path)) stop("Missing run-level matrix: ", path)
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("sample", "country") %in% names(d))) stop("Missing sample/country in ", path)
  if (anyNA(d$sample) || anyDuplicated(d$sample)) stop("Missing or duplicate SRRs in ", path)
  cls <- grep("^Cluster_", names(d), value = TRUE)
  if (!length(cls)) stop("No clusters in ", path)
  X <- as.matrix(d[, cls, drop = FALSE]); suppressWarnings(storage.mode(X) <- "numeric")
  if (any(!is.finite(X)) || any(X < 0)) stop("Invalid cluster matrix in ", path)
  X <- (X > 0) * 1L
  good <- rowSums(X) > 0
  d <- d[good, , drop = FALSE]
  d$country <- trimws(as.character(d$country))
  list(data=d, X=X[good, , drop = FALSE], cls=cls)
}

subsets <- function(m, province_map) {
  d <- m$data
  china <- !is.na(d$country) & tolower(d$country) == "china"
  list(
    country=list(X=m$X, sample=d$sample, group=d$country),
    china_province=list(X=m$X[china, , drop=FALSE], sample=d$sample[china],
                        group=unname(province_map[as.character(d$sample[china])]))
  )
}

eligible_subset <- function(obj) {
  g <- trimws(as.character(obj$group))
  ok <- !is.na(g) & nzchar(g) & toupper(g) != "UNKNOWN"
  g <- g[ok]; X <- obj$X[ok, , drop=FALSE]; sam <- obj$sample[ok]
  n_pre <- length(g)
  tabs <- table(g); valid <- names(tabs[tabs >= MIN_N]);keep <- g %in% valid
  g <- factor(g[keep]); X <- X[keep, , drop=FALSE]; sam <- sam[keep]
  if (nlevels(g)<2L || nrow(X)<=nlevels(g)) return(NULL)
  # All-zero rows are excluded prior to distance calculation. No all-zero rows should
  # appear because read_run() already filtered them.
  if (any(rowSums(X)==0)) stop("Unexpected all-zero SRR in PERMANOVA analysis")
  list(X=X,g=g,sample=sam,SRRs_before_filter=n_pre)
}

save_tsv <- function(d, path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  write.table(d,path,sep="\t",quote=FALSE,row.names=FALSE,na="")
}
