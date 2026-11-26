#!/usr/bin/env Rscript
# Objective: compare BINARY detection of ANY retained viral cluster (0/1)
# across countries and Chinese provinces using Kruskal-Wallis, pooled-rank
# Dunn pairwise tests, and Fisher exact pairwise sensitivity tests.
# IMPORTANT: this does NOT test each individual viral cluster separately.
# Inputs: full run matrices, feature headers of filtered matrices, province map.
# Outputs: overall KW, pairwise Dunn/Fisher, significance and audit TSVs.
# Dependencies: base R. Run from repository root:
#   Rscript scripts/statistical_analysis/kruskal_dunn.R
# BH correction: within condition/scope for pairs, across conditions within
# scope for omnibus KW. Zero-detection runs remain in denominators.
# Responsible contributor: [confirm before submission]
# BIOF3001 Project 7: run-based BINARY viral detection, not richness.
# Comparisons: country vs country; Chinese province vs Chinese province.
# Each SRR is 1 if it detects >=1 retained viral cluster, otherwise 0.
# ALL original SRRs remain (including all-zero runs).
# Conditions: 90/70% x no_singletons/min5/min6.
# Kruskal-Wallis + tie-corrected Dunn; Fisher exact sensitivity check.
# No external R packages required.
# Run from project root: Rscript scripts/statistical_analysis/kruskal_dunn.R

# STEP 1 — Configure paths and group-size threshold. All original strict-
# positive SRRs are kept, even if none detects a retained cluster.
ROOT <- "results"
OUT <- file.path(ROOT, "statistical_analysis", "run_based_analysis", "kw_dunn_binary_detection", "all_runs")
MIN_N <- as.integer(Sys.getenv("MIN_GROUP_N", "5"))
if (is.na(MIN_N) || MIN_N < 2) stop("MIN_GROUP_N must be >= 2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, nrows = -1L) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             nrows = nrows, sep = "\t")
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE,
              col.names = TRUE, quote = FALSE, na = "")
}
clean_group <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | toupper(x) %in% c("", "NA", "N/A", "NONE", "UNKNOWN", "NOT AVAILABLE")] <- NA_character_
  x
}
# STEP 2 — Locate the complete run matrix and the filtered matrix header.
# The latter defines which clusters count as retained for this condition.
condition_paths <- function(id, filter) {
  full <- file.path(ROOT, "viral_clustering", "matrices",
                    sprintf("%s_aS60_aL60_run_with_singletons.tsv", id))
  filtered <- if (filter == "no_singletons") {
    file.path(ROOT, "statistical_analysis", "run_based_analysis", sprintf("%s_aS60_aL60", id),
              "informative_runs_only.tsv")
  } else {
    file.path(ROOT, "statistical_analysis", "run_based_analysis", paste0(filter, "_observations"),
              sprintf("%s_aS60_aL60_informative_runs_%s.tsv", id, filter))
  }
  list(full = full, filtered = filtered)
}
# STEP 3 — Convert selected cluster columns into ONE binary outcome per SRR:
# 1 = detects at least one retained cluster; 0 = detects none.
# This is NOT a per-cluster significance test.
read_detection <- function(id, filter) {
  paths <- condition_paths(id, filter)
  full <- read_tsv(paths$full)
  header <- names(read_tsv(paths$filtered, nrows = 0L))
  features <- grep("^Cluster_", header, value = TRUE)
  if (!length(features)) stop("No retained cluster features in ", paths$filtered)
  if (!all(c("sample", "country") %in% names(full))) stop("Missing sample/country in original matrix")
  if (anyNA(full$sample) || anyDuplicated(full$sample)) stop("Missing/duplicate SRR IDs")
  if (!all(features %in% names(full))) stop("Filtered cluster IDs absent from original matrix: ", paths$filtered)
  mat <- as.matrix(full[, features, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(!is.finite(mat)) || any(!(mat %in% c(0, 1)))) stop("Expected 0/1 cluster matrix")
  data.frame(sample = as.character(full$sample), country = clean_group(full$country),
             detected = as.integer(rowSums(mat) > 0), stringsAsFactors = FALSE)
}
# Dunn ranks are pooled across ALL groups for the given omnibus comparison.
# Tie correction is important because the response consists only of 0s and 1s.
# STEP 4 — Rank pooled binary observations for tie-corrected Dunn
# comparisons. Also compute Fisher exact tests for the 2×2 detection tables;
# Fisher is informative because the outcome is binary.
make_pairs <- function(y, g) {
  g <- factor(g)
  N <- length(y)
  ranks <- rank(y, ties.method = "average")
  ties <- as.numeric(table(y))
  v <- N * (N + 1) / 12 - sum(ties^3 - ties) / (12 * (N - 1))
  combos <- combn(levels(g), 2, simplify = FALSE)
  rows <- lapply(combos, function(groups) {
    a <- y[g == groups[1]]
    b <- y[g == groups[2]]
    ra <- ranks[g == groups[1]]
    rb <- ranks[g == groups[2]]
    z <- if (v > 0) (mean(ra) - mean(rb)) /
      sqrt(v * (1 / length(a) + 1 / length(b))) else 0
    # Fisher is especially appropriate for individual comparisons of binary detection.
    counts <- matrix(c(sum(a), length(a) - sum(a), sum(b), length(b) - sum(b)),
                     nrow = 2, byrow = TRUE)
    fp <- fisher.test(counts)$p.value
    data.frame(group1 = groups[1], group2 = groups[2],
               n1 = length(a), detected1 = sum(a), not_detected1 = length(a) - sum(a),
               prevalence_pct1 = 100 * mean(a),
               n2 = length(b), detected2 = sum(b), not_detected2 = length(b) - sum(b),
               prevalence_pct2 = 100 * mean(b),
               prevalence_difference_pp = 100 * (mean(a) - mean(b)),
               Dunn_Z = z, Dunn_p = 2 * pnorm(-abs(z)), Fisher_p = fp,
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  d$Dunn_p_BH <- p.adjust(d$Dunn_p, method = "BH")
  d$Fisher_p_BH <- p.adjust(d$Fisher_p, method = "BH")
  d
}

# STEP 5 — Join SRRs to provinces, checking that each SRR has one
# unambiguous province assignment.
province_file <- file.path(ROOT, "statistical_analysis", "run_based_analysis", "china_province_audit", "china_srr_provinces.tsv")
prov <- read_tsv(province_file)
if (!all(c("sample", "province") %in% names(prov))) stop("Province map requires sample, province")
if (anyNA(prov$sample) || anyDuplicated(prov$sample)) stop("Duplicate/empty samples in province map")
province_by_sample <- setNames(clean_group(prov$province), as.character(prov$sample))

kw_list <- list(); pair_list <- list(); audit_list <- list()
# STEP 6 — Iterate over 90%/70% clustering and three frequency
# thresholds; test countries and China provinces separately.
for (id in c(90, 70)) for (filter in c("no_singletons", "min5", "min6")) {
  runs <- read_detection(id, filter)
  runs$province <- unname(province_by_sample[runs$sample])
  prefix <- sprintf("%s_%s", id, filter)
  write_tsv(runs, file.path(OUT, paste0(prefix, "_all_SRR_detection.tsv")))
  audit_list[[length(audit_list) + 1L]] <- data.frame(
    identity = id, filter = filter, total_original_SRRs = nrow(runs),
    runs_with_any_retained_cluster = sum(runs$detected),
    zero_detection_runs = sum(runs$detected == 0),
    stringsAsFactors = FALSE)
  for (scope in c("country", "china_province")) {
    if (scope == "country") {
      work <- runs
      work$group <- work$country
    } else {
      work <- runs[!is.na(runs$country) & tolower(runs$country) == "china", , drop = FALSE]
      work$group <- clean_group(work$province)
    }
    work <- work[!is.na(work$group), , drop = FALSE]
    group_counts <- table(work$group)
    work <- work[work$group %in% names(group_counts[group_counts >= MIN_N]), , drop = FALSE]
    work$group <- factor(work$group)
    if (nlevels(work$group) < 2L) {
      message("Skipping ", prefix, " ", scope, ": fewer than 2 eligible groups")
      next
    }
    y <- work$detected
    # STEP 7 — Perform Kruskal–Wallis omnibus comparison of group
    # detection rates, with a neutral result if every SRR has the same outcome.
    kw <- if (length(unique(y)) > 1L) kruskal.test(y, work$group) else NULL
    kw_h <- if (is.null(kw)) 0 else as.numeric(kw$statistic)
    kw_p <- if (is.null(kw)) 1 else kw$p.value
    kw_list[[length(kw_list) + 1L]] <- data.frame(
      identity = id, filter = filter, scope = scope,
      response = "any_retained_viral_cluster_detected_0_1",
      included_SRRs = nrow(work), n_groups = nlevels(work$group),
      KW_H = kw_h, KW_df = nlevels(work$group) - 1L, KW_p = kw_p,
      stringsAsFactors = FALSE)
    p <- make_pairs(y, work$group)
    p$identity <- id; p$filter <- filter; p$scope <- scope
    pair_list[[length(pair_list) + 1L]] <- p
    message(prefix, " ", scope, ": N=", nrow(work), ", positive=", sum(y),
            ", KW p=", signif(kw_p, 4), ", Dunn BH hits=", sum(p$Dunn_p_BH < 0.05))
  }
}
if (!length(kw_list) || !length(pair_list)) stop("No valid analyses")
KW <- do.call(rbind, kw_list)
DN <- do.call(rbind, pair_list)
# Correct the six omnibus tests separately for country and province scopes.
# STEP 8 — Correct for multiple comparisons. BH is applied
# separately to the country and province omnibus test families.
KW$KW_p_BH_scope <- ave(KW$KW_p, KW$scope,
                         FUN = function(x) p.adjust(x, method = "BH"))
combine_key <- function(x) paste(x$identity, x$filter, x$scope, sep = "|")
idx <- match(combine_key(DN), combine_key(KW))
if (anyNA(idx)) stop("Internal omnibus/Dunn join failed")
DN$KW_p <- KW$KW_p[idx]
DN$KW_p_BH_scope <- KW$KW_p_BH_scope[idx]
DN$significant_KW_and_Dunn_BH <- DN$KW_p_BH_scope < 0.05 & DN$Dunn_p_BH < 0.05
DN$significant_Fisher_BH <- DN$Fisher_p_BH < 0.05
# STEP 9 — Export all omnibus/pairwise results, selected significant
# comparisons, a denominator audit, and a precise methods note.
write_tsv(KW, file.path(OUT, "ALL_kruskal.tsv"))
write_tsv(DN, file.path(OUT, "ALL_dunn_fisher_pairs.tsv"))
write_tsv(DN[DN$significant_KW_and_Dunn_BH, , drop = FALSE],
          file.path(OUT, "SIGNIFICANT_KW_DUNN_PAIRS.tsv"))
write_tsv(DN[DN$significant_Fisher_BH, , drop = FALSE],
          file.path(OUT, "SIGNIFICANT_FISHER_PAIRS.tsv"))
write_tsv(do.call(rbind, audit_list), file.path(OUT, "INPUT_AUDIT.tsv"))
writeLines(c(
  "UNIT: individual SRR. RESPONSE: 1 if at least one retained viral cluster detected, else 0.",
  "ALL ORIGINAL SRRs INCLUDED (including zero-detection runs); only geography groups with >= MIN_GROUP_N runs are tested.",
  "Retained features are the exact columns from existing no_singletons/min5/min6 filtered matrix headers.",
  "Each condition/scope uses a tie-corrected omnibus KW and pooled-rank Dunn pairs.",
  "Dunn and Fisher BH adjustments are within each identity x filtering condition x geography scope.",
  "Omnibus KW BH adjustment is across identity/filter conditions separately within country and province scopes.",
  "Fisher exact tests are recommended as sensitivity/primary pairwise analysis for a binary outcome.",
  "Zero means no retained viral cluster detected in the data, not proof of virus absence.",
  "The tests compare detection prevalence only; viral identity and community composition are not tested."
), file.path(OUT, "METHODS.txt"))
cat("DONE. Results at: ", OUT, "\n", sep = "")
