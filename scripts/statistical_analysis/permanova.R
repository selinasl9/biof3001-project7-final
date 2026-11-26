#!/usr/bin/env Rscript
# Objective: overall and pairwise PERMANOVA of Bray-Curtis binary viral-cluster
# profiles among countries and Chinese provinces.
# Input: filtered run matrices and China province map; common.R functions.
# Outputs: overall PERMANOVA, group counts, pairwise tables.
# Dependencies: R, vegan. Usage from repository root:
#   Rscript scripts/statistical_analysis/permanova.R
# Parameters: N_PERM, MIN_GROUP_N, SEED (environment variables).
# Statistical caution: geography and BioProject can be confounded;
# unrestricted permutations do not control for such confounding.
# Responsible contributor: [confirm before submission]
# Original binary Bray-Curtis PERMANOVA: overall and pairwise.
# p_BH is within (identity,filter,analysis). p_BH_global is across all pairwise tests.
# IMPORTANT: Batch/BioProject and geography may be confounded; permutations are
# unrestricted here only because project stratification cannot separate them.
if (!requireNamespace('vegan',quietly=TRUE)) stop('Install R package vegan')
source('scripts/statistical_analysis/common.R')
set.seed(SEED)
OUT <- file.path(BASE,'pairwise_permanova_braycurtis')
GLOBAL_OUT <- file.path(BASE,'permanova_braycurtis_lingoes')
all_pairs <- list(); all_overall <- list();all_counts <- list()

adonis_pair <- function(X,g) {
  # Preserve original sample-level Bray-Curtis, not NMDS/t-SNE coordinates.
  D <- vegan::vegdist(X,method='bray')
  dd <- data.frame(group=g)
  a <- as.data.frame(vegan::adonis2(D~group,data=dd,permutations=PERM))
  c(F=as.numeric(a[1,'F']), R2=as.numeric(a[1,'R2']), p=as.numeric(a[1,'Pr(>F)']))
}

for (i in seq_len(nrow(CONDITIONS))) {
  id <- CONDITIONS$identity[i]; fil <- CONDITIONS$filter[i]
  m <- read_run(id,fil); p <- province_lookup()
  for (analysis in names(subsets(m,p))) {
    d <- eligible_subset(subsets(m,p)[[analysis]])
    if (is.null(d)) { message('SKIP: ',id,' ',fil,' ',analysis);next }
    overall <- adonis_pair(d$X,d$g)
    all_overall[[length(all_overall)+1L]] <- data.frame(
      identity=id,filter=fil,analysis=analysis,n_SRRs=nrow(d$X),
      n_groups=nlevels(d$g),n_clusters=ncol(d$X),permutations=PERM,
      pseudo_F=overall['F'],R2=overall['R2'],p=overall['p'])
    for (group in levels(d$g))
      all_counts[[length(all_counts)+1L]] <- data.frame(
        identity=id,filter=fil,analysis=analysis,
        group=group,n_SRRs=sum(d$g==group))
    pairs <- combn(levels(d$g),2,simplify=FALSE)
    condition_rows <- list()
    for (pair in pairs) {
      keep <- d$g %in% pair
      group <- droplevels(d$g[keep]); X <- d$X[keep,,drop=FALSE]
      r <- adonis_pair(X,group)
      condition_rows[[length(condition_rows)+1L]] <- data.frame(
        identity=id,filter=fil,analysis=analysis,
        group1=pair[1],group2=pair[2],
        n_group1=sum(group==pair[1]),n_group2=sum(group==pair[2]),
        n_clusters=ncol(X),pseudo_F=r['F'],R2=r['R2'],p=r['p'])
    }
    a <- do.call(rbind,condition_rows)
    a$p_BH <- p.adjust(a$p,method='BH')
    all_pairs[[length(all_pairs)+1L]] <- a
    message(id,' ',fil,' ',analysis,': ',nrow(a),' pairwise comparisons')
  }
}
if (!length(all_overall)) stop('No analyses succeeded')
overall <- do.call(rbind,all_overall)
overall$p_BH <- p.adjust(overall$p,method='BH')
save_tsv(overall,file.path(GLOBAL_OUT,'overall_permanova.tsv'))
save_tsv(do.call(rbind,all_counts),file.path(GLOBAL_OUT,'group_counts.tsv'))
if (length(all_pairs)) {
  pairs <- do.call(rbind,all_pairs)
  pairs$p_BH_global <- p.adjust(pairs$p,method='BH')
  save_tsv(pairs,file.path(OUT,'pairwise_permanova_results.tsv'))
}
message('PERMANOVA complete: ',OUT)
