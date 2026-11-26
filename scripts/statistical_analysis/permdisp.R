#!/usr/bin/env Rscript
# Objective: test within-group multivariate dispersion overall and by
# country/province pairs using Bray-Curtis distance to median.
# Input: filtered run matrices and China province map; common.R helpers.
# Output: overall and pairwise PERMDISP tables.
# Dependencies: R, vegan. Usage from repository root:
#   Rscript scripts/statistical_analysis/permdisp.R
# Uses Lingoes correction in betadisper; unlike PERMANOVA above,
# this is a dispersion test, not a test for differences in centres.
# Responsible contributor: [confirm before submission]
# PERMDISP on ORIGINAL Bray-Curtis dissimilarities with Lingoes correction
# applied internally for distance-to-median calculation, via betadisper(add).
if (!requireNamespace('vegan',quietly=TRUE)) stop('Install R package vegan')
source('scripts/statistical_analysis/common.R')
set.seed(SEED)
OUT <- file.path(BASE,'pairwise_permanova_braycurtis')
GLOBAL_OUT <- file.path(BASE,'permanova_braycurtis_lingoes')
all_overall <- list(); all_pairs <- list()

disperse <- function(X,g) {
  D <- vegan::vegdist(X,method='bray')
  fit <- vegan::betadisper(D,g,type='median',bias.adjust=TRUE,add='lingoes')
  perm <- vegan::permutest(fit,permutations=PERM)
  tab <- as.data.frame(perm$tab)
  c(F=as.numeric(tab[1,'F']),p=as.numeric(tab[1,'Pr(>F)']))
}

for(i in seq_len(nrow(CONDITIONS))) {
  id <- CONDITIONS$identity[i];fil <- CONDITIONS$filter[i]
  m <- read_run(id,fil);p <- province_lookup()
  for(analysis in names(subsets(m,p))) {
    d <- eligible_subset(subsets(m,p)[[analysis]])
    if (is.null(d))next
    overall <- disperse(d$X,d$g)
    all_overall[[length(all_overall)+1L]] <- data.frame(
      identity=id,filter=fil,analysis=analysis,n_SRRs=nrow(d$X),
      n_groups=nlevels(d$g),PERMDISP_F=overall['F'],PERMDISP_p=overall['p'])
    res <- list()
    for(pair in combn(levels(d$g),2,simplify=FALSE)) {
      keep <- d$g %in% pair;group <- droplevels(d$g[keep])
      X <- d$X[keep,,drop=FALSE];r <- disperse(X,group)
      res[[length(res)+1L]] <- data.frame(
        identity=id,filter=fil,analysis=analysis,
        group1=pair[1],group2=pair[2],
        n_group1=sum(group==pair[1]),n_group2=sum(group==pair[2]),
        PERMDISP_F=r['F'],PERMDISP_p=r['p'])
    }
    tmp <- do.call(rbind,res)
    tmp$PERMDISP_p_BH <- p.adjust(tmp$PERMDISP_p,method='BH')
    all_pairs[[length(all_pairs)+1L]] <- tmp
    message(id,' ',fil,' ',analysis,': ',nrow(tmp),' pairwise PERMDISP comparisons')
  }
}
if(!length(all_overall))stop('No PERMDISP analyses succeeded')
overall <- do.call(rbind,all_overall)
overall$PERMDISP_p_BH <- p.adjust(overall$PERMDISP_p,method='BH')
save_tsv(overall,file.path(GLOBAL_OUT,'overall_permdisp_lingoes.tsv'))
if(length(all_pairs)) {
  d <- do.call(rbind,all_pairs)
  d$PERMDISP_p_BH_global <- p.adjust(d$PERMDISP_p,method='BH')
  save_tsv(d,file.path(OUT,'pairwise_permdisp_lingoes.tsv'))
}
message('PERMDISP complete: ',OUT)
