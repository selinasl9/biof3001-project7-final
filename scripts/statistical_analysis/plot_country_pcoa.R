#!/usr/bin/env Rscript
# Objective: descriptive country-level PCoA of viral-cluster prevalence
# profiles (NOT the run-level PERMANOVA input).
# Input: country prevalence matrices for 90%/70% aS60_aL60 clustering.
# Output: country PCoA coordinate TSVs and PNG figures.
# Dependencies: R, vegan, ggplot2. Run from repository root:
#   Rscript scripts/statistical_analysis/plot_country_pcoa.R
# Lingoes correction is applied to the Bray-Curtis dissimilarity.
# Responsible contributor: [confirm before submission]
# Country-level PCoA using prevalence fractions for country-singleton-filtered
# viral clusters. Descriptive country-level visualization, NOT SRR-level inference.
# STEP 1 — Validate R plotting/ordination dependencies, then create
# the output folder for descriptive country-level figures.
if (!requireNamespace('vegan',quietly=TRUE))stop('Install vegan')
if (!requireNamespace('ggplot2',quietly=TRUE))stop('Install ggplot2')
out <- 'results/statistical_analysis/country_pcoa_final'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
# STEP 2 — Load the 90% and 70% COUNTRY prevalence matrices in turn.
# Each row is a country, rather than an individual sequencing run.
for(id in c(90,70)) {
  cond <- paste0(id,'_aS60_aL60')
  path <- file.path('results/viral_clustering/matrices',paste0(cond,'_country_with_singletons.tsv'))
  if(!file.exists(path))stop('Missing ',path)
  df <- read.delim(path,check.names=FALSE)
  cls <- grep('^Cluster_',names(df),value=TRUE)
  # A shared cluster must appear in at least two countries BEFORE min-n filtering.
  # STEP 3 — Retain clusters present in at least two countries BEFORE
  # excluding countries with fewer than five strict-positive SRRs.
  shared <- cls[colSums(df[cls]>0)>=2]
  d <- df[df$n_samples>=5,,drop=FALSE]
  X <- as.matrix(d[,shared,drop=FALSE]);storage.mode(X) <- 'double'
  keep <- rowSums(X)>0; X <- X[keep,,drop=FALSE];d <- d[keep,,drop=FALSE]
  if(nrow(X)<3L || ncol(X)<2L){message('Skipping ',cond,': insufficient information');next}
  # STEP 4 — Compute Bray–Curtis dissimilarity between country
  # prevalence profiles and apply Lingoes correction for PCoA.
  bc <- vegan::vegdist(X,method='bray')
  pc <- vegan::wcmdscale(bc,k=2,eig=TRUE,add='lingoes')
  # STEP 5 — Save coordinates for audit/replotting and annotate the
  # country plot with explained positive-eigenvalue proportions.
  coords <- data.frame(country=d$country,n_SRRs=d$n_samples,
                       PCoA1=pc$points[,1],PCoA2=pc$points[,2])
  write.table(coords,file.path(out,paste0(cond,'_pcoa_coordinates.tsv')),
              sep='\t',row.names=FALSE,quote=FALSE)
  pct <- round(100*pc$eig[1:2]/sum(pc$eig[pc$eig>0]),1)
  # STEP 6 — Render country labels and the first two PCoA axes.
  plot <- ggplot2::ggplot(coords,ggplot2::aes(PCoA1,PCoA2,label=country)) +
    ggplot2::geom_point(size=3,color='#299fa5')+
    ggplot2::geom_text(vjust=-.8,size=3.3,check_overlap=TRUE)+
    ggplot2::labs(title=paste(id,'% country prevalence PCoA'),
                  x=paste0('PCoA1 (',pct[1],'%)'),y=paste0('PCoA2 (',pct[2],'%)'))+
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out,paste0(cond,'_pcoa.png')),plot,
                  width=10,height=7,dpi=250,bg='white')
}
