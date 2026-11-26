#!/usr/bin/env Rscript
#
# render_virus_tree_ggtree.R
#
# Draw the virus tree with ggtree, midpoint rooted, with colour strips beside the tips.
#
# Inputs:
#   TREEFILE       the virus tree, defaults to the tree_gt08 tree
#   ANNOTATION     contig_annotation.tsv, country, province, host_genus, host_group, cluster
#   REFS           refs_moderate_only_msa.fasta + blast_dedup.fa, for the reference names
#
# Outputs:
#   PREFIX_ggtree.pdf and _ggtree.png
#
# Flags:
#   --china        China contigs only
#   --all-labels   name every downloaded reference as well

suppressPackageStartupMessages({
  library(grid)
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(ggtree)
  library(ggnewscale)
})

# Step 1: read the arguments + the two flags
args <- commandArgs(trailingOnly = TRUE)
china_only <- "--china" %in% args
all_labels <- "--all-labels" %in% args     # name every downloaded reference as well
args <- args[!(args %in% c("--china", "--all-labels"))]

tree_file <- if (length(args) > 0) args[1] else "results/msa/complete/tree_gt08/tree_gt08.treefile"
annot_file <- if (length(args) > 1) args[2] else "results/msa/complete/tree_gt08/contig_annotation.tsv"
out_prefix <- if (length(args) > 2) args[3] else sub("\\.treefile$", "", tree_file)
# Set the reference FASTA files, the tips to drop + the colours for each strip
# The colours mirror common/palette.py, which R cannot read, so a change has to be made twice
refs_fasta <- c("results/msa/complete/refs_moderate_only_msa.fasta",
                "results/msa/tree_v2/blast_dedup.fa")

# These three sit on very long branches and drag the midpoint root onto the diatom virus
DROP <- c("NC_040589.1", "NC_040555.1", "NC_003996.1")

# Countries match Cammy's balanced beta analysis, only those in the tree reach the legend
COUNTRY_COLOURS <- c("China" = "#E6194B", "Australia" = "#3CB44B", "USA" = "#4363D8",
                     "Finland" = "#F58231", "Brazil" = "#911EB4", "Trinidad and Tobago" = "#42D4F4",
                     "South Korea" = "#F032E6", "Egypt" = "#9A6324", "Canada" = "#808000",
                     "Democratic Republic of the Congo" = "#469990", "Viet Nam" = "#000075",
                     "other" = "#555555")
HOST_COLOURS <- c("Culex" = "#1F78B4", "Aedes" = "#E31A1C", "unknown" = "#E8E8E8")
# Provinces are assigned in alphabetical order, so the figure is reproducible
PROVINCE_PALETTE <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948",
                      "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#1B9E77", "#D95F02",
                      "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
# Cluster colours recycle over the 28 clusters, as the cluster number is on the label too
CLUSTER_PALETTE <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                     "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

# Step 2: read the tree
tree <- read.tree(tree_file)
# Strip _R_, which MAFFT adds to anything it reverse complemented, not part of the name
tree$tip.label <- sub("^_R_", "", tree$tip.label)

# Step 3: read the annotation table, split by hand as read.delim trips over a lost tab
lines <- readLines(annot_file)
lines <- lines[nzchar(trimws(lines))]
header <- strsplit(sub("\\s+$", "", lines[1]), "\t", fixed = TRUE)[[1]]
if (length(header) == 1) header <- strsplit(trimws(lines[1]), "[ \t]+")[[1]]
rows <- lapply(strsplit(lines[-1], "\t", fixed = TRUE), function(x) {
  length(x) <- length(header)
  trimws(x)
})
annot <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(annot) <- header
cat(nrow(annot), "rows in", annot_file, "\n")
annot$contig <- sub("^_R_", "", annot$contig)

# Step 4: fill in the columns an older annotation table did not carry
# Work the genus out here, as older tables carry only host_organism
if (!"host_genus" %in% names(annot)) {
  sp <- if ("host_species" %in% names(annot)) annot$host_species else
        if ("host_organism" %in% names(annot)) annot$host_organism else ""
  genus <- sub(" .*", "", sp)
  known <- c("Culex", "Aedes", "Ochlerotatus", "Anopheles", "Armigeres",
             "Mansonia", "Coquillettidia", "Culiseta", "Uranotaenia")
  genus[!genus %in% known] <- "unknown"
  annot$host_genus <- genus
}
if (!"host_group" %in% names(annot)) {
  # Pool Ochlerotatus into Aedes, now treated as a subgenus of it
  annot$host_group <- ifelse(annot$host_genus %in% c("Aedes", "Ochlerotatus"), "Aedes",
                      ifelse(annot$host_genus == "Culex", "Culex", "unknown"))
}
if ("province" %in% names(annot)) {
  # Fix the Beijing samples, which are recorded by district
  annot$province <- sub("^Dongcheng District$", "Beijing", annot$province)
  annot$province[annot$country != "China"] <- ""
}
rownames(annot) <- annot$contig

# Step 5: drop the three long-branch references, then prune to China if asked
drop <- intersect(DROP, tree$tip.label)
if (length(drop)) tree <- drop.tip(tree, drop)
cat("dropped", length(drop), "sequences:", paste(drop, collapse = " "), "\n")

if (china_only) {
  # Keep every reference, but only the contigs collected in China
  contigs <- intersect(tree$tip.label, annot$contig)
  away <- contigs[annot[contigs, "country"] != "China"]
  tree <- drop.tip(tree, away)
  cat("China only:", length(away), "contigs from other countries dropped,",
      length(tree$tip.label), "tips left\n")
  out_prefix <- paste0(out_prefix, "_china")
}

# Step 6: midpoint root the tree, now the long branches are gone
tree <- midpoint(tree)
cat("midpoint rooted,", length(tree$tip.label), "tips\n")

# Step 7: build the reference labels, accession + the virus name from the FASTA header
ref_name <- c()
for (f in refs_fasta[file.exists(refs_fasta)]) {
  heads <- grep("^>", readLines(f), value = TRUE)
  acc <- sub("^>(\\S+).*", "\\1", heads)
  nm <- sub("^>\\S+\\s*", "", heads)
  nm <- sub(",.*$", "", nm)
  nm <- sub("\\s+(strain|isolate|complete|genomic|polyprotein).*$", "", nm)
  ref_name <- c(ref_name, setNames(nm, acc))
}
ref_name <- ref_name[!duplicated(names(ref_name))]

# Step 8: work out which tips are ours + build every tip's label
is_contig <- tree$tip.label %in% annot$contig
if (!any(is_contig)) {
  cat("tree tips:      ", paste(head(tree$tip.label, 3), collapse = ", "), "\n")
  cat("table contigs:  ", paste(head(annot$contig, 3), collapse = ", "), "\n")
  cat("table columns:  ", paste(names(annot), collapse = ", "), "\n")
  stop("no tip name matches the annotation table, check ", annot_file)
}
cat(sum(is_contig), "of", length(tree$tip.label), "tips are our contigs\n")
place <- annot[tree$tip.label, "country"]
if ("province" %in% names(annot)) {
  prov <- annot[tree$tip.label, "province"]
  place <- ifelse(!is.na(prov) & nzchar(prov), paste0(place, ": ", prov), place)
}
label <- ifelse(is_contig,
                paste0("c", annot[tree$tip.label, "cluster"], " ", place, " ", tree$tip.label),
                ifelse(tree$tip.label %in% names(ref_name),
                       paste0(ref_name[tree$tip.label], " (", tree$tip.label, ")"),
                       tree$tip.label))
# Name only our contigs and the RefSeq references, as 370 downloads drown the big tree out
is_refseq <- grepl("^NC_", tree$tip.label)
busy <- length(tree$tip.label) > 300 && !all_labels
if (busy) {
  label[!is_contig & !is_refseq] <- ""
  cat("labels hidden for", sum(!is_contig & !is_refseq),
      "downloaded references, use --all-labels to show them\n")
}
tipinfo <- data.frame(label = tree$tip.label, tiplab = label,
                      kind = ifelse(is_contig, "contig", "reference"),
                      stringsAsFactors = FALSE)

# Step 9: build one data frame per strip, rownames the tip labels so gheatmap lines them up
strip <- function(column, values = NULL) {
  v <- rep(NA_character_, length(tree$tip.label))
  names(v) <- tree$tip.label
  if (column %in% names(annot) && any(is_contig)) {
    value <- as.character(annot[match(tree$tip.label[is_contig], annot$contig), column])
    v[is_contig] <- value
  } else {
    cat("no", column, "column in the table, that strip is left blank\n")
  }
  if (!is.null(values)) v[!is.na(v) & !(v %in% values)] <- "other"
  data.frame(x = v, row.names = tree$tip.label)
}

country <- strip("country", names(COUNTRY_COLOURS))
host_genus <- strip("host_genus")
host_group <- strip("host_group", names(HOST_COLOURS))
province <- strip("province")
province$x[!is.na(province$x) & !nzchar(province$x)] <- NA
cluster <- strip("cluster")
cluster$x <- factor(cluster$x, levels = as.character(sort(unique(as.integer(cluster$x)))))
# Keep only the countries in the tree, so the legend is not full of empty entries
present <- unique(na.omit(country$x))
COUNTRY_COLOURS <- COUNTRY_COLOURS[names(COUNTRY_COLOURS) %in% present]
country_of <- country$x
cluster_colours <- setNames(rep(CLUSTER_PALETTE, length.out = nlevels(cluster$x)),
                            levels(cluster$x))
province_levels <- sort(unique(na.omit(province$x)))
province_colours <- setNames(rep(PROVINCE_PALETTE, length.out = length(province_levels)),
                             province_levels)

names(country) <- "country"
names(host_genus) <- "host genus"
names(host_group) <- "Culex/Aedes"
names(province) <- "province"
names(cluster) <- "cluster"

# Step 10: draw the tree, tip labels coloured by country so it reads without the strip
tipinfo$labelcol <- ifelse(is_contig, country_of, "reference")
label_colours <- c(COUNTRY_COLOURS, reference = "#8A8A8A")

p <- ggtree(tree, size = 0.3) %<+% tipinfo

span <- max(p$data$x, na.rm = TRUE)
colw <- 0.045                      # heatmap column width, a fraction of the tree width
gap <- span * 0.012
step <- span * colw * 1.25
# Place the strips right after the tips, the labels after the strips
lab_off <- gap + 4 * step + span * 0.01

p <- p +
  # Mark every one of our contigs, so they are findable among the references
  geom_tippoint(aes(colour = labelcol, subset = kind == "contig"), size = 0.5) +
  geom_tiplab(aes(label = tiplab, colour = labelcol), size = 1.6, offset = lab_off,
              align = TRUE, linesize = 0.08, linetype = "dotted") +
  scale_colour_manual(values = label_colours, na.value = "#8A8A8A", guide = "none") +
  new_scale_fill()

# Step 11: add the four strips left to right, gheatmap drawing one column of boxes each
add_strip <- function(p, df, scale_fill, offset) {
  gheatmap(p, df, offset = offset, width = colw, color = "white",
           colnames = TRUE, colnames_position = "bottom", colnames_angle = 90,
           hjust = 1, colnames_offset_y = -1, font.size = 2.5) +
    scale_fill +
    new_scale_fill()
}

p <- add_strip(p, country,
               scale_fill_manual(values = COUNTRY_COLOURS, na.value = "#F2F2F2",
                                 name = "country", na.translate = FALSE),
               gap)
p <- add_strip(p, province,
               scale_fill_manual(values = province_colours, na.value = "#F2F2F2",
                                 name = "province (China)", na.translate = FALSE),
               gap + step)
p <- add_strip(p, host_group,
               scale_fill_manual(values = HOST_COLOURS, na.value = "#F2F2F2",
                                 name = "Culex / Aedes", na.translate = FALSE),
               gap + 2 * step)
p <- gheatmap(p, cluster, offset = gap + 3 * step, width = colw, color = "white",
              colnames = TRUE, colnames_position = "bottom", colnames_angle = 90,
              hjust = 1, colnames_offset_y = -1, font.size = 2.5) +
  scale_fill_manual(values = cluster_colours, na.value = "#F2F2F2",
                    name = "CD-HIT cluster", na.translate = FALSE,
                    guide = guide_legend(ncol = 2, order = 3))

# Step 12: set the legend, title + theme
p <- p + theme_tree2() +
  xlim(0, span * 2.1) +
  theme(legend.position = "left",
        legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 10),
        plot.margin = margin(5, 5, 20, 5)) +
  ggtitle(if (china_only) "China contigs and RefSeq genomes, midpoint rooted"
          else "Clusters 0-27 contigs and RefSeq genomes, midpoint rooted")

# Step 13: save it, tall enough that the tip labels do not overlap
height <- max(6, 0.085 * length(tree$tip.label))
ggsave(paste0(out_prefix, "_ggtree.pdf"), p, width = 11, height = height, limitsize = FALSE)
ggsave(paste0(out_prefix, "_ggtree.png"), p, width = 11, height = height, dpi = 200,
       limitsize = FALSE, bg = "white")
cat("wrote", paste0(out_prefix, "_ggtree.pdf"), "and .png\n")
