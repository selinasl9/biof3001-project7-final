# palette.R
#
# Hold the figure colours for the R figures, the same values as common/palette.py.
#
# Used by: figures/render_virus_tree_ggtree.R
#
# Set the country, reference, host, province, cluster and virus colours
# R cannot read the Python module, so a colour change has to be made in both files;
# the nine countries match Cammy's balanced beta analysis
COUNTRY_COLOURS <- c(
  "China" = "#332288",
  "South Korea" = "#CC6677",
  "Australia" = "#117733",
  "Finland" = "#DDCC77",
  "USA" = "#88CCEE",
  "Egypt" = "#882255",
  "Canada" = "#999933",
  "Democratic Republic of the Congo" = "#44AA99",
  "Viet Nam" = "#AA4499"
)

OTHER_COUNTRY_COLOUR <- "#333333"

REFERENCE_COLOUR <- "#999999"

HOST_COLOURS <- c("Culex" = "#E69F00", "Aedes" = "#009E73", "unknown" = "#DDDDDD")

CHINA_NO_PROVINCE_COLOUR <- "#B0B0B0"

VIRUS_PALETTE <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#469990",
                   "#F032E6", "#808000", "#42D4F4", "#9A6324", "#000075", "#BFEF45",
                   "#800000")

# Provinces are assigned in alphabetical order, so the figure is reproducible
PROVINCE_PALETTE <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948",
                      "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#1B9E77", "#D95F02",
                      "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

# Cluster colours recycle over the 28 clusters, as the number is on the tip label too
CLUSTER_PALETTE <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                     "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

# A contig 87-89% identical to its nearest named reference is not that virus,
# so clusters 4 and 26 are relabelled and get their own colour
RELABEL_FROM <- "Mosquito flavivirus"
RELABEL_CLUSTERS <- c("4", "26")
RELABEL_TO <- "Unclassified ISF (Xizang lineage)"
RELABEL_COLOURS <- c("Unclassified ISF (Xizang lineage)" = "#911EB4")

# Distance beyond which a contig is treated as a new species
NEW_SPECIES_DISTANCE <- 0.15


country_colour <- function(country) {
  # Colour for a country, or the catch-all for anything outside the nine
  ifelse(country %in% names(COUNTRY_COLOURS),
         COUNTRY_COLOURS[country], OTHER_COUNTRY_COLOUR)
}


relabel_lineage <- function(nearest_reference_name, cluster) {
  # The lineage name for a contig, after the relabelling rule
  ifelse(nearest_reference_name == RELABEL_FROM &
           as.character(cluster) %in% RELABEL_CLUSTERS,
         RELABEL_TO, nearest_reference_name)
}


recycle_palette <- function(palette, levels) {
  # One colour per level, recycled if there are more levels than colours
  setNames(rep(palette, length.out = length(levels)), levels)
}
