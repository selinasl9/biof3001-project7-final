#!/usr/bin/env python3
"""
palette.py

Hold the figure colours in one place, so every figure colours the same thing the same way.

Used by: figures/plot_cluster_summary.py, figures/plot_host_vs_virus.py,
figures/render_virus_tree.py, tables/make_result_tables.py
"""

# Set the country, reference, host, cluster and virus colours + the relabelling rule
# The nine countries match Cammy's balanced beta analysis, and a contig 87-89% identical to
# its nearest named reference is not that virus, so clusters 4 and 26 are relabelled
COUNTRY_COLOURS = {
    "China": "#332288",
    "South Korea": "#CC6677",
    "Australia": "#117733",
    "Finland": "#DDCC77",
    "USA": "#88CCEE",
    "Egypt": "#882255",
    "Canada": "#999933",
    "Democratic Republic of the Congo": "#44AA99",
    "Viet Nam": "#AA4499",
}

OTHER_COUNTRY_COLOUR = "#333333"

REFERENCE_COLOUR = "#999999"

HOST_COLOURS = {"Culex": "#E69F00", "Aedes": "#009E73", "unknown": "#DDDDDD"}

CHINA_NO_PROVINCE_COLOUR = "#B0B0B0"

VIRUS_PALETTE = ["#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#469990",
                 "#F032E6", "#808000", "#42D4F4", "#9A6324", "#000075", "#BFEF45",
                 "#800000"]

RELABEL = [("Mosquito flavivirus", {"4", "26"}, "Unclassified ISF (Xizang lineage)")]
RELABEL_COLOURS = {"Unclassified ISF (Xizang lineage)": "#911EB4"}

NEW_SPECIES_DISTANCE = 0.15


def relabel_lineage(nearest_reference_name, cluster):
    """Return the lineage name for a contig, after the relabelling rule."""
    for reference, clusters, new_name in RELABEL:
        if nearest_reference_name == reference and str(cluster) in clusters:
            return new_name
    return nearest_reference_name


def country_colour(country):
    """Return the colour for a country, or the catch-all for anything outside the nine."""
    return COUNTRY_COLOURS.get(country, OTHER_COUNTRY_COLOUR)
