#!/usr/bin/env python3
"""
plot_cluster_summary.py

Draw how far each cluster sits from a known virus, and where its contigs come from.

Inputs:
    ANNOTATION   contig_annotation.tsv, one row per contig
    NEAREST      cluster_nearest_neighbour.tsv, cluster, nearest_accession, distance

Outputs:
    OUTDIR/novelty_before_after.png and .pdf
    OUTDIR/cluster_by_place.png and .pdf
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.append(str(Path(__file__).resolve().parents[1]))

from common.palette import CHINA_NO_PROVINCE_COLOUR, NEW_SPECIES_DISTANCE

# Set the place fix-up + the figure colours
# Dongcheng is a Beijing district, so the BioSample names it instead of a province
DISTRICT_TO_PROVINCE = {"Dongcheng District": "Beijing"}

DIVERGENT_COLOUR = "#D7263D"
CLOSE_COLOUR = "#2E86AB"
OLD_POINT_EDGE = "#8A8A8A"
LINK_COLOUR = "#C9C9C9"


def read_contigs(path):
    """Return the annotation rows, with province blanked outside China."""
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for row in rows:
        province = (row.get("province") or "").strip()
        province = DISTRICT_TO_PROVINCE.get(province, province)
        row["province"] = province if row["country"] == "China" else ""
    return rows


def place_of(row):
    """Return the place a contig is counted under: a Chinese province, or a country."""
    if row["country"] != "China":
        return row["country"]
    return f"China: {row['province']}" if row["province"] else "China: province n/a"


def place_colours(rows):
    """Return a colour-per-row function, provinces on one colormap and countries on another."""
    provinces = sorted({row["province"] for row in rows if row["province"]})
    countries = sorted({row["country"] for row in rows if row["country"] != "China"})

    turbo, tab10 = plt.get_cmap("turbo"), plt.get_cmap("tab10")
    province_colour = {name: turbo(i / max(len(provinces) - 1, 1))
                       for i, name in enumerate(provinces)}
    country_colour = {name: tab10(i % 10) for i, name in enumerate(countries)}

    def colour_of(row):
        if row["country"] != "China":
            return country_colour[row["country"]]
        return province_colour.get(row["province"], CHINA_NO_PROVINCE_COLOUR)

    return colour_of


def plot_novelty(by_cluster, nearest, outdir):
    """Draw the distance to the nearest reference, before and after the GenBank genomes."""
    # Step 1: take the old distance as the closest any contig sat in the RefSeq-only tree
    before = {}
    for cluster, contigs in by_cluster.items():
        distances = []
        for contig in contigs:
            try:
                distances.append(float(contig["tree_distance"]))
            except (ValueError, KeyError):
                pass
        if distances:
            before[cluster] = min(distances)

    # Step 2: order the clusters, most divergent at the top
    order = sorted(by_cluster, key=lambda cluster: -nearest[cluster][1])

    # Step 3: draw one row per cluster, a line from the old distance to the new one
    figure, axes = plt.subplots(figsize=(11, 7))
    for i, cluster in enumerate(order):
        after = nearest[cluster][1]

        if cluster in before:
            axes.plot([before[cluster], after], [i, i], color=LINK_COLOUR, lw=2.2,
                      zorder=1, solid_capstyle="round")
            axes.plot([before[cluster]], [i], "o", ms=7, mfc="white",
                      mec=OLD_POINT_EDGE, mew=1.4, zorder=2)

        axes.plot([after], [i], "o", ms=8,
                  mfc=DIVERGENT_COLOUR if after > NEW_SPECIES_DISTANCE else CLOSE_COLOUR,
                  mec="white", mew=0.8, zorder=3)
        axes.text(-0.012, i, f"cluster {cluster}  (n={len(by_cluster[cluster])})",
                  ha="right", va="center", fontsize=8)

        furthest = max(after, before.get(cluster, after))
        axes.text(furthest + 0.014, i, nearest[cluster][0], ha="left", va="center",
                  fontsize=7, color="#555")

    # Step 4: mark the species cut-off, so a reader can see which clusters are still novel
    axes.axvline(NEW_SPECIES_DISTANCE, color=DIVERGENT_COLOUR, ls=(0, (4, 3)), lw=1.0, zorder=0)
    axes.text(NEW_SPECIES_DISTANCE + 0.005, len(order) - 0.4,
              "roughly where a new species starts", fontsize=7.5, color=DIVERGENT_COLOUR)

    # Step 5: set the axes, legend + title, with the cluster names standing in for y ticks
    axes.set_ylim(-1, len(order))
    axes.set_yticks([])
    axes.invert_yaxis()
    axes.set_xlim(-0.17, 0.60)
    axes.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    axes.set_xlabel("patristic distance to the nearest reference (substitutions per site)",
                    fontsize=9)
    for side in ("top", "right", "left"):
        axes.spines[side].set_visible(False)

    legend = [
        plt.Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec=OLD_POINT_EDGE,
                   mew=1.4, label="nearest RefSeq virus, 91-reference tree"),
        plt.Line2D([], [], ls="", marker="o", ms=8, mfc=CLOSE_COLOUR, mec="white",
                   label="nearest genome after adding 370 GenBank references"),
        plt.Line2D([], [], ls="", marker="o", ms=8, mfc=DIVERGENT_COLOUR, mec="white",
                   label=f"still divergent (> {NEW_SPECIES_DISTANCE})"),
    ]
    axes.legend(handles=legend, loc="lower right", bbox_to_anchor=(1.0, 0.03),
                fontsize=8, frameon=False)
    axes.set_title("Adding GenBank genomes removes most of the apparent novelty",
                   fontsize=12, loc="left")

    # Step 6: save a PNG for the slides + a PDF for the report
    figure.tight_layout()
    figure.savefig(f"{outdir}/novelty_before_after.png", dpi=200)
    figure.savefig(f"{outdir}/novelty_before_after.pdf")


def plot_cluster_by_place(by_cluster, colour_of, outdir):
    """Draw one stacked bar per cluster, split by where its contigs came from."""
    # Step 1: order the clusters, biggest at the top
    order = sorted(by_cluster, key=lambda cluster: (-len(by_cluster[cluster]), cluster))

    # Step 2: draw one bar per cluster, one segment per place, collecting the legend as we go
    figure, axes = plt.subplots(figsize=(11, 7.5))
    seen = {}

    for i, cluster in enumerate(order):
        contigs = by_cluster[cluster]
        counts = Counter(place_of(contig) for contig in contigs)
        colours = {place_of(contig): colour_of(contig) for contig in contigs}

        left = 0
        for place, count in counts.most_common():
            axes.barh(i, count, left=left, color=colours[place], edgecolor="white",
                      lw=0.8, height=0.78)
            seen[place] = colours[place]
            left += count

        # Call out how many places a cluster spans, as a single-place cluster is the point
        n_places = len(counts)
        axes.text(left + 0.3, i, f"{n_places} place{'s' if n_places > 1 else ''}",
                  va="center", fontsize=7.5,
                  color="#333" if n_places > 1 else "#2E7D32",
                  fontweight="normal" if n_places > 1 else "bold")
        axes.text(-0.4, i, f"cluster {cluster}", ha="right", va="center", fontsize=8)

    # Step 3: set the axes, the note on what a place is, legend + title
    axes.set_ylim(-1, len(order))
    axes.invert_yaxis()
    axes.set_yticks([])
    axes.set_xlabel("contigs", fontsize=9)
    axes.text(0, len(order) + 0.4,
              "a place is a Chinese province, or a country outside China; "
              "contigs with no province recorded are pooled as \"province n/a\"",
              fontsize=7.5, color="#666", va="top")
    axes.set_xlim(-6, 30)
    axes.set_xticks([0, 5, 10, 15, 20, 25, 30])
    for side in ("top", "right", "left"):
        axes.spines[side].set_visible(False)

    axes.legend(handles=[Patch(facecolor=colour, label=place)
                         for place, colour in sorted(seen.items())],
                loc="lower right", fontsize=7.2, frameon=False, ncol=2,
                title="province / country", title_fontsize=8)
    axes.set_title("Almost every cluster comes from one place\n"
                   "26 of 28 clusters are single-country (Cramer's V 0.85); "
                   "10 of the 16 clusters with province data are single-province (V 0.69)",
                   fontsize=11.5, loc="left")

    # Step 4: save a PNG for the slides + a PDF for the report
    figure.tight_layout()
    figure.savefig(f"{outdir}/cluster_by_place.png", dpi=200)
    figure.savefig(f"{outdir}/cluster_by_place.pdf")


def main():
    annotation, nearest_path, outdir = sys.argv[1], sys.argv[2], sys.argv[3]

    # Step 1: read the contigs, with the province tidied up
    rows = read_contigs(annotation)

    # Step 2: read the nearest reference + its distance, one row per cluster
    nearest = {int(row["cluster"]): (row["nearest_accession"],
                                     float(row["patristic_distance"]))
               for row in csv.DictReader(open(nearest_path), delimiter="\t")}

    # Step 3: group the contigs by cluster, which is what both figures draw
    by_cluster = defaultdict(list)
    for row in rows:
        by_cluster[int(row["cluster"])].append(row)

    # Step 4: draw both figures
    plot_novelty(by_cluster, nearest, outdir)
    plot_cluster_by_place(by_cluster, place_colours(rows), outdir)
    print("wrote novelty_before_after and cluster_by_place")


if __name__ == "__main__":
    main()
