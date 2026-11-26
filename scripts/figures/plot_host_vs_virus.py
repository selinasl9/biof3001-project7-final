#!/usr/bin/env python3
"""
plot_host_vs_virus.py

Draw the host COI tree, with the flaviviruses each mosquito species carries beside it.

Inputs:
    TREE          host_coi.treefile, the barcode-trimmed host COI tree
    ANNOTATION    contig_annotation.tsv, for the carriage counts
    COLOURS       optional virus_colours.tsv, from render_virus_tree.py --dump-colours

Outputs:
    OUT.pdf       the tree + stacked bars
"""

import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, __file__.rsplit("/", 2)[0])
from common import newick
from common.palette import RELABEL_COLOURS, VIRUS_PALETTE, relabel_lineage

# Set the tips to drop, the support cutoff + the figure's fixed text
# Drosophila is the outgroup, and it lands oddly on COI alone
DROP_TIPS = {"Drosophila melanogaster"}

MIN_SUPPORT = 50

BRANCH_COLOUR = "#2B2B2B"
TITLE = "Mosquito phylogeny (COI) and flavivirus distribution across host species"


def read_host_tree(path):
    """Return the COI tree with the outgroup dropped and the branches ladderised."""
    root = prune(newick.parse(path))
    newick.link_parents(root)
    newick.ladderise(root)
    return root


def prune(node):
    """Drop the tips in DROP_TIPS, collapsing any stem that is left behind."""
    if not node["children"]:
        # Tidy the tip name, as the COI tips are "Culex_pipiens|KX..." with _ for spaces
        node["name"] = node["name"].split("|")[0].replace("_", " ")
        return None if node["name"] in DROP_TIPS else node

    children = [kept for kept in (prune(child) for child in node["children"]) if kept]
    if not children:
        return None
    if len(children) == 1:
        children[0]["length"] += node["length"]
        return children[0]

    node["children"] = children
    return node


def read_carriage(path):
    """Return {host species: {virus lineage: contigs}} from the annotation table."""
    carriage = defaultdict(lambda: defaultdict(int))
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            species = (row.get("host_species") or "").strip()
            if not row.get("contig") or not species:
                continue
            lineage = relabel_lineage(row["nearest_reference_name"], row["cluster"])
            carriage[species][lineage] += 1
    return carriage


def virus_colours(carriage, path=None):
    """Return the lineages + a colour each, with the big tree's own colours taking over."""
    # Sort by contigs then name, as ties would otherwise re-order and move colours between runs
    lineages = sorted({v for species in carriage for v in carriage[species]},
                      key=lambda v: (-sum(carriage[s].get(v, 0) for s in carriage), v))
    colours = {v: VIRUS_PALETTE[i % len(VIRUS_PALETTE)] for i, v in enumerate(lineages)}

    if path:
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["virus"] in colours:
                    colours[row["virus"]] = row["colour"]
    colours.update({k: v for k, v in RELABEL_COLOURS.items() if k in colours})
    return lineages, colours


def place(root):
    """Give every node an x (distance from the root) and a y, and return the tips in order."""
    order = []

    def walk(node, depth):
        node["x"] = depth + node["length"]
        if not node["children"]:
            node["y"] = len(order)
            order.append(node)
            return
        for child in node["children"]:
            walk(child, node["x"])
        node["y"] = sum(child["y"] for child in node["children"]) / len(node["children"])

    walk(root, 0.0)
    return order


def draw_tree(ax, node, span):
    """Draw the tree itself, with a support value above each node that has one."""
    for child in node["children"]:
        ax.plot([node["x"], child["x"]], [child["y"], child["y"]],
                color=BRANCH_COLOUR, lw=1.4, zorder=2)
        draw_tree(ax, child, span)

    if not node["children"]:
        return
    ax.plot([node["x"], node["x"]],
            [node["children"][0]["y"], node["children"][-1]["y"]],
            color=BRANCH_COLOUR, lw=1.4, zorder=2)

    support = node.get("support", "")
    if support.replace(".", "").isdigit() and float(support) >= MIN_SUPPORT:
        ax.text(node["x"] - span * 0.008, node["y"] + 0.14, support, ha="right",
                va="bottom", fontsize=6.5, color="#777777")


def draw_bars(ax, order, carriage, lineages, colours, span):
    """Draw a species name, a dotted leader + one stacked bar per host species."""
    name_x = span * 1.06
    bar_x = name_x + span * 0.95
    bar_width = span * 1.5
    most = max(sum(carriage.get(tip["name"], {}).values()) for tip in order) or 1

    for tip in order:
        y = tip["y"]
        ax.plot([tip["x"], name_x - span * 0.01], [y, y], color="#D5D5D5", lw=0.4, zorder=1)
        ax.text(name_x, y, tip["name"], va="center", fontsize=8.5, style="italic")

        x = bar_x
        for lineage in lineages:
            contigs = carriage.get(tip["name"], {}).get(lineage, 0)
            if not contigs:
                continue
            width = bar_width * contigs / most
            ax.add_patch(Rectangle((x, y - 0.3), width, 0.6, facecolor=colours[lineage],
                                   edgecolor="white", lw=0.5, zorder=3))
            # Write the count only where the segment is wide enough to hold it
            if width > span * 0.09:
                ax.text(x + width / 2, y, str(contigs), ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold", zorder=4)
            x += width

        total = sum(carriage.get(tip["name"], {}).values())
        if total:
            ax.text(x + span * 0.02, y, f"n={total}", va="center", fontsize=7, color="#555555")

    return bar_x + bar_width


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    tree_path, annotation_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
    colour_path = sys.argv[4] if len(sys.argv) > 4 else None

    # Step 1: read the COI tree, outgroup dropped and ladderised
    root = read_host_tree(tree_path)

    # Step 2: count how many contigs of each virus lineage each host species carries
    carriage = read_carriage(annotation_path)

    # Step 3: pick a colour per lineage, matching the big tree where it has one
    lineages, colours = virus_colours(carriage, colour_path)

    # Step 4: place every node + size a figure tall enough for the tips
    order = place(root)
    span = max(tip["x"] for tip in order)
    tips = len(order)
    fig, ax = plt.subplots(figsize=(11, max(3.2, 0.52 * tips) + 1.6))

    # Step 5: draw the tree on the left + the stacked bars on the right
    draw_tree(ax, root, span)
    right = draw_bars(ax, order, carriage, lineages, colours, span)

    # Step 6: add the scale bar, title + legend, with the axes themselves hidden
    ax.set_xlim(-span * 0.03, right + span * 0.22)
    ax.set_ylim(-1.4, tips - 0.2)
    ax.invert_yaxis()
    ax.axis("off")

    scale = span * 0.1
    ax.plot([0, scale], [tips - 0.6, tips - 0.6], color=BRANCH_COLOUR, lw=1.2)
    ax.text(scale / 2, tips - 0.45, f"{scale:.2f} substitutions per site (COI)",
            ha="center", va="top", fontsize=7)
    ax.text(0, -1.1, TITLE, fontsize=11, va="bottom")
    ax.legend(handles=[Rectangle((0, 0), 1, 1, facecolor=colours[v], label=v)
                       for v in lineages],
              loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=4, frameon=False,
              fontsize=7.5, handlelength=1.1, handleheight=1.1)

    # Step 7: save it, leaving room under the axes for the legend
    fig.subplots_adjust(bottom=0.16)
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out, f"({tips} host species, {len(lineages)} viruses)")


if __name__ == "__main__":
    main()
