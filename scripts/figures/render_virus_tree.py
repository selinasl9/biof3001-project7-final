#!/usr/bin/env python3
"""
render_virus_tree.py

Draw the annotated ISF tree, with clade blocks, coloured branches + heatmap strips.

Inputs:
    TREE           the virus tree .treefile, for topology and branch lengths
    ANNOTATION     contig_annotation.tsv, cluster, country, province, host, length, pident
    REFS           optional reference FASTA, for the reference names
    --ref-meta     reference_country.tsv, to colour the reference tips too

Outputs:
    OUT.pdf        the tree, plus a .png beside it
    --dump-colours virus -> colour, read by plot_host_vs_virus.py

Flags:
    --prune                drop the reference tips that carry no contigs
    --china single|province|region   how China is broken up in the place column
    --by cluster|place|clade|clade+place   what the branch colours mean
    --simple-heat          two heatmap columns instead of four
    --no-bands             no shaded bands behind the rows
    --split-at D           branch length at which a clade is cut (default 0.06)
    --clades N             keep only the N largest clades
    --flat-clades          one tone per clade instead of shading within it
    --sisters              label the sister group of each clade
    --collapse-refs        collapse a run of reference tips into one triangle
    --virus-palette        colour by virus name rather than by clade
    --dump-colours FILE    write the virus -> colour table
    --zoom A,B,C           draw only these clades
    --context N            tips of context to keep around a zoom (default 2)
    --slide                slide proportions, bigger text, fewer labels
    --height H             figure height in inches
    --min-label N          smallest clade that gets a label (default 3)

"""

import csv, re, sys
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

sys.path.insert(0, __file__.rsplit("/", 2)[0])
from common import newick

TREE, ANNOT, OUT, TITLE = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
REFS = sys.argv[5] if len(sys.argv) > 5 else None

# Step 1: read the tree, re-root it + put the tips in a sensible order
# The parser, the tip walk and the ladderising live in common/newick.py
root = newick.parse(TREE)
# Strip _R_, which MAFFT adds to anything it reverse complemented, not part of the name
for tip in newick.tips(root):
    tip["name"] = tip["name"].replace("_R_", "")

tips = newick.tips
path_to_root = newick.path_to_root
dist = newick.patristic_distance

# Midpoint root here rather than with newick.midpoint_root, as the cut point is also
# what the clade cutting below measures from
all_tips = tips(root)
a = max(all_tips, key=lambda t: dist(all_tips[0], t))
b = max(all_tips, key=lambda t: dist(a, t))
total = dist(a, b)
def find_cut(start, other):
    oth = set(id(x) for x in path_to_root(other)); walked = 0.0
    for n in path_to_root(start):
        if id(n) in oth: return None, 0.0
        if walked + n["length"] >= total / 2: return n, total / 2 - walked
        walked += n["length"]
    return None, 0.0
cut, remain = find_cut(a, b)
if cut is None: cut, remain = find_cut(b, a)

def reroot(cut, below):
    chain = []; n = cut
    while n["parent"] is not None: chain.append(n); n = n["parent"]
    old_root = n
    for child, par in zip(chain, chain[1:] + [old_root]):
        par["children"].remove(child); child["children"].append(par); par["parent"] = child
    above = cut["children"].pop()
    new = {"children": [cut, above], "name": "", "support": "", "length": 0.0, "parent": None}
    cut["length"], above["length"] = below, max(cut["length"] - below, 0.0)
    cut["parent"] = new; above["parent"] = new
    return new
root = reroot(cut, remain)
newick.collapse_single_children(root)
newick.ladderise(root)

# Set the province to region lookup, used by --china region and the place labels
REGION = {"Beijing": "North", "Tianjin": "North", "Hebei": "North", "Shanxi": "North",
          "Inner Mongolia": "North",
          "Heilongjiang": "Northeast", "Jilin": "Northeast", "Liaoning": "Northeast",
          "Shanghai": "East", "Jiangsu": "East", "Zhejiang": "East", "Anhui": "East",
          "Fujian": "East", "Jiangxi": "East", "Shandong": "East",
          "Henan": "Central", "Hubei": "Central", "Hunan": "Central",
          "Guangdong": "South", "Guangxi": "South", "Hainan": "South",
          "Chongqing": "Southwest", "Sichuan": "Southwest", "Guizhou": "Southwest",
          "Yunnan": "Southwest", "Xizang": "Southwest",
          "Shaanxi": "Northwest", "Gansu": "Northwest", "Qinghai": "Northwest",
          "Ningxia": "Northwest", "Xinjiang": "Northwest"}

# Step 2: read the annotation table, tidy the provinces + work out which tips to keep
# Some BioSamples give a site or a city where a province belongs
SITE_TO_PROVINCE = {"xiangyun": "Yunnan", "dali": "Yunnan", "kunming": "Yunnan",
                    "jiujiang": "Jiangxi", "ganzhou": "Jiangxi", "nanchang": "Jiangxi",
                    "dongcheng": "Beijing", "xicheng": "Beijing", "chaoyang": "Beijing",
                    "haidian": "Beijing", "fengtai": "Beijing",
                    "pudong": "Shanghai", "minhang": "Shanghai", "jiading": "Shanghai",
                    "wuhan": "Hubei", "qingdao": "Shandong", "harbin": "Heilongjiang",
                    "sanya": "Hainan", "haikou": "Hainan"}

def norm_province(prov, country):
    prov = (prov or "").strip()
    if country != "China" or not prov: return ""
    if prov in REGION: return prov
    low = prov.lower()
    for key, province in SITE_TO_PROVINCE.items():
        if key in low: return province
    # Strip the suffix, as in "Xizang Autonomous Region" or "Beijing City"
    short = re.sub(r"\s+(Province|City|Municipality|Autonomous Region|District|County)$",
                   "", prov).strip()
    return short if short in REGION else prov

annot = {}
with open(ANNOT) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        r["contig"] = r["contig"].replace("_R_", "")
        r["province"] = norm_province(r.get("province"), r.get("country"))
        annot[r["contig"]] = r

REFMETA = sys.argv[sys.argv.index("--ref-meta") + 1] if "--ref-meta" in sys.argv else None
refcountry, refhost = {}, {}
if REFMETA:
    import csv as _csv
    for r in _csv.DictReader(open(REFMETA), delimiter="\t"):
        if r.get("country"): refcountry[r["accession"]] = r["country"]
        if r.get("host"): refhost[r["accession"]] = r["host"]
    print(f"{len(refcountry)} references have a country")

refname = {}
with open(ANNOT) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        acc, nm = (r.get("nearest_reference") or "").strip(), (r.get("nearest_reference_name") or "").strip()
        if acc and nm: refname.setdefault(acc, nm)
if REFS:
    for path_ in REFS.split(","):
        try: lines = open(path_)
        except FileNotFoundError: continue
        for line in lines:
            if line.startswith(">"):
                acc, *rest = line[1:].split()
                nm = " ".join(rest).split(",")[0]
                nm = re.split(r"\s+(?:strain|isolate|complete|genomic|polyprotein|gene)\b", nm)[0]
                refname.setdefault(acc, nm.strip())

ZOOM = []
if "--zoom" in sys.argv:
    ZOOM = [int(z) for z in sys.argv[sys.argv.index("--zoom") + 1].split(",")]
CONTEXT = int(sys.argv[sys.argv.index("--context") + 1]) if "--context" in sys.argv else 2
if ZOOM:
    want = [t for t in tips(root) if t["name"] in annot
            and int(annot[t["name"]]["cluster"]) in ZOOM]
    paths = [set(id(x) for x in path_to_root(w)) for w in want]
    mrca = None
    for nd in path_to_root(want[0]):
        if all(id(nd) in pth for pth in paths): mrca = nd; break
    for _ in range(CONTEXT):
        if mrca["parent"] is not None: mrca = mrca["parent"]
    root = mrca; root["parent"] = None; root["length"] = 0.0
    all_tips = tips(root)

PRUNE = "--prune" in sys.argv
if PRUNE:
    def nearest_set():
        keep = set()
        for t in tips(root):
            if t["name"] in annot:
                keep.add(t["name"])
        refs = [t for t in tips(root) if t["name"] not in keep]
        for t in tips(root):
            if t["name"] in annot:
                d, who = min(((dist(t, r), r["name"]) for r in refs), key=lambda z: z[0])
                keep.add(who)
        return keep
    KEEP = nearest_set()
    def prune(nd):
        nd["children"] = [k for k in nd["children"] if prune(k)]
        if not nd["children"]:
            return nd["name"] in KEEP
        return True
    prune(root)
    newick.collapse_single_children(root)
    while len(root["children"]) == 1:
        root = root["children"][0]; root["parent"] = None
    newick.ladderise(root)

COLLAPSE = "--collapse-refs" in sys.argv
if COLLAPSE:
    def tagc(nd):
        if not nd["children"]:
            nd["nc"] = 1 if nd["name"] in annot else 0
            return nd["nc"]
        nd["nc"] = sum(tagc(k) for k in nd["children"])
        return nd["nc"]
    tagc(root)
    def mark(nd):
        for k in nd["children"]:
            if k["nc"] == 0 and k["children"]:
                k["collapsed"] = True
                k["nref"] = len(tips(k))
            else:
                mark(k)
    mark(root)

order = []
def deepest(nd, d):
    d += nd["length"]
    return d if not nd["children"] else max(deepest(k, d) for k in nd["children"])

def assign(n, depth):
    n["x"] = depth + n["length"]
    if n.get("collapsed"):
        n["xmax"] = deepest(n, depth)
        n["y"] = len(order); order.append(n)
    elif not n["children"]:
        n["y"] = len(order); order.append(n)
    else:
        for k in n["children"]: assign(k, n["x"])
        n["y"] = sum(k["y"] for k in n["children"]) / len(n["children"])
assign(root, 0.0)
span = max(max(t["x"] for t in order), max((t.get("xmax", 0) for t in order), default=0))
n = len(order)
if n <= 45:   FS_TIP, FS_CL, FS_HDR, MS, LW = 8.5, 11.0, 9.0, 8.0, 2.4
elif n <= 110: FS_TIP, FS_CL, FS_HDR, MS, LW = 6.0, 8.5, 7.0, 6.0, 1.6
else:          FS_TIP, FS_CL, FS_HDR, MS, LW = 3.6, 5.4, 4.4, 3.0, 1.0

# --slide gives an overview for a slide, no tip names, only the bigger clades labelled
SLIDE = "--slide" in sys.argv
HEIGHT = float(sys.argv[sys.argv.index("--height") + 1]) if "--height" in sys.argv else None
MINLAB = int(sys.argv[sys.argv.index("--min-label") + 1]) if "--min-label" in sys.argv else 3
if SLIDE:
    FS_CL, FS_HDR, MS, LW = 7.0, 7.0, 2.6, 1.1

# Measure every contig to its nearest reference tip, on this tree
reftips = [t for t in order if t["name"] not in annot]
nearest = {}
for t in order:
    if t["name"] in annot:
        if reftips:
            d, who = min(((dist(t, r), r["name"]) for r in reftips), key=lambda z: z[0])
        else:
            d, who = None, ""
        nearest[t["name"]] = (d, who)

# Step 3: set the colours + the heatmap columns, following Cammy's figures
# --china province | region | single sets how much detail the China tips carry
CHINA = sys.argv[sys.argv.index("--china") + 1] if "--china" in sys.argv else "province"
TITLE_P = {"single": "China / other countries",
           "region": "Country / province",
           "province": "Country / province"}[CHINA]

REGION_ABBR = {"North": "N", "Northeast": "NE", "East": "E", "Central": "C",
               "South": "S", "Southwest": "SW", "Northwest": "NW"}

def short_place(r):
    """Return "Yunnan (SW)" for a Chinese province, the country for anywhere else."""
    if r.get("country") != "China":
        return r.get("country", "")
    pv = (r.get("province") or "").strip()
    if not pv: return "China"
    ab = REGION_ABBR.get(REGION.get(pv, ""), "")
    return f"{pv} ({ab})" if ab else pv

REGION_ORDER = ["North", "Northeast", "East", "Central", "South", "Southwest", "Northwest"]
# Colour China warm and everywhere else cool, one hue per region shaded per province,
# so the two never read as the same family and neighbouring provinces still look related
REGION_COL = {"Southwest": "#E6194B",   # red
              "South":     "#F58231",   # orange
              "East":      "#FFE119",   # yellow
              "Central":   "#9A6324",   # brown
              "Northeast": "#808000",   # olive
              "North":     "#F032E6",   # magenta
              "Northwest": "#F032E6"}
COOL_PAL = ["#4363D8",   # blue
            "#3CB44B",   # green
            "#911EB4",   # purple
            "#42D4F4",   # cyan
            "#000075",   # navy
            "#469990",   # teal
            "#BFEF45"]   # lime

def vivid(c):
    """Saturate a colour, as a band reads better saturated than dark."""
    h, sat, v = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
    return matplotlib.colors.hsv_to_rgb((h, min(1.0, sat * 1.45), min(1.0, v * 1.35)))

def shade(hexcol, i, n):
    """Return the i-th of n variants of one region colour, hue nudged + lightness stepped."""
    h, sat, v = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(hexcol))
    if n > 1:
        spread = 0.075
        h = (h + spread * (i - (n - 1) / 2) / max((n - 1) / 2, 1)) % 1.0
    # Walk the lightness from well below the base colour to well above it
    f = 1.0 if n == 1 else 0.55 + 0.90 * i / (n - 1)
    if f >= 1:
        v = min(1.0, v * (1 + (f - 1) * 0.30))
        sat = max(0.30, sat * (1 - (f - 1) * 0.70))
    else:
        v = max(0.32, v * (1 - (1 - f) * 0.95))
        sat = min(1.0, sat * (1 + (1 - f) * 0.30))
    return matplotlib.colors.hsv_to_rgb((h, sat, v))

provs = sorted({r["province"] for r in annot.values() if r["province"]})
# Shade provinces within their own region, biggest first so the commonest keeps the base
by_region = {}
for pv in provs:
    by_region.setdefault(REGION.get(pv, "Northwest"), []).append(pv)
count = Counter(r["province"] for r in annot.values() if r["province"])
prov_col = {}
for reg, members in by_region.items():
    base = REGION_COL.get(reg, "#9E9D24")
    # Shade darkest for the commonest province of the region, lightest for the rarest
    ordered = sorted(members, key=lambda x: (-count[x], x))
    for i, pv in enumerate(ordered):
        prov_col[pv] = shade(base, i, len(ordered))

others = sorted({r["country"] for r in annot.values() if r["country"] != "China"})
ctry_col = {c: COOL_PAL[i % len(COOL_PAL)] for i, c in enumerate(others)}
CN_ONE = "#C62828"
CN_NA = "#B0B0B0"

def place_label(r):
    if r["country"] != "China":
        return r["country"]
    if CHINA == "single":
        return "China"
    if CHINA == "region":
        return "China: " + REGION.get(r["province"], "region n/a")
    return "China: " + r["province"] if r["province"] else "China: province n/a"

def place_colour(r):
    if r["country"] != "China":
        return ctry_col[r["country"]]
    if CHINA == "single":
        return CN_ONE
    if CHINA == "region":
        return REGION_COL.get(REGION.get(r["province"], ""), CN_NA)
    return prov_col.get(r["province"], CN_NA)

def place_order():
    if CHINA == "single":
        keys = ["China"]
    elif CHINA == "region":
        keys = ["China: " + g for g in REGION_ORDER] + ["China: region n/a"]
    else:
        keys = ["China: " + p for p in provs] + ["China: province n/a"]
    return keys, list(others)

HOST = {"Culex": "#1F78B4", "Aedes": "#E31A1C", "unknown": "#FFFFFF"}
# Give the odd genus outside Culex / Aedes a marker of its own
for _g in sorted({(r.get("host_group") or "unknown") for r in annot.values()}):
    HOST.setdefault(_g, {"Coquillettidia": "#33A02C", "Haemagogus": "#6A3D9A"}.get(_g, "#B15928"))
CLUST_PAL = ["#1F78B4", "#E31A1C", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928",
             "#A6CEE3", "#FB9A99", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#D4A017",
             "#17BECF", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#2CA02C",
             "#D62728", "#9467BD", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
             "#66A61E", "#E6AB02", "#A6761D", "#000075"]
GREY = "#2B2B2B"
def dark(c):
    r, g, b = matplotlib.colors.to_rgb(c)
    return (r * 0.70, g * 0.70, b * 0.70)
all_clusters = sorted({int(r["cluster"]) for r in annot.values()})
clust_col = {c: CLUST_PAL[i % len(CLUST_PAL)] for i, c in enumerate(all_clusters)}
clusters = sorted({int(annot[t["name"]]["cluster"]) for t in order if t["name"] in annot})

def num(r, key):
    try: return float(r.get(key) or "")
    except ValueError: return None
def year(r):
    m = re.match(r"(\d{4})", r.get("collection_date") or "")
    return float(m.group(1)) if m else None

from matplotlib.colors import LinearSegmentedColormap
WHITE_TEAL = LinearSegmentedColormap.from_list("wt", ["#FFFFFF", "#17627A"])
WHITE_RUST = LinearSegmentedColormap.from_list("wr", ["#FFFFFF", "#A33B20"])
SIMPLE = "--simple-heat" in sys.argv
HEAT = [("genome length (kb)", "viridis", lambda r, c: (num(r, "length_bp") or 0) / 1000 or None),
        ("% identity to\nnearest virus", "plasma", lambda r, c: num(r, "pident")),
        ("distance to nearest\nreference (this tree)", "magma_r", lambda r, c: nearest[c][0]),
        ("collection year", "cividis", lambda r, c: year(r))]
if SIMPLE:
    HEAT = [("genome length (kb)", WHITE_TEAL, lambda r, c: (num(r, "length_bp") or 0) / 1000 or None),
            ("collection year", WHITE_RUST, lambda r, c: year(r))]
heat_vals = []
for title_, cmapname, fn in HEAT:
    vs = [fn(annot[t["name"]], t["name"]) for t in order if t["name"] in annot]
    vs = [v for v in vs if v is not None]
    heat_vals.append((min(vs), max(vs)) if vs else (0.0, 1.0))

# Step 4: cut the tree into clades at --split-at, so the blocks come from the tree itself
# Keep splitting whichever group reaches deepest, as you would read genotypes off by eye
NCLADE = int(sys.argv[sys.argv.index("--clades") + 1]) if "--clades" in sys.argv else 99
# Split a clade into sister clades only when it is this divergent inside
SPLIT_AT = float(sys.argv[sys.argv.index("--split-at") + 1]) if "--split-at" in sys.argv else 0.06

def spread(n):
    """Return how far the deepest tip sits below this node, in substitutions per site."""
    return max(t["x"] for t in tips(n)) - n["x"]

def cut_clades(k):
    """Split where the tree branches, not by how many tips a group holds."""
    parts = [root]
    while len(parts) < k:
        cand = [p for p in parts if p["children"] and spread(p) > SPLIT_AT]
        if not cand: break
        pick = max(cand, key=spread)
        parts.remove(pick)
        parts.extend(pick["children"])
    parts.sort(key=lambda p: min(t["y"] for t in tips(p)))
    return parts

clade_of = {}
# Order the palette so clades next to each other never share a hue family
CLADE_PAL = ["#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#469990",
             "#F032E6", "#808000", "#42D4F4", "#9A6324", "#000075", "#BFEF45",
             "#800000", "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E7298A",
             "#A6761D", "#386CB0", "#F0027F", "#BF5B17", "#66C2A5", "#FC8D62"]

def tone(hexcol, j):
    """Return the j-th variant of one colour, so sister groups separate without a new hue."""
    if j == 0: return hexcol
    h, sat, v = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(hexcol))
    steps = [(1.00, 1.00), (0.42, 1.30), (1.15, 0.58), (0.65, 1.15), (1.25, 0.78),
             (0.30, 1.40), (1.30, 0.44)]
    fs, fv = steps[j % len(steps)]
    return matplotlib.colors.hsv_to_rgb((h, min(1.0, max(0.18, sat * fs)),
                                         min(1.0, max(0.25, v * fv))))

# --virus-palette gives a clade its nearest virus's colour, matching the host figure
VIRUS_PALETTE = "--virus-palette" in sys.argv
VIRUS_PAL = ["#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4",
             "#F032E6", "#BFEF45", "#469990", "#9A6324", "#800000", "#000075",
             "#808000", "#FABED4", "#DCBEFF", "#AAFFC3"]
_vcount = Counter(r["nearest_reference_name"] for r in annot.values())
VIRUS_COL = {v: VIRUS_PAL[i % len(VIRUS_PAL)]
             for i, (v, _) in enumerate(_vcount.most_common())}

clade_col, clade_label, clade_base = {}, {}, {}
parts = cut_clades(min(NCLADE, len(CLADE_PAL)))
final = []
for i, part in enumerate(parts):
    base = CLADE_PAL[i % len(CLADE_PAL)]
    # --flat-clades gives one colour per clade, no sister sub-clade tones
    kids = None if "--flat-clades" in sys.argv else (
        part["children"] if (part["children"] and len(tips(part)) >= 4) else None)
    if kids and len(kids) >= 2:
        for j, k in enumerate(kids):
            final.append((k, tone(base, j), base))
    else:
        final.append((part, base, base))
if VIRUS_PALETTE:
    # Give each clade its commonest virus, repeats of a virus taking the next tone
    rank = Counter()
    recol = []
    for part, col, base in final:
        vs = Counter(annot[t["name"]]["nearest_reference_name"]
                     for t in tips(part) if t["name"] in annot)
        if vs:
            v = vs.most_common(1)[0][0]
            b = VIRUS_COL.get(v, "#808080")
            recol.append((part, tone(b, rank[v]), b)); rank[v] += 1
        else:
            recol.append((part, "#A9A9A9", "#A9A9A9"))
    final = recol

for i, (part, col, base) in enumerate(final):
    for t in tips(part):
        clade_of[t["name"]] = i
    clade_col[i] = col
    clade_base[i] = base
    clade_label[i] = ""

# Step 5: give every internal node the group its tips belong to, which colours a branch
GROUP = sys.argv[sys.argv.index("--by") + 1] if "--by" in sys.argv else "cluster"
NOBANDS = "--no-bands" in sys.argv
PLACE_PAL = {}
# --by clade follows the tree's own clades, clade+place keeps each tip's place colour
CLADE_MODE = GROUP in ("clade", "clade+place")
def group_key(name):
    r = annot.get(name)
    if not r: return None
    if CLADE_MODE:
        return clade_of.get(name)
    if GROUP == "place":
        return place_label(r)
    return int(r["cluster"])
def group_colour(k):
    if CLADE_MODE:
        return dark(clade_col.get(k, "#A9A9A9"))
    if GROUP != "place":
        return dark(clust_col[k])
    if k not in PLACE_PAL:
        for t in order:
            if t["name"] in annot and place_label(annot[t["name"]]) == k:
                PLACE_PAL[k] = dark(place_colour(annot[t["name"]])); break
        else:
            PLACE_PAL[k] = GREY
    return PLACE_PAL[k]
def group_name(k, n):
    if CLADE_MODE:
        return f"n={n}"
    return (f"{k}  (n={n})" if GROUP == "place"
            else f"CD-HIT cluster {k}  (n={n})")

def tip_cluster(t):
    return group_key(t["name"])
def tag(nd):
    if not nd["children"]:
        nd["cl"] = tip_cluster(nd); nd["ncontig"] = 1 if nd["cl"] is not None else 0; return
    cs, nc = set(), 0
    for k in nd["children"]:
        tag(k); nc += k["ncontig"]
        if k["ncontig"]: cs.add(k["cl"])
    nd["ncontig"] = nc; nd["cl"] = cs.pop() if len(cs) == 1 else None
tag(root)
def col_of(nd):
    return group_colour(nd["cl"]) if (nd["cl"] is not None and nd["ncontig"]) else GREY

# Step 6: lay the figure out, tips, labels, heatmap strips, clade blocks + legend
fig_h = max(4.0, 0.105 * n) if n > 45 else max(4.0, 0.28 * n)
if HEIGHT: fig_h = HEIGHT
fig = plt.figure(figsize=(13, fig_h + 2.0))
gs = fig.add_gridspec(2, 1, height_ratios=[fig_h, 2.0], hspace=0.02)
ax = fig.add_subplot(gs[0]); lax = fig.add_subplot(gs[1]); lax.axis("off")
colw = span * 0.030
gap = span * 0.030
# Place the names straight off the tree, then the markers + the matrix
tip_x = span + gap
labelw = span * (0.02 if SLIDE else (0.30 if n > 110 else 0.52))
x_origin = tip_x + labelw
x_host = x_origin + colw * 1.25
x_heat = [x_host + colw * (1.6 + 1.15 * i) for i in range(4)]
x_block_end = x_heat[-1] + colw
bar_x = x_block_end + span * 0.020
lab_x = bar_x + span * 0.012
right = lab_x + span * 0.24

# Draw the clade blocks behind everything, like the section bands
groups = []
for t in order:
    k = group_key(t["name"])
    if k is not None and k not in groups: groups.append(k)
groups = sorted(groups, key=lambda z: str(z))
blocks = []
for c in groups:
    ys = sorted(t["y"] for t in order if group_key(t["name"]) == c)
    runs, start, prev = [], ys[0], ys[0]
    for y in ys[1:]:
        if y - prev > 1.5: runs.append((start, prev)); start = y
        prev = y
    runs.append((start, prev))
    blocks.append((c, ys, runs))
placed = []
for c, ys, runs in sorted(blocks, key=lambda z: -len(z[1])):
    col = group_colour(c)
    band = vivid(col)
    for y0, y1 in (() if NOBANDS else runs):
        ax.add_patch(Rectangle((-span * 0.015, y0 - 0.5), x_block_end + span * 0.015,
                               y1 - y0 + 1.0, facecolor=band, alpha=0.28, edgecolor="none",
                               zorder=-5))
    # Draw one bracket per contiguous run, so a bar never spans tips of another place
    rows = [annot[t["name"]] for t in order if group_key(t["name"]) == c]
    for y0, y1 in runs:
        ax.plot([bar_x, bar_x], [y0 - 0.45, y1 + 0.45], color=col, lw=LW * 2.6,
                solid_capstyle="butt", zorder=4)
        tick = span * 0.004
        for yy in (y0 - 0.45, y1 + 0.45):
            ax.plot([bar_x, bar_x + tick], [yy, yy], color=col, lw=LW * 1.1,
                    solid_capstyle="butt", zorder=4)
    extra = ("" if GROUP == "place" or CLADE_MODE else
             "  " + ", ".join(f"{p} {k}" for p, k in
                              Counter(r["province"] or r["country"]
                                      for r in rows).most_common(2)))
    if CLADE_MODE:
        places = Counter(r["province"] or r["country"] for r in rows)
        extra = "  " + ", ".join(f"{p} {k}" for p, k in places.most_common(3))
    # Label every run of 3 or more, smaller ones sharing the biggest
    runs_sorted = sorted(runs, key=lambda r: -(r[1] - r[0]))
    first = True
    for y0, y1 in runs_sorted:
        n_here = sum(1 for t in order if group_key(t["name"]) == c and y0 <= t["y"] <= y1)
        # Skip the one-offs on a slide, as there is room only for the clades worth naming
        if SLIDE and n_here < MINLAB:
            continue
        mid = (y0 + y1) / 2
        step = n * (0.030 if SLIDE else 0.009)
        # Nudge off any label already at this height rather than drop it
        for _ in range(60):
            if all(abs(mid - q) >= step for q in placed): break
            mid += step
        placed.append(mid)
        if CLADE_MODE:
            here = [annot[t["name"]] for t in order
                    if group_key(t["name"]) == c and y0 <= t["y"] <= y1]
            pl = Counter(short_place(r) for r in here)
            txt = ", ".join(f"{q} {k}" for q, k in pl.most_common(3)) or f"n={n_here}"
            # Put the country once at the front, as province names read bare without it
            cn_provs = {(r.get("province") or "").strip() for r in here
                        if r.get("country") == "China" and (r.get("province") or "").strip()}
            if cn_provs:
                txt = "China \u2013 " + txt
        else:
            txt = (group_name(c, len(ys)) + extra if first
                   else f"{group_name(c, n_here).split('  (')[0]}  (n={n_here})")
        ax.text(lab_x, mid, txt, va="center", ha="left",
                fontsize=FS_CL if n_here >= 3 else FS_CL * 0.82, color=col,
                fontweight="bold" if first and len(ys) >= 5 else "normal")
        first = False

def draw(nd):
    if nd.get("collapsed"):
        return
    for k in nd["children"]:
        c = col_of(k); lw = LW * 1.15 if c != GREY else LW * 0.5
        ax.plot([nd["x"], k["x"]], [k["y"], k["y"]], color=c, lw=lw,
                solid_capstyle="round", zorder=2 if c != GREY else 1)
        if k.get("collapsed"):
            from matplotlib.patches import Polygon
            ax.add_patch(Polygon([(k["x"], k["y"]), (k["xmax"], k["y"] - 0.42),
                                  (k["xmax"], k["y"] + 0.42)], closed=True,
                                 facecolor="#D9D9D9", edgecolor="#9A9A9A", lw=0.3, zorder=1))
        draw(k)
    if nd["children"]:
        ys = [k["y"] for k in nd["children"]]
        c = col_of(nd)
        ax.plot([nd["x"], nd["x"]], [min(ys), max(ys)], color=c,
                lw=LW if c != GREY else LW * 0.45, zorder=2 if c != GREY else 1)
draw(root)

for t in order:
    r = annot.get(t["name"]); y = t["y"]
    if t.get("collapsed"):
        if not SLIDE:
            ax.text(tip_x, y, f"{t['nref']} references", va="center", fontsize=FS_TIP,
                    color="#7A7A7A", style="italic")
        continue
    if not SLIDE:
        ax.plot([t["x"], tip_x - span * 0.004], [y, y], color="#D5D5D5", lw=0.2,
                ls=(0, (1, 2)), zorder=0)
    if not r:
        if not SLIDE:
            acc = t["name"]
            nm = refname.get(acc)
            txt = f"{nm} ({acc})" if nm else acc
            if refcountry.get(acc):
                txt += f"  [{refcountry[acc]}]"
            ax.text(tip_x, y, txt, va="center", fontsize=FS_TIP, color="#777777",
                    style="italic")
        continue
    place = short_place(r)
    if not SLIDE:
        ax.text(tip_x, y, f"{t['name']}  {place}  (c{r['cluster']})", va="center",
                fontsize=FS_TIP,
                color=(dark(place_colour(r)) if GROUP == "clade+place"
                       else group_colour(group_key(t["name"]))))
    # Mark the origin, a triangle for a Chinese province, a circle for any other country
    mk = "^" if r["country"] == "China" else "o"
    fc = place_colour(r)
    ax.plot([x_origin], [y], mk, ms=MS, mfc=fc, mec="#404040", mew=0.25, zorder=5)
    hg = r.get("host_group") or "unknown"
    ax.plot([x_host], [y], "o", ms=MS, mfc=HOST.get(hg, "#FFFFFF"), mec="#404040",
            mew=0.25, zorder=5)
    for (title_, cmapname, fn), (lo, hi), x0 in zip(HEAT, heat_vals, x_heat):
        v = fn(r, t["name"])
        if v is None:
            ax.add_patch(Rectangle((x0, y - 0.45), colw, 0.9, facecolor="#FFFFFF",
                                   edgecolor="#DDDDDD", lw=0.2, hatch="///", zorder=3))
            continue
        cm = cmapname if hasattr(cmapname, "__call__") else plt.get_cmap(cmapname)
        col = cm(Normalize(lo, hi)(v))
        ax.add_patch(Rectangle((x0, y - 0.45), colw, 0.9, facecolor=col,
                               edgecolor="white", lw=0.15, zorder=3))


# --sisters marks one pair of sister clades, for explaining what a clade is
if "--sisters" in sys.argv:
    best = []
    def scan(nd):
        if len(nd["children"]) == 2:
            a, b = nd["children"]
            ca = {group_key(t["name"]) for t in tips(a)} - {None}
            cb = {group_key(t["name"]) for t in tips(b)} - {None}
            if ca and cb and not (ca & cb):
                best.append((min(len(ca), len(cb)), len(tips(a)) + len(tips(b)),
                             nd, a, b))
        for k in nd["children"]: scan(k)
    scan(root)
    if best:
        _, _, mrca, ka, kb = max(best, key=lambda z: (z[0], z[1]))
        ax.plot([mrca["x"]], [mrca["y"]], "o", ms=MS * 1.5, mfc="#FFFFFF",
                mec="#2B2B2B", mew=1.1, zorder=6)
        ax.annotate("their most recent\ncommon ancestor", xy=(mrca["x"], mrca["y"]),
                    xytext=(mrca["x"] - span * 0.22, mrca["y"] + n * 0.09),
                    fontsize=FS_CL * 0.95, color="#2B2B2B", ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.8), zorder=6)
        ys_a = [t["y"] for t in tips(ka)]; ys_b = [t["y"] for t in tips(kb)]
        bx = right + span * 0.03
        for ys_, nm in ((ys_a, "clade"), (ys_b, "sister clade")):
            ax.plot([bx, bx], [min(ys_), max(ys_)], color="#2B2B2B", lw=1.3)
            for yy in (min(ys_), max(ys_)):
                ax.plot([bx, bx - span * 0.012], [yy, yy], color="#2B2B2B", lw=1.3)
            ax.text(bx + span * 0.015, sum(ys_) / len(ys_), nm, va="center", ha="left",
                    fontsize=FS_CL, color="#2B2B2B", fontweight="bold")
        right = right + span * 0.20
ax.set_ylim(n + 1, -3)
ax.set_xlim(-span * 0.015, right)
ax.set_yticks([])
for side in ("top", "right", "left"): ax.spines[side].set_visible(False)
if "--axis" in sys.argv:
    ax.set_xlabel("substitutions per site", fontsize=7)
    ax.tick_params(axis="x", labelsize=6)
else:
    ax.spines["bottom"].set_visible(False)
    ax.set_xticks([])
    nice = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    step = min(nice, key=lambda v: abs(v - span / 6))
    y0 = n + 0.5
    ax.plot([0, step], [y0, y0], color="#333333", lw=1.6, solid_capstyle="butt",
            clip_on=False)
    ax.text(step / 2, y0 + 0.9, f"{step:g} substitutions per site", ha="center", va="top",
            fontsize=6.5, color="#333333")

for x0, h in zip([x_origin, x_host] + x_heat,
                 [TITLE_P, "mosquito host"] + [t_ for t_, _, _ in HEAT]):
    ax.text(x0 + (colw / 2 if x0 in x_heat else 0), -1.0, h.replace("\n", " "),
            rotation=90, ha="center", va="bottom", fontsize=FS_HDR)

shown = {place_label(annot[t["name"]]) for t in order if t["name"] in annot}
swatch = {}
for t in order:
    if t["name"] in annot:
        r = annot[t["name"]]
        swatch[place_label(r)] = (place_colour(r), "^" if r["country"] == "China" else "o")
cn_keys, other_keys = place_order()
# Draw the legends + scales in the panel underneath
L = [Line2D([], [], ls="", marker=swatch[k][1], ms=5, mfc=swatch[k][0], mec="#404040",
            mew=0.3, label=k) for k in cn_keys + other_keys if k in shown]
leg1 = lax.legend(handles=L, loc="upper left", bbox_to_anchor=(0.0, 0.80), fontsize=5.6,
                  title=TITLE_P, title_fontsize=6.5, frameon=False,
                  ncol=3 if len(L) > 10 else 1, handletextpad=0.4, columnspacing=1.0,
                  labelspacing=0.35)
h2 = [Line2D([], [], ls="", marker="o", ms=5, mfc=v, mec="#404040", mew=0.3, label=k)
      for k, v in HOST.items()]
h2.append(Line2D([], [], ls="-", lw=1.4, color=GREY, label="reference sequence / mixed clade"))
h2.append(Line2D([], [], ls="-", lw=1.6, color="#1F78B4",
                 label=("branch colour = the clade it belongs to" if CLADE_MODE else
                        "branch colour = every contig below it is from one place" if GROUP == "place"
                        else "branch colour = every contig below it is one CD-HIT cluster")))
h2.append(Line2D([], [], ls="", marker="s", ms=5, mfc="white", mec="#BBBBBB",
                 label="hatched cell = value not recorded"))
leg2 = lax.legend(handles=h2, loc="upper left", bbox_to_anchor=(0.40, 0.80), fontsize=5.6,
                  title="mosquito host / branches", title_fontsize=6.5, frameon=False,
                  handletextpad=0.4, labelspacing=0.35)
lax.add_artist(leg1)
for i, ((title_, cmapname, fn), (lo, hi)) in enumerate(zip(HEAT, heat_vals)):
    cax = lax.inset_axes([0.78, 0.72 - 0.20 * i, 0.16, 0.035])
    cb = fig.colorbar(ScalarMappable(Normalize(lo, hi), cmap=cmapname), cax=cax,
                      orientation="horizontal")
    cb.set_label(title_.replace("\n", " "), fontsize=5.6, labelpad=2)
    cb.ax.tick_params(labelsize=5.0, length=1.8, pad=1.2)
    cb.set_ticks([lo, hi])
    cb.ax.set_xticklabels([f"{lo:g}", f"{hi:g}"])
ax.set_title(TITLE, fontsize=9.5, loc="left")
# Step 7: write the virus -> colour table, so the host figure matches this one
if "--dump-colours" in sys.argv:
    path = sys.argv[sys.argv.index("--dump-colours") + 1]
    seen = {}
    for t in order:
        r = annot.get(t["name"])
        if not r: continue
        v = r["nearest_reference_name"]
        if VIRUS_PALETTE:
            seen.setdefault(v, Counter())[VIRUS_COL.get(v, "#808080")] += 1
        else:
            k = group_key(t["name"])
            seen.setdefault(v, Counter())[matplotlib.colors.to_hex(
                clade_base.get(k, "#808080"))] += 1
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["virus", "colour"])
        for v, c in sorted(seen.items()): w.writerow([v, c.most_common(1)[0][0]])
    print("wrote", path)

# Step 8: save a PDF for the report + a PNG for the slides
fig.savefig(OUT, bbox_inches="tight")
fig.savefig(OUT.replace(".pdf", ".png"), dpi=170, bbox_inches="tight")
mono = sum(1 for c, ys, runs in blocks if len(runs) == 1)
print("wrote", OUT, f"({n} tips, {len(clusters)} clusters, {mono} contiguous, "
                    f"{len(provs)} provinces, {len(others)} other countries)")
