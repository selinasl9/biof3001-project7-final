#!/usr/bin/env python3
"""
newick.py

Read an IQ-TREE .treefile and measure distances in it, with a node as a dict of
children, name, length, support and parent.

Used by: tables/build_contig_annotation.py, figures/plot_host_vs_virus.py,
figures/render_virus_tree.py, tables/make_result_tables.py

render_virus_tree.py keeps its own re-rooting, as it measures its clades from the cut
point it re-roots at.
"""

import re

# Set the patterns for a node label + a branch length
LABEL = re.compile(r"[^(),:;]*")
LENGTH = re.compile(r":(-?[0-9.eE+-]+)")


def parse(path):
    """Return the root node of the tree in `path`."""
    # Step 1: read the file, dropping the trailing ";" as it would stop the parser
    text = open(path).read().strip().rstrip(";")

    # Keep position in a list so the nested function can move it
    position = [0]

    def node():
        children = []
        if text[position[0]] == "(":
            position[0] += 1
            while True:
                children.append(node())
                if text[position[0]] == ",":
                    position[0] += 1
                else:
                    break
            position[0] += 1

        label = LABEL.match(text[position[0]:]).group(0)
        position[0] += len(label)

        length = 0.0
        if position[0] < len(text) and text[position[0]] == ":":
            match = LENGTH.match(text[position[0]:])
            length = float(match.group(1))
            position[0] += len(match.group(0))

        # Take an internal node's label as support, where IQ-TREE writes the bootstrap value
        return {"children": children,
                "name": "" if children else label,
                "support": label if children else "",
                "length": length,
                "parent": None}

    # Step 2: parse from the outermost bracket down, then fill the parents in
    root = node()
    link_parents(root)
    return root


def link_parents(node, parent=None):
    """Fill in each node's parent, so a tip can be walked back to the root."""
    node["parent"] = parent
    for child in node["children"]:
        link_parents(child, node)


def tips(node):
    """Return every tip at or below this node, left to right."""
    if not node["children"]:
        return [node]
    return [tip for child in node["children"] for tip in tips(child)]


def tip_dict(root):
    """Return {name: tip}, with MAFFT's _R_ prefix stripped so lookups match the contig ids."""
    found = {}
    for tip in tips(root):
        found[tip["name"]] = tip
        if tip["name"].startswith("_R_"):
            found[tip["name"][3:]] = tip
    return found


def collapse_single_children(node):
    """Drop nodes with one child, adding their branch length to that child."""
    for child in node["children"]:
        collapse_single_children(child)
    for index, child in enumerate(node["children"]):
        if len(child["children"]) == 1:
            only = child["children"][0]
            only["length"] += child["length"]
            only["parent"] = node
            node["children"][index] = only


def path_to_root(node):
    """Return the node, its parent, and so on up to the root."""
    walk, current = [], node
    while current is not None:
        walk.append(current)
        current = current["parent"]
    return walk


def patristic_distance(a, b):
    """Return the branch length between two nodes, through their common ancestor."""
    up_a, up_b = path_to_root(a), path_to_root(b)
    shared = set(id(n) for n in up_b)

    distance = 0.0
    for node in up_a:
        if id(node) in shared:
            ancestor = node
            break
        distance += node["length"]
    else:
        raise ValueError("the two nodes are not in the same tree")

    for node in up_b:
        if node is ancestor:
            break
        distance += node["length"]
    return distance


def nearest(tip, candidates):
    """Return the candidate tip closest to `tip`, and that distance."""
    best = min(candidates, key=lambda other: patristic_distance(tip, other))
    return best, patristic_distance(tip, best)


def ladderise(node):
    """Order each node's children by how many tips they hold, smallest first."""
    for child in node["children"]:
        ladderise(child)
    node["children"].sort(key=lambda child: len(tips(child)))


def reroot_above(tip):
    """Put a new root on the branch above `tip`, splitting that branch in half."""
    walk = path_to_root(tip)
    lengths = [node["length"] for node in walk]

    # Flip every edge between the tip and the old root
    for child, parent in zip(walk, walk[1:]):
        parent["children"].remove(child)
        child["children"].append(parent)
    for index in range(1, len(walk)):
        walk[index]["length"] = lengths[index - 1]

    old_root_side = walk[1]
    tip["children"].remove(old_root_side)
    half = lengths[0] / 2
    tip["length"] = old_root_side["length"] = half

    root = {"children": [tip, old_root_side], "name": "", "support": "",
            "length": 0.0, "parent": None}
    link_parents(root)
    return root


def midpoint_root(root):
    """Re-root halfway along the longest tip-to-tip path."""
    # Step 1: find the two tips furthest apart
    all_tips = tips(root)
    far_a, far_b, longest = all_tips[0], all_tips[0], -1.0
    for i, a in enumerate(all_tips):
        for b in all_tips[i + 1:]:
            span = patristic_distance(a, b)
            if span > longest:
                far_a, far_b, longest = a, b, span

    # Step 2: walk up from one of them until half the longest path is behind
    walk = path_to_root(far_a)
    shared = set(id(n) for n in path_to_root(far_b))
    travelled = 0.0
    for node in walk:
        if id(node) in shared:
            break
        if travelled + node["length"] >= longest / 2:
            above = node["length"] - (longest / 2 - travelled)
            return _split_branch(root, node, longest / 2 - travelled, above)
        travelled += node["length"]
    return root


def _split_branch(root, node, below, above):
    """Put a new root on `node`'s branch, `below` away from node and `above` from its parent."""
    parent = node["parent"]
    parent["children"].remove(node)
    node["length"] = below

    # Turn the old root's half of the tree around, as it is on the other side
    walk = path_to_root(parent)
    lengths = [n["length"] for n in walk]
    for child, grandparent in zip(walk, walk[1:]):
        grandparent["children"].remove(child)
        child["children"].append(grandparent)
    for index in range(1, len(walk)):
        walk[index]["length"] = lengths[index - 1]
    parent["length"] = above

    new_root = {"children": [node, parent], "name": "", "support": "",
                "length": 0.0, "parent": None}
    link_parents(new_root)
    collapse_single_children(new_root)
    return new_root
