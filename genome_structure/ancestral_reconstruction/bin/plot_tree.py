#!/usr/bin/env python3
"""
Plot a Newick species tree with tip and internal node labels.

Usage:
    python3 plot_tree.py <tree.nwk> [output.pdf]

Requires only matplotlib (no ete3/biopython needed).
"""

import sys
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict


class Node:
    """Simple tree node."""
    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.x = 0.0  # horizontal (distance from root)
        self.y = 0.0  # vertical (leaf ordering)


def parse_newick(s):
    """Parse a Newick string into a Node tree."""
    s = s.strip().rstrip(';')
    pos = [0]

    def parse_node():
        node = Node()
        if s[pos[0]] == '(':
            pos[0] += 1  # skip '('
            while True:
                child = parse_node()
                node.children.append(child)
                if s[pos[0]] == ',':
                    pos[0] += 1
                elif s[pos[0]] == ')':
                    pos[0] += 1
                    break
        # Read name
        name = []
        while pos[0] < len(s) and s[pos[0]] not in (',', ')', ':', ';'):
            name.append(s[pos[0]])
            pos[0] += 1
        node.name = ''.join(name)
        # Read branch length
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            bl = []
            while pos[0] < len(s) and s[pos[0]] not in (',', ')', ';'):
                bl.append(s[pos[0]])
                pos[0] += 1
            node.length = float(''.join(bl))
        return node

    return parse_node()


def assign_positions(node, x=0.0, leaf_counter=None):
    """Assign x (cumulative branch length) and y (leaf order) positions."""
    if leaf_counter is None:
        leaf_counter = [0]

    node.x = x
    if not node.children:
        # Leaf
        node.y = leaf_counter[0]
        leaf_counter[0] += 1
    else:
        for child in node.children:
            assign_positions(child, x + child.length, leaf_counter)
        # Internal node y = midpoint of children
        node.y = (node.children[0].y + node.children[-1].y) / 2.0


def count_leaves(node):
    """Count total leaves in tree."""
    if not node.children:
        return 1
    return sum(count_leaves(c) for c in node.children)


def draw_tree(node, ax):
    """Recursively draw the tree."""
    for child in node.children:
        # Horizontal line from parent to child
        ax.plot([node.x, child.x], [child.y, child.y], 'k-', lw=1.2)
        # Vertical line connecting children
        ax.plot([node.x, node.x], [node.children[0].y, node.children[-1].y],
                'k-', lw=1.2)
        draw_tree(child, ax)


def collect_nodes(node, tips, internals):
    """Collect tip and internal nodes."""
    if not node.children:
        tips.append(node)
    else:
        internals.append(node)
        for child in node.children:
            collect_nodes(child, tips, internals)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_tree.py <tree.nwk> [output.pdf]")
        sys.exit(1)

    tree_file = sys.argv[1]
    output = sys.argv[2] if len(sys.argv) > 2 else tree_file.rsplit('.', 1)[0] + '_tree.pdf'

    with open(tree_file) as f:
        newick_str = f.read().strip()

    root = parse_newick(newick_str)
    assign_positions(root)

    n_leaves = count_leaves(root)
    tips, internals = [], []
    collect_nodes(root, tips, internals)

    # Figure size scales with number of tips
    fig_height = max(4, n_leaves * 0.5)
    fig_width = max(8, fig_height * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    draw_tree(root, ax)

    # Tip labels (right-aligned to tips)
    max_x = max(t.x for t in tips) if tips else 0
    for tip in tips:
        ax.text(tip.x + max_x * 0.01, tip.y, f"  {tip.name}",
                va='center', ha='left', fontsize=11, fontstyle='italic')

    # Internal node labels (circles + text)
    for nd in internals:
        if nd.name:
            ax.plot(nd.x, nd.y, 'o', color='steelblue', markersize=6, zorder=5)
            ax.text(nd.x, nd.y + 0.15, nd.name,
                    va='bottom', ha='center', fontsize=8,
                    fontweight='bold', color='steelblue')

    # Root label
    if root.name:
        ax.plot(root.x, root.y, 'o', color='steelblue', markersize=6, zorder=5)
        ax.text(root.x, root.y + 0.15, root.name,
                va='bottom', ha='center', fontsize=8,
                fontweight='bold', color='steelblue')

    ax.set_xlabel('Substitutions per site', fontsize=12)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_title('Species Tree with Node Labels', fontsize=14, pad=15)

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Tree saved to {output}")


if __name__ == "__main__":
    main()
