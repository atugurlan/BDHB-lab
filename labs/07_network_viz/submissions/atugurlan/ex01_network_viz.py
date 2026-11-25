"""
Exercise 1 — Visualization of Co-Expression Networks + Hub Genes

TODO:
- Load the expression matrix and module mapping from Lab 6
- Rebuild (or load) the adjacency matrix
- Construct the graph from adjacency
- Color nodes by module
- Compute hub genes (top degree)
- Visualize and export the network figure (.png)
- Export hub genes to CSV (optional)
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# --------------------------
# Config — complete with your values
# --------------------------
HANDLE = "atugurlan"

# Input files
EXPR_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
MODULES_CSV = Path(f"labs/06_wgcna/submissions/{HANDLE}/modules_{HANDLE}.csv")

# Optional: if you saved adjacency in Lab 6, load it here
PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Parameters for adjacency reconstruction (if PRECOMPUTED_ADJ_CSV is None)
CORR_METHOD = "spearman"   # TODO: choose "pearson" or "spearman"
VARIANCE_THRESHOLD = 0.1
USE_ABS_CORR = True        # TODO: use absolute correlations?
ADJ_THRESHOLD = 0.6        # TODO: correlation threshold
WEIGHTED = False           # TODO: use weighted adjacency or binary?

# Visualization parameters
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 60
EDGE_ALPHA = 0.15

# Outputs
OUT_DIR = Path(f"labs/07_network_viz/submissions/{HANDLE}")
OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# Utils
# --------------------------
def ensure_exists(path: Path) -> None:
    """Check that a file exists."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")


def read_expression_matrix(path: Path) -> pd.DataFrame:
    """
    TODO:
    - read CSV
    - set index to gene names
    """
    return pd.read_csv(path, index_col=0)


def read_modules_csv(path: Path) -> Dict[str, int]:
    """
    TODO:
    - read CSV with columns: Gene, Module
    - return dict: gene -> module_id
    """
    modules = pd.read_csv(path)
    return dict(zip(modules["Gene"], modules["Module"]))


def correlation_to_adjacency(expr: pd.DataFrame,
                             method: str,
                             use_abs: bool,
                             threshold: float,
                             weighted: bool) -> pd.DataFrame:
    """
    TODO:
    - compute correlation matrix on expr
    - optionally apply abs()
    - apply threshold to build adjacency
    - remove diagonal
    """
    logged = np.log2(expr + 1)
    variance = logged.var(axis=1)
    filtered = logged.loc[variance > VARIANCE_THRESHOLD]
    corr = filtered.T.corr(method=method)
    if use_abs:
        corr = corr.abs()

    if weighted:
        adjacency = corr.where(corr >= threshold, 0)
    else:
        adjacency = (corr >= threshold).astype(int)
    for gene in adjacency.columns:
        adjacency.loc[gene, gene] = 0
    return adjacency


def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    """
    TODO:
    - convert adjacency DataFrame to NetworkX graph
    - remove isolated nodes
    """
    G = nx.from_pandas_adjacency(A)
    G.remove_nodes_from(list(nx.isolates(G)))
    return G


def color_map_from_modules(nodes: Iterable[str], gene2module: Dict[str, int]) -> Dict[str, tuple]:
    """
    TODO:
    - assign a color to each node based on its module
    - use matplotlib 'tab10' or another palette
    """
    cmap = plt.get_cmap("tab10")
    colors = {}
    for gene in nodes:
        colors[gene] = cmap((gene2module[gene] - 1) % 10)
    return colors


def compute_hubs(G: nx.Graph, topk: int) -> pd.DataFrame:
    """
    TODO:
    - compute degree for every node
    - (optional) compute betweenness centrality
    - return top-k genes
    """
    degrees = dict(G.degree())
    hubs = sorted(degrees.items(), key=lambda item: (-item[1], item[0]))[:topk]
    return pd.DataFrame(hubs, columns=["Gene", "Degree"])


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # 1: Verify input files exist
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 2: Load expression matrix and module mapping
    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    # 3: Load or reconstruct adjacency
    if PRECOMPUTED_ADJ_CSV is not None:
        ensure_exists(PRECOMPUTED_ADJ_CSV)
        A = pd.read_csv(PRECOMPUTED_ADJ_CSV, index_col=0)
    else:
        A = correlation_to_adjacency(
            expr, CORR_METHOD, USE_ABS_CORR, ADJ_THRESHOLD, WEIGHTED)

    # 4: Build graph
    G = graph_from_adjacency(A)
    if set(G.nodes()) != set(gene2module):
        raise ValueError("The reconstructed graph does not match the Lab 6 module genes.")
    print(f"Graph created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

    # 5: Compute colors by module
    colors = color_map_from_modules(G.nodes(), gene2module)
    node_colors = []
    for gene in G.nodes():
        node_colors.append(colors[gene])

    # 6: Compute hub genes
    hubs_df = compute_hubs(G, TOPK_HUBS)
    hub_nodes = set(hubs_df["Gene"])
    node_sizes = []
    for gene in G.nodes():
        if gene in hub_nodes:
            node_sizes.append(NODE_BASE_SIZE * 3)
        else:
            node_sizes.append(NODE_BASE_SIZE)
    print(hubs_df.to_string(index=False))

    # 7: Compute layout and draw graph
    pos = nx.spring_layout(G, seed=SEED)
    plt.figure(figsize=(13, 10))
    nx.draw_networkx_edges(G, pos, alpha=EDGE_ALPHA, width=0.5)
    nx.draw_networkx_nodes(
        G, pos, node_color=node_colors, node_size=node_sizes,
        edgecolors="black", linewidths=0.3)
    ax = plt.gca()
    label_y = 0.9
    for gene in hubs_df["Gene"]:
        ax.annotate(
            gene, xy=pos[gene], xytext=(1.02, label_y),
            textcoords="axes fraction", fontsize=9, annotation_clip=False,
            arrowprops=dict(arrowstyle="-", color="gray", linewidth=0.6))
        label_y -= 0.055
    plt.title("Gene co-expression network - spring layout")
    plt.axis("off")

    cmap = plt.get_cmap("tab10")
    for module in sorted(set(gene2module.values())):
        plt.scatter([], [], c=[cmap((module - 1) % 10)], label=f"Module {module}")
    plt.legend(title="Node colors", loc="upper left", fontsize=8)

    # 8: Save network figure
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=200)
    plt.close()

    # 9: Save hub genes to CSV
    hubs_df.to_csv(OUT_HUBS, index=False)
