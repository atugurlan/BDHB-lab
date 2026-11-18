"""
Exercise — Gene Co-Expression Networks (GCEs) — Building the network and detecting modules

Objective:
- Build a co-expression network from an RNA-Seq expression matrix
- Detect gene modules (communities) using a Louvain-like algorithm (or alternative)

Instructions (in the lab):
1) Data preparation
   - Download and prepare the expression matrix (e.g., GSE115469) into a CSV with:
     * rows = genes (index), columns = samples (sample IDs)
   - Save the file to: data/work/<handle>/lab06/expression_matrix.csv

2) Preprocessing
   - log2(x + 1)
   - filter low-variance genes

3) Correlation → Adjacency
   - complete the `correlation_matrix` function
   - the function `adjacency_from_correlation` is already implemented

4) Graph + Modules
   - build the graph with NetworkX
   - detect modules (Louvain or alternative)
   - export the gene → module mapping to submissions/<handle>/modules_<handle>.csv

Note:
- Document in <github_handle>_notes.md: correlation metric, threshold, short observations.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import networkx as nx


# --------------------------
# Config — fill in as needed
# --------------------------
HANDLE = "atugurlan"
INPUT_CSV = Path(f"data/work/{HANDLE}/lab06/expression_matrix.csv")
OUTPUT_DIR = Path(f"labs/06_wgcna/submissions/{HANDLE}")
OUTPUT_CSV = OUTPUT_DIR / f"modules_{HANDLE}.csv"

CORR_METHOD = "spearman"
VARIANCE_THRESHOLD = 0.1
ADJ_THRESHOLD = 0.6
USE_ABS_CORR = True
MAKE_UNDIRECTED = True


def read_expression_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def log_and_filter(df: pd.DataFrame,
                   variance_threshold: float) -> pd.DataFrame:
    """
    Preprocessing:
    - apply log2(x+1)
    - filter genes with low variance
    """
    logged = np.log2(df + 1)
    variance = logged.var(axis=1)
    return logged.loc[variance > variance_threshold]


def correlation_matrix(df: pd.DataFrame,
                       method: str = "spearman",
                       use_abs: bool = True) -> pd.DataFrame:
    """
    TODO: compute the gene–gene correlation matrix (rows).
    """
    corr = df.T.corr(method=method)
    if use_abs:
        corr = corr.abs()
    return corr


def adjacency_from_correlation(corr: pd.DataFrame,
                               threshold: float,
                               weighted: bool = False) -> pd.DataFrame:
    """
    Build adjacency matrix from correlations.
    - binary: A_ij = 1 if corr_ij >= threshold, else 0
    - weighted: A_ij = corr_ij if corr_ij >= threshold, else 0
    """
    if weighted:
        adjacency = corr.where(corr >= threshold, 0)
    else:
        adjacency = (corr >= threshold).astype(int)

    for gene in adjacency.columns:
        adjacency.loc[gene, gene] = 0
    return adjacency


def graph_from_adjacency(A: pd.DataFrame,
                         undirected: bool = True) -> nx.Graph:
    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(A, create_using=nx.DiGraph)
    isolates = list(nx.isolates(G))
    if isolates:
        G.remove_nodes_from(isolates)
    return G


def detect_modules_louvain_or_greedy(G: nx.Graph) -> Dict[str, int]:
    """
    TODO: detect communities (modules) and return a dict gene -> module_id.
    Options:
      - try louvain_communities(G, seed=42) if available
      - otherwise greedy_modularity_communities(G)
    """
    if G.number_of_edges() == 0:
        return {}

    if hasattr(nx.community, "louvain_communities"):
        communities = nx.community.louvain_communities(G, seed=42)
    else:
        communities = nx.community.greedy_modularity_communities(G)

    mapping = {}
    for module_id, genes in enumerate(communities, start=1):
        for gene in genes:
            mapping[gene] = module_id
    return mapping


def save_modules_csv(mapping: Dict[str, int], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df_modules = (
        pd.DataFrame({"Gene": list(mapping.keys()), "Module": list(mapping.values())})
        .sort_values(["Module", "Gene"])
    )
    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":
    expression = read_expression_matrix(INPUT_CSV)
    filtered = log_and_filter(expression, VARIANCE_THRESHOLD)
    corr = correlation_matrix(filtered, CORR_METHOD, USE_ABS_CORR)
    adjacency = adjacency_from_correlation(corr, ADJ_THRESHOLD)
    G = graph_from_adjacency(adjacency, MAKE_UNDIRECTED)
    print(f"Graph created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

    gene_to_module = detect_modules_louvain_or_greedy(G)
    print(f"Detected {len(set(gene_to_module.values()))} modules.")

    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Saved gene→module mapping to: {OUTPUT_CSV}")
