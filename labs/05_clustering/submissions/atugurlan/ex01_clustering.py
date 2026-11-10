"""
Exercise 6 — Clustering on breast cancer data (toy dataset)

Instructions:
1. Load the WDBC breast cancer dataset from the UCI Repository.
2. Preprocess the data: remove irrelevant columns and convert diagnosis to numeric values.
3. Standardize the data.
4. Implement and visualize clustering using:
   - Hierarchical clustering (dendrogram),
   - K-means (K=2, PCA visualization),
   - DBSCAN (PCA visualization).
5. Save results in submissions/<handle>/:
   - clusters_<handle>.csv
   - hierarchical_<handle>.png
   - kmeans_<handle>.png
   - dbscan_<handle>.png
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# Optional: you may import necessary functions now
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA

if __name__ == "__main__":
    # TODO 1: Load dataset
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
    columns = ["ID", "Diagnosis"] + [f"Feature_{i}" for i in range(1, 31)]
    df = pd.read_csv(url, header=None, names=columns)

    # TODO 2: Preprocessing
    # - remove the ID column
    df = df.drop(columns=["ID"])
    df["Diagnosis"] = df["Diagnosis"].map({"B": 0, "M": 1})

    # TODO 3: Standardization
    X = df.drop(columns=["Diagnosis"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Directory for results
    handle = "atugurlan"
    output_dir = Path(f"labs/05_clustering/submissions/{handle}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # TODO 4: Hierarchical Clustering
    # - use linkage(X_scaled, method="average")
    # - visualize with dendrogram()
    # - save figure as hierarchical_<handle>.png
    Z = linkage(X_scaled, method="average")
    plt.figure(figsize=(10, 6))
    dendrogram(Z, no_labels=True)
    plt.title("Hierarchical clustering (average linkage)")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig(output_dir / f"hierarchical_{handle}.png", dpi=160)
    plt.close()

    # TODO 5: K-means Clustering
    # - apply KMeans with K=2
    # - add labels in df["KMeans_Cluster"]
    # - reduce dimensionality with PCA(n_components=2)
    # - visualize and save plot kmeans_<handle>.png
    kmeans = KMeans(n_clusters=2, random_state=0, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_scaled)
    df["KMeans_Cluster"] = kmeans_labels

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=kmeans_labels, cmap="viridis", s=30, alpha=0.7)
    plt.title("K-means clustering (K=2)")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.colorbar(ticks=[0, 1], label="Cluster")
    plt.tight_layout()
    plt.savefig(output_dir / f"kmeans_{handle}.png", dpi=160)
    plt.close()

    # TODO 6: DBSCAN Clustering
    # - apply DBSCAN (e.g., eps=1.5, min_samples=5)
    # - add labels in df["DBSCAN_Cluster"]
    # - visualize and save plot dbscan_<handle>.png
    dbscan = DBSCAN(eps=1.5, min_samples=5)
    dbscan_labels = dbscan.fit_predict(X_scaled)
    df["DBSCAN_Cluster"] = dbscan_labels

    plt.figure(figsize=(8, 6))
    for label in sorted(set(dbscan_labels)):
        points = X_pca[dbscan_labels == label]
        if label == -1:
            plt.scatter(points[:, 0], points[:, 1], c="gray", s=30, alpha=0.6, label="Noise (-1)")
        else:
            plt.scatter(points[:, 0], points[:, 1], s=30, alpha=0.7, label=f"Cluster {label}")
    plt.title("DBSCAN clustering (eps=1.5, min_samples=5)")
    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / f"dbscan_{handle}.png", dpi=160)
    plt.close()

    # TODO 7: Save results
    # save a CSV with columns ["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]
    # to clusters_<handle>.csv
    df[["Diagnosis", "KMeans_Cluster", "DBSCAN_Cluster"]].to_csv(
        output_dir / f"clusters_{handle}.csv", index=False)
    print("K-means labels:")
    print(df["KMeans_Cluster"].value_counts().sort_index())
    print("DBSCAN labels (-1 means noise):")
    print(df["DBSCAN_Cluster"].value_counts().sort_index())
    print("Saved results in:", output_dir)
