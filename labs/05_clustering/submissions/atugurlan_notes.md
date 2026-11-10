## Most appropriate method

K-means (K=2) was the most useful of the tested settings for summarizing the [Wisconsin Diagnostic Breast Cancer dataset](https://archive.ics.uci.edu/dataset/17/breast-cancer-wisconsin-diagnostic) into two broad groups. It assigned all 569 samples to clusters of 375 and 194 without using the diagnosis labels. 

With `eps=1.5` and `min_samples=5`, DBSCAN found one cluster of 19 samples and labeled 550 as noise. 

The hierarchical dendrogram is useful for inspecting how groups merge at different distances. The clusters are based on standardized features, while PCA is used only to display them.

## How does clustering compare to phylogenetic trees?

Clustering groups samples with similar measured features. A phylogenetic tree shows how organisms or sequences are thought to be related through evolution. In a rooted tree, the branching points represent inferred common ancestors. A clustering dendrogram may look similar, but its branches show how groups merge based on feature similarity, not ancestry.
