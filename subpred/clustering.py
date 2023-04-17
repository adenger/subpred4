# Functions related to unsupervised learning and clustering analysis

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA
import seaborn as sns
import pandas as pd
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
)
from scipy.cluster import hierarchy
from sklearn.preprocessing import scale


def get_linkage(
    feature_data,
    standardize: bool = True,
    method: str = "ward",
    metric: str = "euclidean",
):
    if standardize:
        feature_data = scale(feature_data, copy=True)
    linkage = hierarchy.linkage(feature_data, method=method, metric=metric)
    return linkage


def dendrogram_plot(linkage, max_leaves=15):
    return hierarchy.dendrogram(linkage, truncate_mode="lastp", p=max_leaves)


def get_cluster_labels(linkage, n_clusters=2, index: list = None):
    labels = hierarchy.fcluster(linkage, criterion="maxclust", t=n_clusters)
    labels = pd.Series(data=labels, index=index)
    return labels


def clustering_quality_plots(
    feature_data, min_clusters: int = 2, max_clusters: int = 20
):
    """
    Feature_data is a pandas dataframe with features as columns and samples as rows.
    - Elbow plot/Inertia: the more it looks like an elbow, the better. the optimal number of clusters is at the "joint of the arm"
    - Silhouette: Number between -1 and 1. The higher the better. Measures tightness/overlap of clusters
    - CH: Higher score is better
    - DB: The closer to 0 the better. Average similarity between any cluster and its closest cluster.
    """
    scores = list()

    for n_clusters in range(min_clusters, max_clusters + 1):
        pipe = make_pipeline(
            StandardScaler(), KMeans(n_clusters=n_clusters, random_state=0)
        )
        pipe.fit(feature_data)
        scores.append([n_clusters, "Inertia", pipe["kmeans"].inertia_])
        scores.append(
            [
                n_clusters,
                "Silhouette Coefficient",
                silhouette_score(feature_data, pipe["kmeans"].labels_),
            ]
        )
        scores.append(
            [
                n_clusters,
                "Calinski Harabasz Index",
                calinski_harabasz_score(feature_data, pipe["kmeans"].labels_),
            ]
        )
        scores.append(
            [
                n_clusters,
                "Davies Bouldin Index",
                davies_bouldin_score(feature_data, pipe["kmeans"].labels_),
            ]
        )

    scores_wide = pd.DataFrame.from_records(
        scores, columns=["k", "metric", "score"]
    ).pivot(index="k", columns="metric", values="score")
    return scores_wide.plot(
        subplots=True, layout=(2, 2), figsize=(15, 10), xticks=scores_wide.index
    )
