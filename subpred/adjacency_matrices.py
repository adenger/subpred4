import networkx as nx
import pandas as pd
from subpred.util import load_df

def get_adjacency_matrix(graph, labels: list, edges_filter: set = {"is_a"}):
    subgraph = graph.edge_subgraph(
        [edge for edge in graph.edges(keys=True) if edge[2] in edges_filter]
    )
    assert len(labels) == len(set(labels)), "labels should only contain unique elements"
    subgraph = subgraph.subgraph(labels)

    # scipy sparse matrix
    df_adjacency_matrix = nx.adjacency_matrix(G=subgraph, nodelist=labels)

    df_adjacency_matrix = pd.DataFrame(
        df_adjacency_matrix.todense(), columns=labels, index=labels
    )
    return df_adjacency_matrix


def get_go_adjacency_matrix(
    df_uniprot_goa: pd.DataFrame,
    datasets_location: str = "../data/datasets/",
    edges_filter: set = {"is_a"},
):
    graph_go = load_df("go_obo", folder_path=datasets_location)
    go_ids = sorted(df_uniprot_goa.go_id_ancestor.unique())
    df_adj_matrix_go = get_adjacency_matrix(
        graph_go, labels=go_ids, edges_filter=edges_filter
    )
    return df_adj_matrix_go


def get_chebi_adjacency_matrix(
    df_go_chebi: pd.DataFrame,
    datasets_location: str = "../data/datasets/",
    edges_filter: set = {"is_a"},
    primary_substrate_only: bool = True,
):
    graph_chebi = load_df("chebi_obo", folder_path=datasets_location)
    df_go_chebi_copy = (
        df_go_chebi[df_go_chebi.chebi_go_relation == "has_primary_input"].copy()
        if primary_substrate_only
        else df_go_chebi.copy()
    )
    chebi_id_primary = sorted(
        df_go_chebi_copy.chebi_id.unique()
    )
    df_adj_matrix_chebi = get_adjacency_matrix(
        graph_chebi.copy(), labels=chebi_id_primary, edges_filter=edges_filter
    )
    return df_adj_matrix_chebi