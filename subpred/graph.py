from copy import deepcopy
from subpred.util import load_df
import obonet
import networkx as nx


def get_protein_dataset(
    organism_ids: set = set(), evidence_at_protein_level=True, reviewed: bool = True
):
    df_uniprot = load_df("uniprot", "data/datasets")
    if organism_ids:
        df_uniprot = df_uniprot[df_uniprot.organism_id.isin(organism_ids)]
    if evidence_at_protein_level:
        df_uniprot = df_uniprot[df_uniprot.protein_existence == 1]
    if reviewed:
        df_uniprot = df_uniprot[df_uniprot.reviewed]
    return df_uniprot


def get_go_annotations(proteins: set, include_iea: bool = True):
    df_goa_uniprot = load_df("go", "data/datasets")
    df_goa_uniprot = df_goa_uniprot[df_goa_uniprot.Uniprot.isin(proteins)]
    if not include_iea:
        df_goa_uniprot = df_goa_uniprot[df_goa_uniprot.evidence_code != "IEA"]
    df_goa_uniprot = df_goa_uniprot.reset_index(drop=True)
    return df_goa_uniprot


def get_graph_go():
    graph_go = obonet.read_obo("data/raw/ontologies/go.obo", ignore_obsolete=True)
    # TODO filter?
    return graph_go


def add_ancestors(df_uniprot_goa, graph_go):
    df_uniprot_goa = df_uniprot_goa.assign(
        ancestors=[
            nx.descendants(graph_go, go_id) | {go_id} for go_id in df_uniprot_goa.go_id
        ]
    )

    df_uniprot_goa = (
        df_uniprot_goa.explode("ancestors")
        .drop("go_id", axis=1)
        .rename(columns={"ancestors": "go_id"})
        .drop_duplicates()
        .reset_index(drop=True)
    )

    return df_uniprot_goa


def get_substrate_graph(organism_ids: set):
    df_uniprot = get_protein_dataset(
        organism_ids=organism_ids, evidence_at_protein_level=True, reviewed=True
    )
    df_uniprot_goa = get_go_annotations(proteins=set(df_uniprot.index.tolist()))

    graph_go = get_graph_go()
    go_id_to_name = {id: data["name"] for id, data in graph_go.nodes(data=True)}
    go_name_to_id = {name: id for id, name in go_id_to_name.items()}

    graph_go_mf = graph_go.subgraph(
        nodes=[
            node
            for node, data in graph_go.nodes(data=True)
            if data["namespace"] == "molecular_function"
        ]
    )
    graph_go_mf = graph_go_mf.edge_subgraph(
        edges={edge for edge in graph_go_mf.edges(keys=True) if edge[2] == "is_a"}
    )
    # converting the immutable view returned by subgraph into a new graph
    graph_go_mf = deepcopy(graph_go_mf.copy())

    df_uniprot_goa_mf = df_uniprot_goa[df_uniprot_goa.aspect == "F"]

    # keeping go terms in graph that are not in goa, because they could be parts of paths
    # removing go annotations that are not in graph
    df_uniprot_goa_mf = df_uniprot_goa_mf[
        df_uniprot_goa_mf.go_id.isin(set(graph_go_mf.nodes()))
    ]
    df_uniprot_goa_mf = df_uniprot_goa_mf.reset_index(drop=True)

    df_uniprot_goa_mf_ancestors = add_ancestors(
        df_uniprot_goa=df_uniprot_goa_mf.copy(deep=True), graph_go=graph_go_mf
    )

    # df_goa_uniprot_ecoli[~df_goa_uniprot_ecoli.go_id.isin(set(graph_go.nodes()))].shape[0]

    # TODO only keep goa and graph GO terms that are in intersection?

    # TODO well-defined goal

    pass


import os

print(os.getcwd())

get_substrate_graph({83333})
