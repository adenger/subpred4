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


def get_substrate_graph(organism_ids: set, go_obo_path: str, chebi_obo_path: str):
    # Read go graph
    graph_go = obonet.read_obo(go_obo_path, ignore_obsolete=True)
    # Filter graph by organism
    df_uniprot = get_protein_dataset(
        organism_ids=organism_ids, evidence_at_protein_level=True, reviewed=True
    )
    df_uniprot_goa = get_go_annotations(
        proteins=set(df_uniprot.index.tolist()), include_iea=True
    )
    df_uniprot_goa_ancestors = add_ancestors(
        df_uniprot_goa=df_uniprot_goa, graph_go=graph_go
    )
    print(len(graph_go.nodes()))
    graph_go = graph_go.subgraph(nodes=set(df_uniprot_goa_ancestors.go_id.unique()))
    print(len(graph_go.nodes()))

    # Filter graph by aspect/namespace
    graph_go = graph_go.subgraph(
        nodes=[
            node
            for node, data in graph_go.nodes(data=True)
            if data["namespace"] == "molecular_function"
        ]
    )
    print(len(graph_go.nodes()))

    # Filter graph by relations
    graph_go = graph_go.edge_subgraph(
        edges={edge for edge in graph_go.edges(keys=True) if edge[2] == "is_a"}
    )
    print(len(graph_go.nodes()))

    # Filter graph by function
    go_name_to_id = {data["name"]: id for id, data in graph_go.nodes(data=True)}
    tmtp_ancestors = nx.ancestors(
        graph_go, go_name_to_id["transmembrane transporter activity"]
    )
    graph_go = graph_go.subgraph(
        tmtp_ancestors | {go_name_to_id["transmembrane transporter activity"]}
    )
    print(len(graph_go.nodes()))

    # Read go-chebi mapping
    df_go_to_chebi = load_df("go_chebi", "data/datasets")

    ## Filter for primary input substrates: 
    df_go_to_chebi = (
        df_go_to_chebi[df_go_to_chebi.relation == "has_primary_input"]
        .reset_index(drop=True)
        .drop("relation", axis=1)
    )
    ## Filter by go terms in graph:
    print(df_go_to_chebi.shape[0])
    df_go_to_chebi = df_go_to_chebi[df_go_to_chebi.go_id.isin(graph_go.nodes())]
    print(df_go_to_chebi.shape[0])

    # dict_go_id_to_chebi_ids = (
    #     df_go_to_chebi_primary[["go_id", "chebi_id"]]
    #     .groupby("go_id")
    #     .apply(lambda x: set(x.chebi_id))
    #     .to_dict()
    # )

    # Filter go-chebi mapping by nodes in filtered go graph

    # Read chebi ontology
    graph_chebi = obonet.read_obo(chebi_obo_path, ignore_obsolete=True)

    ## Filter by stars, or other measure?
    
    ## Filter by substrates

    ## Add ancestors to get full graph

    # Stats and Plots, options for iltering

    ## create overlap heatmap between substrates

    ## Calculate protein annotation overlaps for go and chebi, create plot


    pass


import os

os.getcwd()

# TODO parameters
# TODO refactor

get_substrate_graph(
    {83333},
    go_obo_path="data/raw/ontologies/go.obo",
    chebi_obo_path="data/raw/ontologies/chebi.obo",
)
