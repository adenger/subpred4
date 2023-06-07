from copy import deepcopy
from subpred.util import load_df
import obonet
import networkx as nx


def get_protein_dataset(
    datasets_folder_path: str,
    organism_ids: set = set(),
    evidence_at_protein_level=True,
    reviewed: bool = True,
):
    df_uniprot = load_df("uniprot", datasets_folder_path)
    if organism_ids:
        df_uniprot = df_uniprot[df_uniprot.organism_id.isin(organism_ids)]
    if evidence_at_protein_level:
        df_uniprot = df_uniprot[df_uniprot.protein_existence == 1]
    if reviewed:
        df_uniprot = df_uniprot[df_uniprot.reviewed]
    return df_uniprot


def get_go_annotations(
    datasets_folder_path: str, proteins: set, include_iea: bool = True
):
    df_goa_uniprot = load_df("go", datasets_folder_path)
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


def get_substrate_graph(
    organism_ids: set, datasets_folder_path: str, go_obo_path: str, chebi_obo_path: str
):
    # Read go graph
    graph_go = obonet.read_obo(go_obo_path, ignore_obsolete=True)
    ## Filter graph by organism
    df_uniprot = get_protein_dataset(
        datasets_folder_path=datasets_folder_path,
        organism_ids=organism_ids,
        evidence_at_protein_level=True,
        reviewed=True,
    )
    df_uniprot_goa = get_go_annotations(
        datasets_folder_path=datasets_folder_path,
        proteins=set(df_uniprot.index.tolist()),
        include_iea=True,
    )
    ## Add ancestors to goa, i.e. more abstract terms
    df_uniprot_goa = add_ancestors(df_uniprot_goa=df_uniprot_goa, graph_go=graph_go)
    print(len(graph_go.nodes()))
    graph_go = graph_go.subgraph(nodes=set(df_uniprot_goa.go_id.unique()))
    print(len(graph_go.nodes()))

    ## Filter graph by aspect/namespace
    graph_go = graph_go.subgraph(
        nodes=[
            node
            for node, data in graph_go.nodes(data=True)
            if data["namespace"] == "molecular_function"
        ]
    )
    print(len(graph_go.nodes()))

    ## Filter graph by relations
    graph_go = graph_go.edge_subgraph(
        edges={edge for edge in graph_go.edges(keys=True) if edge[2] == "is_a"}
    )
    print(len(graph_go.nodes()))

    ## Filter graph by function
    go_name_to_id = {data["name"]: id for id, data in graph_go.nodes(data=True)}
    tmtp_ancestors = nx.ancestors(
        graph_go, go_name_to_id["transmembrane transporter activity"]
    )
    graph_go = graph_go.subgraph(
        tmtp_ancestors | {go_name_to_id["transmembrane transporter activity"]}
    )
    print(len(graph_go.nodes()))

    ## Filter goa dataset by go terms in filtered graph
    df_uniprot_goa = df_uniprot_goa[df_uniprot_goa.go_id.isin(graph_go.nodes())]

    # Annotate graph with proteins from dataset
    # TODO map go_id to uniprot set

    # Annotate with chebi primary input substrates
    ## Read go-chebi mapping
    df_go_to_chebi = load_df("go_chebi", datasets_folder_path)

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

    substrates_set = set(df_go_to_chebi.chebi_term.unique())
    print(substrates_set)
    # dict_go_id_to_chebi_ids = (
    #     df_go_to_chebi_primary[["go_id", "chebi_id"]]
    #     .groupby("go_id")
    #     .apply(lambda x: set(x.chebi_id))
    #     .to_dict()
    # )

    # Read chebi ontology
    # graph_chebi = obonet.read_obo(chebi_obo_path, ignore_obsolete=True)

    ## Filter by stars, or other measure?

    ## Filter by substrates

    ## Add ancestors to get full graph

    # Stats and Plots, options for iltering

    ## create overlap heatmap between substrates

    ### merge uniprot to go with go to chebi

    # TODO instead, annotate graph with protein set and substrates

    df_uniprot_go_chebi = df_uniprot_goa.merge(
        df_go_to_chebi, how="left", on="go_id"
    )

    ## Calculate protein annotation overlaps for go and chebi, create plot

    pass


# TODO update go terms, chebi ids?
# TODO add ancestors to go graph?

import os

os.getcwd()

# TODO parameters
# TODO refactor

get_substrate_graph(
    {83333},
    datasets_folder_path="data/datasets",
    go_obo_path="data/raw/ontologies/go.obo",
    chebi_obo_path="data/raw/ontologies/chebi.obo",
)
