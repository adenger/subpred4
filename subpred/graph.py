from copy import deepcopy
from subpred.util import load_df
import obonet
import networkx as nx
import pandas as pd

import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np


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
            nx.descendants(graph_go, go_id) | {go_id}
            if go_id in graph_go.nodes()
            else go_id
            for go_id in df_uniprot_goa.go_id
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


def get_id_update_dict(graph, field="alt_id"):
    dict_update_id = dict()
    for node, alt_ids in graph.nodes(data=field):
        if not alt_ids:
            continue
        for alt_id in alt_ids:
            dict_update_id[alt_id] = node

    return dict_update_id


def update_identifiers(
    identifiers: pd.Series, graph: nx.MultiDiGraph, field: str = "alt_id"
):
    dict_update_ids = get_id_update_dict(graph, field=field)
    return identifiers.apply(
        lambda id: dict_update_ids[id] if id in dict_update_ids.keys() else id
    )


def get_filtered_go_graph(
    graph_go: nx.MultiDiGraph,
    root_node: str,
    protein_subset: set = None,
    inter_go_relations: set = {"is_a"},
    aspects: set = {"molecular_function"},
):
    go_name_to_id = {name: id for id, name in graph_go.nodes(data="name")}

    ## Filter graph by protein dataset
    if protein_subset:
        graph_go = graph_go.subgraph(nodes=protein_subset)
    ## Filter graph by aspect/namespace
    graph_go = graph_go.subgraph(
        nodes=[
            node
            for node, data in graph_go.nodes(data=True)
            if data["namespace"] in aspects
        ]
    )

    ## Filter graph by relations
    graph_go = graph_go.edge_subgraph(
        edges={
            edge for edge in graph_go.edges(keys=True) if edge[2] in inter_go_relations
        }
    )

    ## Filter graph by function
    tmtp_ancestors = nx.ancestors(graph_go, go_name_to_id[root_node])
    graph_go = graph_go.subgraph(tmtp_ancestors | {go_name_to_id[root_node]})

    return graph_go.copy()


def get_pairwise_substrate_overlaps(
    dict_chebi_to_uniprot: dict,
    min_overlap: int = 20,
    max_overlap: int = 150,
):
    dict_chebi_to_uniprot_filtered = {
        x: y
        for x, y in dict_chebi_to_uniprot.items()
        if (len(y) >= min_overlap) and (len(y) <= max_overlap)
    }

    ## Calculate protein annotation overlaps for go and chebi, create plot
    sorted_substrate_order = [
        tup[0]
        for tup in sorted(
            [
                (chebi_id, len(proteins))
                for chebi_id, proteins in dict_chebi_to_uniprot_filtered.items()
            ],
            key=lambda x: x[1],
            reverse=True,
        )
    ]
    df_substrate_overlaps = pd.DataFrame(
        0, columns=sorted_substrate_order, index=sorted_substrate_order
    )

    for chebi_id1 in sorted_substrate_order:
        for chebi_id2 in sorted_substrate_order:
            intersection = (
                dict_chebi_to_uniprot_filtered[chebi_id1]
                & dict_chebi_to_uniprot_filtered[chebi_id2]
            )
            overlap = len(intersection)
            df_substrate_overlaps.at[chebi_id1, chebi_id2] = overlap

    return df_substrate_overlaps


def get_go_chebi_mapping(
    datasets_folder_path: str,
    graph_chebi: nx.MultiDiGraph,
    graph_go: nx.MultiDiGraph,
    allowed_relations: set = {"has_primary_input"},
    include_ancestor_chebi_ids=False,
):
    # map the terms of the two graphs to each other, with the allowed relations. If no relation exists then identifier is removed from the mapping
    ## Read go-chebi mapping
    df_go_to_chebi = load_df("go_chebi", datasets_folder_path)

    # update chebi ids
    df_go_to_chebi.chebi_id = update_identifiers(
        identifiers=df_go_to_chebi.chebi_id, graph=graph_chebi
    )

    # update go ids
    df_go_to_chebi.go_id = update_identifiers(
        identifiers=df_go_to_chebi.go_id, graph=graph_go
    )

    ## Filter for primary input substrates:
    df_go_to_chebi = (
        df_go_to_chebi[df_go_to_chebi.relation.isin(allowed_relations)]
        .reset_index(drop=True)
        .drop("relation", axis=1)
    )
    ## Filter by filtered graphs:
    print(df_go_to_chebi.shape[0])
    df_go_to_chebi = df_go_to_chebi[
        df_go_to_chebi.go_id.isin(graph_go.nodes())
    ].reset_index(drop=True)
    print(df_go_to_chebi.shape[0])
    df_go_to_chebi = df_go_to_chebi[
        df_go_to_chebi.chebi_id.isin(graph_chebi.nodes())
    ].reset_index(drop=True)
    print(df_go_to_chebi.shape[0])

    ## Add ancestors
    if include_ancestor_chebi_ids:
        graph_chebi_isa = graph_chebi.edge_subgraph(
            [edge for edge in list(graph_chebi.edges(keys=True)) if edge[2] == "is_a"]
        )
        go_chebi_original_chebi_ids = set(df_go_to_chebi.chebi_id)
        df_go_to_chebi = df_go_to_chebi.drop("chebi_term", axis=1)
        df_go_to_chebi.chebi_id = [
            nx.descendants(graph_chebi_isa, chebi_id) | {chebi_id}
            if chebi_id in graph_chebi_isa.nodes()
            else {chebi_id}
            for chebi_id in df_go_to_chebi.chebi_id
        ]
        # df_go_to_chebi = df_go_to_chebi[~df_go_to_chebi.chebi_id.isnull()]
        # df_go_to_chebi = df_go_to_chebi[df_go_to_chebi.chebi_id.isin(graph_chebi_isa.nodes())]
        df_go_to_chebi = df_go_to_chebi.explode("chebi_id").reset_index(drop=True)
        chebi_id_to_name = {
            id: data["name"] for id, data in graph_chebi.nodes(data=True)
        }
        # filter for original set of molecules (could remove potential substrate classes)
        # TODO remove this once an automatic mehtod has been implemented
        df_go_to_chebi = df_go_to_chebi[
            df_go_to_chebi.chebi_id.isin(go_chebi_original_chebi_ids)
        ]

        df_go_to_chebi = df_go_to_chebi.assign(
            chebi_term=df_go_to_chebi.chebi_id.map(chebi_id_to_name)
        )

    df_go_to_chebi = df_go_to_chebi.drop_duplicates().reset_index(drop=True)
    return df_go_to_chebi


def create_heatmap(
    df_matrix,
    title: str,
    width: int,
    height: int,
    lower_triangle_only: bool = True,
    annotate_values: bool = True,
    values_format: str = ".0f",
    output_path=None,
):
    plt.figure(figsize=(width, height))
    plt.title(title)

    mask = np.triu(np.ones_like(df_matrix), k=1) if lower_triangle_only else None
    g = sns.heatmap(df_matrix, annot=annotate_values, fmt=values_format, mask=mask)

    if output_path:
        fig = g.get_figure()
        fig.savefig(output_path, bbox_inches="tight", dpi=100)

    return g


def sort_by_sample_count(annotations: list, annotation_to_samples: dict):
    # sort list of annotations by number of mapped samples
    nodes_to_count = {
        annotation: len(samples)
        for annotation, samples in annotation_to_samples.items()
        if annotation in annotations
    }
    nodes_to_count = sorted(nodes_to_count.items(), key=lambda x: x[1], reverse=True)
    annotations = [y[0] for y in nodes_to_count]
    return annotations


def add_paths_as_edges(graph, nodes_subset, relations_paths: set = {"is_a"}):
    # creates a subset of the graph, and adds an edge if there is a path that uses the relations_path relations between two nodes.
    # solves the problem caused by creating a subset of a graph that does not contain connecting nodes between two nodes.

    # might not include all possible values
    relations_paths_values = {
        "has_functional_parent",
        "has_parent_hydride",
        "has_part",
        "has_role",
        "is_a",
        "is_conjugate_acid_of",
        "is_conjugate_base_of",
        "is_enantiomer_of",
        "is_substituent_group_from",
        "is_tautomer_of",
    }
    # filter graph
    graph_subset = graph.subgraph(nodes_subset).copy()
    # add edges to show relations. edge if there is a path
    graph_relations = graph.edge_subgraph(
        [e for e in graph.edges(keys=True) if e[2] in relations_paths]
    ).copy()

    for node1 in nodes_subset:
        for node2 in nodes_subset:
            if node1 == node2:
                continue
            if nx.has_path(graph_relations, node1, node2) and not nx.has_path(
                graph_subset, node1, node2
            ):
                graph_subset.add_edge(node1, node2)
    return graph_subset


def graph_plot(
    graph_to_plot,
    root_node: str,
    label_name: str = "name",
    # label_count: str = None,
    title: str = "",
    rotation: int = 0,
    # size: tuple = (30, 30),
    width: int = 20,
    height: int = 15,
    undirected_edges: list = None,
    graphviz_layout: str = "neato",
    node_size=10000,
    output_path=None,
):
    graph_to_plot = graph_to_plot.copy()

    plt.figure(3, figsize=(width, height))

    # progs: "neato", ‘dot’, ‘twopi’, ‘fdp’, ‘sfdp’, ‘circo’, neato
    layout = nx.nx_agraph.graphviz_layout(
        graph_to_plot,
        prog=graphviz_layout,
        root=root_node,
    )

    nx.draw(graph_to_plot, layout, node_size=node_size, node_color="white")

    labels_name_dict = dict(graph_to_plot.nodes(data=label_name))
    text = nx.draw_networkx_labels(
        graph_to_plot,
        pos=layout,
        labels=labels_name_dict,  # , verticalalignment="bottom" # TODO center if no counts
    )
    for _, t in text.items():
        t.set_rotation(rotation)
        t.set_rotation_mode("anchor")

    # if label_count:
    #     labels_count_dict = dict(graph_to_plot.nodes(data=label_count))
    #     text = nx.draw_networkx_labels(
    #         graph_to_plot, pos=layout, labels=labels_count_dict, verticalalignment="top"
    #     )
    #     for _, t in text.items():
    #         t.set_rotation(rotation)
    #         t.set_rotation_mode("anchor")

    if undirected_edges:
        nodes_set = set(graph_to_plot.nodes())

        undirected_edges = [
            edge
            for edge in undirected_edges
            if (edge[0] in nodes_set and edge[1] in nodes_set)
        ]
        if len(undirected_edges[0]) > 2:
            graph_to_plot.add_weighted_edges_from(undirected_edges)
        else:
            graph_to_plot.add_edges_from(undirected_edges)
        nx.draw_networkx_edges(
            graph_to_plot,
            layout,
            edgelist=undirected_edges,
            alpha=0.8,
            style="dashed",
            arrows=False,
            edge_color="black",
        )
        if len(undirected_edges[0]) > 2:
            # labels present
            nx.draw_networkx_edge_labels(
                graph_to_plot,
                layout,
                edge_labels={(edge[0], edge[1]): edge[2] for edge in undirected_edges},
                font_color="black",
                alpha=0.8,
            )

    plt.title(title)

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", dpi=100)

    return plt.gcf()


def preprocess_data(
    organism_ids: set,
    datasets_folder_path: str,
    root_node: str = "transmembrane transporter activity",
):
    df_uniprot = get_protein_dataset(
        datasets_folder_path=datasets_folder_path,
        organism_ids=organism_ids,
        evidence_at_protein_level=True,
        reviewed=True,
    )

    graph_go = load_df("go_obo", folder_path=datasets_folder_path)#obonet.read_obo(go_obo_path, ignore_obsolete=True)

    df_uniprot_goa = get_go_annotations(
        datasets_folder_path=datasets_folder_path,
        proteins=set(df_uniprot.index.tolist()),
        include_iea=True,
    )

    # update go ids
    df_uniprot_goa.go_id = update_identifiers(df_uniprot_goa.go_id, graph_go)

    ## Add ancestors to goa, i.e. more abstract terms
    df_uniprot_goa = add_ancestors(df_uniprot_goa=df_uniprot_goa, graph_go=graph_go)
    print(len(graph_go.nodes()))

    # Filtering GO graph
    graph_go = get_filtered_go_graph(
        graph_go=graph_go,
        root_node=root_node,
        protein_subset=set(df_uniprot_goa.go_id.unique()),
        inter_go_relations={"is_a"},
        aspects={"molecular_function"},
    )

    # filter chebi
    graph_chebi = load_df("chebi_obo", folder_path=datasets_folder_path)#obonet.read_obo(chebi_obo_path, ignore_obsolete=True)
    ## Filter by manually annotated entries
    print(len(graph_chebi.nodes()))
    graph_chebi = graph_chebi.subgraph(
        [x for x, data in graph_chebi.nodes(data=True) if "3_STAR" in data["subset"]]
    )
    print(len(graph_chebi.nodes()))
    return df_uniprot, df_uniprot_goa, graph_go, graph_chebi


def get_substrate_matrix(
    datasets_folder_path: str,
    graph_chebi,
    graph_go,
    df_uniprot_goa,
    min_overlap=20,
    max_overlap="half",
    include_ancestor_chebi_ids=False,
):
    # GO-Chebi mapping
    df_go_to_chebi = get_go_chebi_mapping(
        datasets_folder_path=datasets_folder_path,
        graph_chebi=graph_chebi,
        graph_go=graph_go,
        include_ancestor_chebi_ids=include_ancestor_chebi_ids,
    )

    # Mapping Uniprot to GO to Chebi
    df_uniprot_go_chebi = df_uniprot_goa[["Uniprot", "go_id"]].merge(
        df_go_to_chebi, how="inner", on="go_id"
    )

    df_uniprot_go_chebi = df_uniprot_go_chebi.drop_duplicates().reset_index(drop=True)

    dict_chebi_to_uniprot = (
        df_uniprot_go_chebi[["chebi_id", "Uniprot"]]
        .groupby("chebi_id")
        .apply(lambda x: set(x.Uniprot))
        .to_dict()
    )

    match max_overlap:
        case None:
            max_overlap = 1e6
        case "half":  # number of proteins /2
            max_overlap = df_uniprot_go_chebi.Uniprot.unique().shape[0] // 2

    # Calculating pairwise overlaps for substrates with more then n substrates
    df_substrate_overlaps = get_pairwise_substrate_overlaps(
        dict_chebi_to_uniprot=dict_chebi_to_uniprot,
        min_overlap=min_overlap,
        max_overlap=max_overlap,
    )
    chebi_id_to_name = {id: data["name"] for id, data in graph_chebi.nodes(data=True)}
    df_substrate_overlaps.columns = df_substrate_overlaps.columns.map(chebi_id_to_name)
    df_substrate_overlaps.index = df_substrate_overlaps.index.map(chebi_id_to_name)
    return df_substrate_overlaps, dict_chebi_to_uniprot


def get_graph_plot(
    df_substrate_overlaps,
    dict_chebi_to_uniprot,
    graph_chebi,
    title: str,
    graph_output_path=None,
    node_size: int = 10000,
    width: int = 20,
    height: int = 15,
):
    # sort nodes by number of samples
    chebi_name_to_id = {data["name"]: id for id, data in graph_chebi.nodes(data=True)}
    substrates = df_substrate_overlaps.index.map(chebi_name_to_id)
    substrates = sort_by_sample_count(
        annotations=substrates, annotation_to_samples=dict_chebi_to_uniprot
    )
    graph_chebi_heatmap = add_paths_as_edges(
        graph=graph_chebi,
        nodes_subset=substrates,
        relations_paths={"is_a", "is_tautomer_of"},
    )

    return graph_chebi_heatmap, graph_plot(
        graph_chebi_heatmap,
        chebi_name_to_id["chemical entity"],
        title=title,
        node_size=node_size,
        output_path=graph_output_path,
        width=width,
        height=height,
    )


# if __name__ == "__main__":
#     organism_ids = {83333}
#     datasets_folder_path = "data/datasets"
#     go_obo_path = "data/raw/ontologies/go.obo"
#     chebi_obo_path = "data/raw/ontologies/chebi.obo"
#     heatmap_output_path = "plots/heatmap_ecoli.png"
#     graph_output_path = "plots/graph_ecoli.png"

#     ################
#     # Read datasets
#     ################

#     # Uniprot, GOA, GO, Chebi
#     df_uniprot, df_uniprot_goa, graph_go, graph_chebi = preprocess_data(
#         organism_ids=organism_ids,
#         datasets_folder_path=datasets_folder_path,
#         go_obo_path=go_obo_path,
#         chebi_obo_path=chebi_obo_path,
#     )

#     ###################
#     # Substrate matrix
#     ###################

#     df_substrate_overlaps, dict_chebi_to_uniprot = get_substrate_matrix(
#         datasets_folder_path=datasets_folder_path,
#         graph_chebi=graph_chebi,
#         graph_go=graph_go,
#         df_uniprot_goa=df_uniprot_goa,
#     )

#     ###############
#     # Heatmap
#     ##############

#     create_heatmap(
#         df_matrix=df_substrate_overlaps,
#         title="Substrate molecular species overlaps for substrates with 20 or more transport proteins",
#         width=15,
#         height=10,
#         lower_triangle_only=True,
#         output_path=heatmap_output_path,
#     )

#     ###########
#     # Graph
#     ###########

#     get_graph_plot(
#         df_substrate_overlaps=df_substrate_overlaps,
#         dict_chebi_to_uniprot=dict_chebi_to_uniprot,
#         graph_chebi=graph_chebi,
#     )
# # if __name__ == "__main__":
# #     for organism_id_individual in [3702, 9606, 83333, 559292]:
# #         organism_ids = {organism_id_individual}
# #         organism_ids_str = "+".join([str(i) for i in organism_ids])
# #         print(organism_ids_str)
# #         create_plots(
# #             organism_ids=organism_ids,
# #             heatmap_output_path=f"plots/heatmap_{organism_ids_str}.png",
# #             graph_output_path=f"plots/graph_{organism_ids_str}.png",
# #         )
