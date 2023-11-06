from subpred.util import load_df
from collections import Counter
import networkx as nx

def get_properties_counts(dataset_path):
    graph_chebi = load_df("chebi_obo", folder_path=dataset_path)

    chebi_id_to_properties = {
        chebi_id: properties_list
        for chebi_id, properties_list in graph_chebi.nodes(data="property_value")
        # if properties_list
    }

    properties_counter = Counter()

    for chebi_term, prop_list in chebi_id_to_properties.items():
        if not prop_list:
            properties_counter["NONE"] +=1
            continue
        else:
            properties_counter["ANY"] +=1
        alreadycounted = set()
        for property_str in prop_list:
            property_type = property_str.split()[0].split("/")[-1]
            if property_type not in alreadycounted:  
                # ignore cases where one molecule has the same property twice
                alreadycounted.add(property_type)
                properties_counter[property_type] +=1

    return properties_counter 


def get_id_update_dict(graph):
    id_update_dict = dict()
    for term, alt_ids in graph.nodes(data="alt_id"):
        if not alt_ids:
            id_update_dict[term] = term
            continue
        for alt_id in alt_ids:
            id_update_dict[alt_id] = term
    for term in graph.nodes():
        id_update_dict[term] = term
    return id_update_dict


def get_go_chebi_annotations(
    dataset_path: str,
    go_ids_subset: set = None,
    go_chebi_relations_subset: set = {"has_primary_input", "has_participant"},
    filter_by_3star: bool = False,
    add_ancestors:bool=False
):
    df_go_chebi = load_df("go_chebi", folder_path=dataset_path)
    graph_chebi = load_df("chebi_obo", folder_path=dataset_path)
    graph_go = load_df("go_obo", folder_path=dataset_path)

    chebi_id_update_dict = get_id_update_dict(graph_chebi)
    go_id_update_dict = get_id_update_dict(graph_go)
    df_go_chebi["go_id"] = df_go_chebi.go_id.map(go_id_update_dict)
    df_go_chebi["chebi_id"] = df_go_chebi.chebi_id.map(chebi_id_update_dict)

    # some ids are null, because versions don't match exactly
    df_go_chebi = df_go_chebi[~df_go_chebi.go_id.isnull()]
    df_go_chebi = df_go_chebi[~df_go_chebi.chebi_id.isnull()]

    df_go_chebi = df_go_chebi.rename(columns={"relation": "chebi_go_relation"})

    if go_chebi_relations_subset:
        df_go_chebi = df_go_chebi[
            df_go_chebi.chebi_go_relation.isin(go_chebi_relations_subset)
        ].reset_index(drop=True)

    if go_ids_subset:
        df_go_chebi = df_go_chebi[df_go_chebi.go_id.isin(go_ids_subset)]

    if filter_by_3star:
        chebi_3star = {
            node
            for node, subset in graph_chebi.nodes(data="subset")
            if "3_STAR" in subset
        }
        df_go_chebi = df_go_chebi[df_go_chebi.chebi_id.isin(chebi_3star)]

    if add_ancestors:
        graph_chebi_isa = graph_chebi.edge_subgraph(
            edges=[
                (source, sink, key)
                for source, sink, key in graph_chebi.edges(keys=True)
                if key == "is_a"
            ]
        )
        # add ancestor chebi ids
        df_go_chebi["chebi_id_ancestor"] = df_go_chebi.chebi_id.transform(
            lambda x: set(nx.descendants(graph_chebi_isa, x) | {x})
        )
        df_go_chebi = df_go_chebi.explode("chebi_id_ancestor")
        chebi_id_to_term = {k: v for k, v in graph_chebi.nodes(data="name")}
        df_go_chebi["chebi_term_ancestor"] = df_go_chebi.chebi_id_ancestor.map(chebi_id_to_term)
        df_go_chebi = df_go_chebi.reset_index(drop=True)

    go_id_to_term = {go_id: go_term for go_id, go_term in graph_go.nodes(data="name")}
    df_go_chebi.insert(1, "go_term", df_go_chebi.go_id.map(go_id_to_term))

    df_go_chebi = df_go_chebi.sort_values(
        ["go_id", "chebi_id", "chebi_go_relation"]
    ).reset_index(drop=True)

    return df_go_chebi
