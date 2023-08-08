from subpred.util import load_df
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

    df_go_chebi = df_go_chebi.sort_values(
        ["go_id", "chebi_id", "chebi_go_relation"]
    ).reset_index(drop=True)

    return df_go_chebi



