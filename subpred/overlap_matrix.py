import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def get_go_overlap_matrix(df_uniprot_goa: pd.DataFrame, exclude_iea: bool):
    df_uniprot_goa_overlap = df_uniprot_goa.copy(deep=True)
    if exclude_iea:
        df_uniprot_goa_overlap = df_uniprot_goa[df_uniprot_goa.evidence_code != "IEA"]
    go_ids_unique = sorted(df_uniprot_goa_overlap.go_id_ancestor.unique())
    go_to_proteins = (
        df_uniprot_goa_overlap[["Uniprot", "go_id_ancestor"]]
        .groupby("go_id_ancestor")
        .apply(lambda x: set(x.Uniprot))
        .to_dict()
    )
    records = list()
    for go_id1 in go_ids_unique:
        set1 = go_to_proteins[go_id1]
        for go_id2 in go_ids_unique:
            set2 = go_to_proteins[go_id2]
            overlap = len(set1 & set2)
            records.append([go_id1, go_id2, overlap])

    df_go_overlaps = pd.DataFrame(records, columns=["go_id1", "go_id2", "overlap"])
    df_go_overlaps = df_go_overlaps.pivot(
        index="go_id1", columns="go_id2", values="overlap"
    )
    return df_go_overlaps


def get_overlap_plot(
    subset, df_go_overlaps, go_id_to_term, filename=None, shorten_labels=False
):
    # the percentage of proteins annotated with GO term 1 (x axis) that is also annotated with GO term2
    best_subset_overlaps = df_go_overlaps.loc[subset, subset].rename(
        columns=go_id_to_term, index=go_id_to_term
    )
    best_subset_overlaps = best_subset_overlaps.loc[
        sorted(best_subset_overlaps.index), sorted(best_subset_overlaps.columns)
    ]

    best_subset_overlaps = best_subset_overlaps / np.diag(best_subset_overlaps)

    if shorten_labels:
        shorten_labels_map = {
            c: c.replace("transmembrane transporter activity", "TTA")
            for c in best_subset_overlaps.columns
        }
        shorten_labels_map.update(
            {
                "ATPase activity, coupled to transmembrane movement of ions, rotational mechanism":\
                      "ATPase activity ... rotational mechanism (GO:0044769)"
            }
        )
        best_subset_overlaps = best_subset_overlaps.rename(
            columns=shorten_labels_map, index=shorten_labels_map
        )

    tmp = best_subset_overlaps.copy()
    np.fill_diagonal(tmp.values, 0)
    print("median", tmp.median(axis=None).round(2))
    print("mean", tmp.mean(axis=None).round(2))

    g = sns.heatmap(best_subset_overlaps, annot=True, fmt=".2f")
    g.set_xlabel("GO term 1")
    g.set_ylabel("GO term 2")
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches="tight")
    return g


def plot_go_overlap_matrix(df_go_overlaps: pd.DataFrame, df_uniprot_goa: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(30, 20))
    go_id_to_term = {
        go_id: go_term
        for go_id, go_term in df_uniprot_goa[["go_id_ancestor", "go_term_ancestor"]]
        .drop_duplicates()
        .to_records(index=False)
    }
    df_go_overlaps_ge20 = df_go_overlaps.loc[
        np.diag(df_go_overlaps) >= 20, np.diag(df_go_overlaps) >= 20
    ]

    df_go_overlaps_ge20 = df_go_overlaps_ge20.rename(
        columns=go_id_to_term, index=go_id_to_term
    )
    return sns.heatmap(df_go_overlaps_ge20, ax=ax, annot=True)


def get_chebi_overlaps(
    df_go_chebi: pd.DataFrame,
    df_uniprot_goa: pd.DataFrame,
    exclude_iea: bool,
    primary_input_only: bool = True,
):

    df_go_chebi_overlaps = df_go_chebi.copy(deep=True)
    df_uniprot_goa_chebi_overlap = df_uniprot_goa.copy(deep=True)
    if exclude_iea:
        df_uniprot_goa_chebi_overlap = df_uniprot_goa_chebi_overlap[
            df_uniprot_goa_chebi_overlap.evidence_code != "IEA"
        ]
    # df_go_chebi_overlaps["proteins"] = df_go_chebi_overlaps.go_id.map(go_to_proteins)

    if primary_input_only:
        df_go_chebi_overlaps = df_go_chebi_overlaps[
            df_go_chebi_overlaps.chebi_go_relation == "has_primary_input"
        ]

    df_uniprot_go_chebi = (
        pd.merge(
            df_uniprot_goa_chebi_overlap[["Uniprot", "go_id_ancestor"]].rename(
                columns={"go_id_ancestor": "go_id"}
            ),
            df_go_chebi_overlaps[["go_id", "chebi_id"]],
            on="go_id",
            how="inner",
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )
    dict_chebi_to_uniprot = (
        df_uniprot_go_chebi[["Uniprot", "chebi_id"]]
        .drop_duplicates()
        .groupby("chebi_id")
        .apply(lambda x: set(x.Uniprot))
        .to_dict()
    )
    chebi_terms = sorted(dict_chebi_to_uniprot.keys())
    df_chebi_overlaps = pd.DataFrame(
        data=[
            [
                len(dict_chebi_to_uniprot[chebi1] & dict_chebi_to_uniprot[chebi2])
                for chebi1 in chebi_terms
            ]
            for chebi2 in chebi_terms
        ],
        columns=chebi_terms,
        index=chebi_terms,
    )
    return df_chebi_overlaps
