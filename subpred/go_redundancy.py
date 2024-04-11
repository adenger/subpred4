import itertools
from copy import copy
import pandas as pd
import numpy as np
import networkx as nx
import os
from subpred.go_annotations import get_go_subgraph
from subpred.util import load_df
from subpred.go_prediction import (
    get_model_evaluation_matrix_parallel,
    process_pairwise_eval_results,
)


def get_proteins(go_set: set, go_id_to_proteins: dict):

    protein_set = set()
    for go_term in go_set:
        protein_set.update(go_id_to_proteins[go_term])
    return protein_set


def get_go_subset(
    df_uniprot_goa,
    root_node="GO:0022857",
    min_samples=20,
    max_samples_percentile=None,
    excluded_terms=None,
):

    graph_go = get_go_subgraph(
        graph_go=load_df("go_obo"),
        root_node=root_node,
        keys={"is_a"},
        namespaces={"molecular_function"},
    )

    go_terms_set = set(graph_go.nodes())
    go_terms_set = set(df_uniprot_goa.go_id_ancestor.unique()) & go_terms_set

    go_term_to_sample_count = (
        df_uniprot_goa[["go_id_ancestor", "Uniprot"]]
        .drop_duplicates()
        .groupby("go_id_ancestor")
        .apply(len)
        .to_dict()
    )

    if min_samples:
        go_terms_set = go_terms_set & set(
            {
                go_term
                for go_term, sample_count in go_term_to_sample_count.items()
                if sample_count >= min_samples
            }
        )

    if max_samples_percentile:
        sample_counts = list(go_term_to_sample_count.values())
        sample_count_percentile = np.percentile(sample_counts, max_samples_percentile)
        go_terms_set = go_terms_set & set(
            {
                go_term
                for go_term, sample_count in go_term_to_sample_count.items()
                if sample_count <= sample_count_percentile
            }
        )

    go_terms_set = {go_term for go_term in go_terms_set if go_term != root_node}
    if excluded_terms:
        go_terms_set = go_terms_set - excluded_terms
    go_terms_list = sorted(list(go_terms_set))
    return go_terms_list


def get_pairwise_test_scores(
    df_sequences,
    df_uniprot_goa,
    dataset_name="",
    min_samples_unique=5,
    exclude_iea_go_terms=False,
    cache_folder=""
):
    if dataset_name:
        file_name = f"ml_models_min{min_samples_unique}_{dataset_name}.pickle"
    else:
        file_name = f"ml_models_min{min_samples_unique}.pickle"

    if cache_folder:
        file_name = f"{cache_folder}/{file_name}"

    if not os.path.isfile(file_name):
        # recalculating with lower threshold to include more terms
        pairwise_eval_results = get_model_evaluation_matrix_parallel(
            df_sequences,  # TODO try embeddings here?
            df_uniprot_goa,
            exclude_iea=exclude_iea_go_terms,
            standardize_samples=True,
            multi_output=True,
            min_samples_per_class=20,
            min_unique_samples_per_class=min_samples_unique,
            model_name="pbest_svc_multi",  # kbest_svc_multi, pca_svc_multi
            feature_selection_parameters=[10, 20, 50],  # TODO try other model here?
            n_jobs=-1,
        )
        df_train_results, df_test_results = process_pairwise_eval_results(
            pairwise_eval_results, df_uniprot_goa, convert_go_ids_to_terms=False
        )
        df_test_scores_new = df_test_results
        df_test_scores_new.to_pickle(file_name)
    else:
        df_test_scores_new = pd.read_pickle(file_name)

    return df_test_scores_new


def get_go_id_to_level(go_terms_list, root_node="GO:0022857"):
    graph_go = get_go_subgraph(
        graph_go=load_df("go_obo"),
        root_node=root_node,
        keys={"is_a"},
        namespaces={"molecular_function"},
    )

    go_id_to_level = {
        go_id: len(nx.shortest_path(graph_go, go_id, root_node))
        for go_id in go_terms_list
    }
    return go_id_to_level


def optimize_subset(
    go_terms_list: list,
    df_test_scores: pd.DataFrame,
    go_id_to_proteins: dict,
    go_id_to_level: dict,
    min_coverage: float = 0.8,
    epsilon_f1: float = 0.0,
    nan_value: float = 0.0,
    prefer_abstract_terms: bool = False,
    verbose: bool = True,
    random_seed: int = 1,
):
    total_proteins_before = len(get_proteins(go_terms_list, go_id_to_proteins))
    current_subset = copy(go_terms_list)
    df_test_scores_subset = df_test_scores.loc[go_terms_list, go_terms_list]
    random_generator = np.random.default_rng(seed=random_seed)
    while True:
        go_term_to_f1_when_removed = list()
        for go_term in current_subset:
            next_potential_subset = list(filter(lambda x: x != go_term, current_subset))
            next_coverage = (
                len(get_proteins(set(next_potential_subset), go_id_to_proteins))
                / total_proteins_before
            )
            if next_coverage < min_coverage:
                continue

            # remove go term, save average and min f1 across all pairs
            df_test_scores_subset_new = df_test_scores_subset.loc[
                next_potential_subset, next_potential_subset
            ]

            # Only keep non-diagonal values, replace NaN with 0 for calculating the average
            df_test_scores_subset_new_long = df_test_scores_subset_new.melt(
                ignore_index=False, value_name="f1_score"
            ).reset_index()
            df_test_scores_subset_new_long = df_test_scores_subset_new_long[
                df_test_scores_subset_new_long.pos_label
                != df_test_scores_subset_new_long.neg_label
            ]
            df_test_scores_subset_new_long.f1_score = df_test_scores_subset_new_long.f1_score.transform(
                # means that not enough unique samples (less than MIN_SAMPLES_UNIQUE)
                lambda s: nan_value if np.isnan(s) else s
            )
            f1_when_removed = (
                df_test_scores_subset_new_long.f1_score.mean()
            )  # TODO try min and max, maybe adjust max_value as well
            go_term_to_f1_when_removed.append((go_term, f1_when_removed))

        if len(go_term_to_f1_when_removed) < 1:
            break

        go_term_to_f1_when_removed = sorted(
            go_term_to_f1_when_removed, key=lambda x: x[1], reverse=True
        )

        max_value = go_term_to_f1_when_removed[0][1]
        possible_go_terms_to_remove = [
            go_term
            for go_term, f1_score in go_term_to_f1_when_removed
            if f1_score >= max_value - epsilon_f1
        ]

        # Tie breaker 1: distance to root node. Select term with highest or lowest
        possible_go_terms_to_remove_levels = [
            (go_id, go_id_to_level[go_id]) for go_id in possible_go_terms_to_remove
        ]
        possible_go_terms_to_remove_levels = sorted(
            possible_go_terms_to_remove_levels,
            key=lambda x: x[1],
            reverse=prefer_abstract_terms,
        )
        selected_level = possible_go_terms_to_remove_levels[0][1]
        possible_go_terms_to_remove = [
            k for k, v in possible_go_terms_to_remove_levels if v == selected_level
        ]
        assert (
            possible_go_terms_to_remove[0] == possible_go_terms_to_remove_levels[0][0]
        )

        # tie breaker 2: random sample.
        go_term_to_remove = random_generator.choice(possible_go_terms_to_remove)
        current_subset = list(filter(lambda x: x != go_term_to_remove, current_subset))

        if verbose:
            print(possible_go_terms_to_remove_levels)
            print(max_value)
            print(go_term_to_f1_when_removed)
            print(possible_go_terms_to_remove)
    return current_subset


def count_nans_nondiag(df):
    return np.isnan(df.values.ravel()).sum() - df.shape[0]


def get_subset_eval(
    go_terms_subset: list,
    all_proteins: list,
    go_id_to_proteins: dict,
    df_test_scores_new: pd.DataFrame,
):
    coverage = len(get_proteins(go_terms_subset, go_id_to_proteins)) / len(
        get_proteins(all_proteins, go_id_to_proteins)
    )
    df_scores_subset = df_test_scores_new.loc[go_terms_subset, go_terms_subset]
    mean = df_scores_subset.mean(axis=None, skipna=True)
    median = df_scores_subset.median(axis=None, skipna=True)
    std = df_scores_subset.values.std(
        axis=None, where=~np.isnan(df_scores_subset.values)
    )
    nans = count_nans_nondiag(df_scores_subset)

    return pd.Series(
        [coverage, mean, median, std, nans, len(go_terms_subset)],
        index=["coverage", "mean", "median", "std", "nans", "subset_length"],
    )


def get_go_id_to_proteins(df_uniprot_goa):

    go_id_to_proteins = (
        df_uniprot_goa[["Uniprot", "go_id_ancestor"]]
        .groupby("go_id_ancestor")
        .agg(set)
        .Uniprot.to_dict()
    )
    return go_id_to_proteins


def subset_pipeline(
    df_uniprot_goa,
    df_sequences,
    root_node="GO:0022857",
    min_samples_per_term=20,
    max_samples_percentile: int = None,
    min_unique_samples_per_term=5,
    min_coverage=0.8,
    epsilon_f1=0.0,
    nan_value=0.0,
    prefer_abstract_terms=False,
    verbose=False,
    excluded_terms=None,
    random_seed=1,
    return_scores=True,
    return_baseline_scores=False,
    # exclude_iea_go_terms=False,
    dataset_name="",
    cache_folder="../data/intermediate/notebooks_cache",
    external_subset=None
):
    """Reduces a set of GO terms to a subset according to specified parameters, using a greedy algorithm.

    Args:
        df_uniprot_goa (pd.DataFrame): GO annotations from get_transmembrane_transporter_dataset or get_go_annotations_subset
        df_sequences (pd.DataFrame): Proteins from from get_transmembrane_transporter_dataset or get_sequence_dataset
        root_node (str, optional): Only GO terms that are ancestors of this root node, with is_a relations, are kept in the dataset. 
            Defaults to "GO:0022857", which corresponds to "transmembrane transporter activity.
        min_samples_per_term (int, optional): Parameter *n* in the paper. GO terms with less than *n* annotated proteins are removed from the dataset. 
            Defaults to 20.
        max_samples_percentile (int, optional): Parameter *p* in the paper. 
            Removes the GO terms that are in the top *p*th percentile in terms of annotated proteins. 
            Defaults to None.
        min_unique_samples_per_term (int, optional): Parameter *m* from paper. ML models are only calculated for pairs GO terms 
            where each GO term has at least *m* unique annotated proteins. A unique protein is a protein that is not part of the intersection set.
            This removes cases where e.g. one GO term is the parent of another, which causes not enough samples to be available for each fold of the 5fCV.
            Higher values can lead to more stable results, but require a lower min_coverage to produce non-redundant results.
            Defaults to 5.
        min_coverage (float, optional): Minimum fraction of the original protein set to be included in the final dataset. Defaults to 0.8.
        epsilon_f1 (float, optional): Two ML evaluation F1 scores a,b are seen as equal if b-epsilon < a < b+epsilon.
            Tries to reduce possible impact of rounding errors on results when sorting the GO terms by their evaluation scores.
            Introduces some deterministic randomness in the algorithm, and can be used to escape local minima. Defaults to 0.0.
        nan_value (float, optional): When calculating the average F1 scores, NaN values can be given a penalty in the calculation,
            making sure that GO terms leading to NaNs are removed as soon as possible. Should be set to 0 or a negative value. Defaults to 0.0.
        prefer_abstract_terms (bool, optional): Preference when choosing between a term closer to the root node and a way further away 
            from the root node that have equal scores. 
            Defaults to False.
        verbose (bool, optional): Print which GO terms are removed during each step. Defaults to False.
        excluded_terms (pd.DataFrame, optional): GO terms that will always be deleted at the start. Defaults to None.
        random_seed (int, optional): The last tie breaker is a random draw with a numpy random generator, this is the random seed for that. Defaults to 1.
        return_scores (bool, optional): If false, only the final subset is returned. If true, returns the subset and the matrix of pairwise F1 scores. Defaults to True.
        return_baseline_scores (bool, optional): Skip all calculations and only return the pairwise F1 scores for the unoptimized dataset. Defaults to False.
        dataset_name (str, optional): Used for distinguishing between the ML results in the cache_folder. 
            Important: Use a different dataset_name when using a different df_sequences/df_uniprot_goa.
            Defaults to "".
        cache_folder (str, optional): Folder for storing the cached pairwise ML results. Defaults to "../data/intermediate/notebooks_cache".
        external_subset: Skip all calculations, use eval functions on this subset instead.
    Returns:
        Optimized GO subset
        Pairwise F1 scores (optional)
    """
    go_id_to_proteins = get_go_id_to_proteins(df_uniprot_goa)

    go_terms_list = get_go_subset(
        df_uniprot_goa=df_uniprot_goa,
        root_node=root_node,
        min_samples=min_samples_per_term,
        max_samples_percentile=max_samples_percentile,
        excluded_terms=excluded_terms,
    )

    go_id_to_level = get_go_id_to_level(go_terms_list, root_node=root_node)

    goa_has_iea_annotations = "IEA" in df_uniprot_goa.evidence_code.unique()

    df_test_scores_new = get_pairwise_test_scores(
        df_sequences=df_sequences,
        df_uniprot_goa=df_uniprot_goa,
        min_samples_unique=min_unique_samples_per_term,
        exclude_iea_go_terms=not goa_has_iea_annotations,
        dataset_name=dataset_name,
        cache_folder=cache_folder,
    )
    if return_baseline_scores:
        scores_before = get_subset_eval(
            go_terms_subset=go_terms_list,
            all_proteins=go_terms_list,
            go_id_to_proteins=go_id_to_proteins,
            df_test_scores_new=df_test_scores_new,
        )
        return go_terms_list, scores_before
    
    if external_subset:
        go_terms_list_len = len(go_terms_list)
        go_terms_list = list(set(external_subset) | set(go_terms_list))
        external_subset_len = len(external_subset)
        external_subset = list(filter(lambda x: x in go_terms_list, external_subset))
        print(len(go_terms_list), go_terms_list_len, len(external_subset), external_subset_len)
        scores_external = get_subset_eval(
            go_terms_subset=external_subset,
            all_proteins=go_terms_list,
            go_id_to_proteins=go_id_to_proteins,
            df_test_scores_new=df_test_scores_new,
        )
        return external_subset, scores_external

    optimized_subset = optimize_subset(
        go_terms_list=go_terms_list,
        df_test_scores=df_test_scores_new,
        go_id_to_proteins=go_id_to_proteins,
        go_id_to_level=go_id_to_level,
        min_coverage=min_coverage,
        epsilon_f1=epsilon_f1,
        nan_value=nan_value,
        prefer_abstract_terms=prefer_abstract_terms,
        verbose=verbose,
        random_seed=random_seed,
    )
    scores_after = get_subset_eval(
        go_terms_subset=optimized_subset,
        all_proteins=go_terms_list,
        go_id_to_proteins=go_id_to_proteins,
        df_test_scores_new=df_test_scores_new,
    )

    if return_scores:
        return optimized_subset, scores_after

    return optimized_subset
