from random import choice, seed
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import SelectKBest, SelectPercentile
from sklearn.preprocessing import MultiLabelBinarizer, LabelEncoder
from sklearn.model_selection import (
    StratifiedKFold,
    GridSearchCV,
)
from sklearn.multioutput import MultiOutputClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import f1_score
from subpred.features import calculate_features
from sklearn.metrics import make_scorer
from joblib import Parallel, delayed


def get_label_to_proteins(
    df_uniprot_goa, min_samples_per_class: int, exclude_iea: bool
):
    df_uniprot_goa_view = df_uniprot_goa.copy()
    if exclude_iea:  # TODO maybe just do this at the beginning?
        df_uniprot_goa_view = df_uniprot_goa_view[
            df_uniprot_goa_view.evidence_code != "IEA"
        ]

    df_uniprot_labels = (
        df_uniprot_goa_view[["Uniprot", "go_id_ancestor"]]
        .drop_duplicates()
        .reset_index(drop=True)
        .rename(columns={"go_id_ancestor": "go_id"})
    )
    label_protein_counts = df_uniprot_labels.groupby("go_id").apply(len)
    labels_enough_proteins = set(
        label_protein_counts[
            label_protein_counts >= min_samples_per_class
        ].index.unique()
    )

    df_uniprot_labels = df_uniprot_labels[
        df_uniprot_labels.go_id.isin(labels_enough_proteins)
    ]
    dict_label_to_proteins = (
        df_uniprot_labels.drop_duplicates()
        .reset_index(drop=True)
        .groupby("go_id")
        .apply(lambda x: set(x.Uniprot))
        .to_dict()
    )
    return dict_label_to_proteins


def transform_labels(
    df_features: pd.DataFrame, df_two_labels: pd.DataFrame, multilabel: bool = True
):
    """Transforms string labels into numbers deterministically

    Args:
        df_features (pd.DataFrame): _description_
        df_two_labels (pd.DataFrame): _description_
        multilabel (bool, optional):
            If true, returns 1-hot encoding of sorted labels (allows overlaps).
            if false, returns index of each label in sorted array
            Defaults to True.
    Returns:
        Numpy arrays for sklearn
    """
    if multilabel:
        series_labels = df_two_labels.groupby("Uniprot").apply(lambda x: list(x.label))
    else:
        # series does not contain duplicated proteins, otherwise verify_integrity throws exception
        series_labels = df_two_labels.set_index("Uniprot", verify_integrity=True)
    sample_names = series_labels.index.values
    X = df_features.loc[sample_names].values
    # automatically uses sorted set of labels as reference
    label_encoder = MultiLabelBinarizer() if multilabel else LabelEncoder()
    y = label_encoder.fit_transform(series_labels.values.ravel())
    class_names = label_encoder.classes_
    assert len(class_names) == 2, class_names
    return X, y, class_names, sample_names


def get_classification_tasks(
    df_features: pd.DataFrame,
    dict_label_to_proteins: dict,
    multi_output: bool,
    # min_samples_per_class: int=20,
    min_unique_samples_per_class: int = 15,
):
    classes = sorted(list(dict_label_to_proteins.keys()))

    for i, go_id1 in enumerate(classes):
        set_proteins1 = dict_label_to_proteins[go_id1]
        df_labels1 = pd.DataFrame({"Uniprot": list(set_proteins1), "label": go_id1})
        for j in range(i, len(classes)):
            go_id2 = classes[j]
            if go_id1 == go_id2:
                continue
            set_proteins2 = dict_label_to_proteins[go_id2]
            df_labels2 = pd.DataFrame({"Uniprot": list(set_proteins2), "label": go_id2})
            df_labels = (
                pd.concat([df_labels1, df_labels2])
                .drop_duplicates()
                .sort_values(["Uniprot", "label"])
                .reset_index(drop=True)
            )
            # TODO if enough samples are available, it could also be viable to turn antiport/symport into a separate class
            intersection_set = set_proteins1 & set_proteins2
            unique_proteins1 = set_proteins1 - intersection_set
            unique_proteins2 = set_proteins2 - intersection_set
            if (
                len(unique_proteins1) < min_unique_samples_per_class
                or len(unique_proteins2) < min_unique_samples_per_class
            ):
                # filter out cases where all proteins from one class are also in the other class (e.g. only (1,1) and (1,0))
                continue
            if not multi_output:
                df_labels = df_labels[
                    df_labels.Uniprot.isin(unique_proteins1 | unique_proteins2)
                ]
            yield transform_labels(df_features, df_labels, multilabel=multi_output)


def stratify_multioutput(element):
    if (element == [0, 1]).all():
        return 1
    elif (element == [1, 0]).all():
        return 0
    elif (element == [1, 1]).all():
        return choice([0, 1])
    else:
        raise ValueError("invalid target value in multi-output binary classification")


def get_stratification_array(arr):
    # solves the problem that stratified kfold does not work on multi-output targets
    seed(1)
    y_stratify = np.array([stratify_multioutput(el) for el in arr])
    return y_stratify


def get_model_scores(
    X,
    y,
    model_name: str = "svc_multi",
    multi_output: bool = True,
    n_threads: int = -1,
    class_names: np.array = np.array([]),
    sample_names: np.array = np.array([]),
):
    models = {
        "svc": make_pipeline(
            StandardScaler(),
            SVC(class_weight="balanced"),
        ),
        "svc_multi": make_pipeline(
            StandardScaler(),
            MultiOutputClassifier(SVC(class_weight="balanced")),
        ),
        "pca_svc_multi": make_pipeline(
            StandardScaler(),
            PCA(),
            MultiOutputClassifier(SVC(class_weight="balanced")),
        ),
        "pbest_svc_multi": make_pipeline(
            StandardScaler(),
            MultiOutputClassifier(
                make_pipeline(SelectPercentile(), SVC(class_weight="balanced"))
            ),
        ),
        "kbest_svc_multi": make_pipeline(
            StandardScaler(),
            MultiOutputClassifier(
                make_pipeline(SelectKBest(), SVC(class_weight="balanced"))
            ),
        ),
        "rf": make_pipeline(
            StandardScaler(),
            RandomForestClassifier(),  # TODO test other model.
        ),
    }
    param_grids = {
        "svc_multi": {
            "multioutputclassifier__estimator__C": [0.1, 1, 10],
            "multioutputclassifier__estimator__gamma": ["scale", "auto"],
        },
        "pca_svc_multi": {
            "multioutputclassifier__estimator__C": [0.1, 1, 10],
            "multioutputclassifier__estimator__gamma": ["scale", "auto"],
            "pca__n_components": [10, 20, None],
        },
        "pbest_svc_multi": {
            "multioutputclassifier__estimator__svc__C": [0.1, 1, 10],
            "multioutputclassifier__estimator__svc__gamma": ["scale", "auto"],
            "multioutputclassifier__estimator__selectpercentile__percentile": [
                10,  # 160
                20,  # 320
                50   # 800
            ],  # percentage of dimensions, 10 default
        },
        "kbest_svc_multi": {
            "multioutputclassifier__estimator__svc__C": [0.1, 1, 10],
            "multioutputclassifier__estimator__svc__gamma": ["scale", "auto"],
            "multioutputclassifier__estimator__selectkbest__k": [
                10,
                40,
                80
            ],  # number of dimensions. 10 is default, 40 is sqrt(1600)
        },
        "svc": {
            "svc__C": [0.1, 1, 10],
            "svc__gamma": ["scale", "auto"],
        },
        "rf": {},
    }
    train_scores_label0 = list()
    train_scores_label1 = list()
    test_scores_label0 = list()
    test_scores_label1 = list()

    # randonly assigns proteins in category "both" to one of them, deterministically
    y_stratify = get_stratification_array(y) if multi_output else y
    for train_index, test_index in StratifiedKFold(n_splits=5).split(X, y_stratify):
        X_train, y_train = X[train_index], y[train_index]
        X_test, y_test = X[test_index], y[test_index]

        if multi_output:
            y_train_stratified = get_stratification_array(y_train)
            grid_search_splits_multioutput_statified = [
                (train_index_inner, test_index_inner)
                for train_index_inner, test_index_inner in StratifiedKFold(
                    n_splits=4
                ).split(X_train, y_train_stratified)
            ]

        gs = GridSearchCV(
            estimator=models[model_name],
            param_grid=param_grids[model_name],
            # optimize for average performance across both labels, punish models that only predict one label
            # get individual scores for each label as well.
            scoring={
                "f1_macro": make_scorer(f1_score, average="macro", zero_division=0),
                "f1_macro_label0": make_scorer(
                    f1_score, labels=[0], average="macro", zero_division=0
                ),
                "f1_macro_label1": make_scorer(
                    f1_score, labels=[1], average="macro", zero_division=0
                ),
            },
            refit="f1_macro",
            n_jobs=n_threads,
            cv=grid_search_splits_multioutput_statified if multi_output else 4,
        )
        gs.fit(X_train, y_train)

        # the average of the two labels is the same as the macro-averaged f1 score
        train_score_label0 = gs.cv_results_["mean_test_f1_macro_label0"][gs.best_index_]
        train_score_label1 = gs.cv_results_["mean_test_f1_macro_label1"][gs.best_index_]
        train_scores_label0.append(train_score_label0)
        train_scores_label1.append(train_score_label1)

        test_scores = f1_score(
            y_true=y_test, y_pred=gs.predict(X_test), average=None, labels=[0, 1]
        )
        test_scores_label0.append(test_scores[0])
        test_scores_label1.append(test_scores[1])

    # return instead of yield to be compatible with joblib parallel
    return (
        class_names,
        np.array(train_scores_label0),
        np.array(train_scores_label1),
        np.array(test_scores_label0),
        np.array(test_scores_label1),
        y,
        sample_names,
    )


def get_model_evaluation_matrix_parallel(
    df_sequences: pd.DataFrame,
    df_uniprot_goa: pd.DataFrame,
    exclude_iea: bool = True,
    standardize_samples: bool = True,
    multi_output: bool = True,
    min_samples_per_class: int = 20,
    min_unique_samples_per_class: int = 15,
    model_name="svc_multi",
    n_jobs: int = -1,
    n_jobs_gridsearch: int = 1,
):
    df_features = calculate_features(
        df_sequences.sequence, standardize_samples=standardize_samples
    )
    dict_label_to_proteins = get_label_to_proteins(
        df_uniprot_goa=df_uniprot_goa,
        min_samples_per_class=min_samples_per_class,
        exclude_iea=exclude_iea,
    )

    # if higher sample count or lower number of classes: paralellize svm or cv or gridsearch instead
    results = Parallel(n_jobs=n_jobs)(
        delayed(get_model_scores)(
            X,
            y,
            model_name=model_name,
            multi_output=multi_output,
            n_threads=n_jobs_gridsearch,
            class_names=class_names,
            sample_names=sample_names,
        )
        for X, y, class_names, sample_names in get_classification_tasks(
            df_features=df_features,
            dict_label_to_proteins=dict_label_to_proteins,
            multi_output=multi_output,
            min_unique_samples_per_class=min_unique_samples_per_class,
        )
    )
    return results


def process_pairwise_eval_results(
    pairwise_eval_results: tuple,
    df_uniprot_goa: pd.DataFrame,
    convert_go_ids_to_terms: bool = True,
):
    go_id_to_term = {
        go_id: go_term
        for go_id, go_term in df_uniprot_goa[["go_id_ancestor", "go_term_ancestor"]]
        .drop_duplicates()
        .to_records(index=False)
    }

    records_train = list()
    records_test = list()

    for (
        class_names,
        train_scores_label0,
        train_scores_label1,
        test_scores_label0,
        test_scores_label1,
        y,
        sample_names,
    ) in pairwise_eval_results:
        go_term0, go_term1 = (
            [go_id_to_term[go_id] for go_id in class_names]
            if convert_go_ids_to_terms
            else list(class_names)
        )
        records_train.append([go_term0, go_term1, train_scores_label0.mean()])
        records_train.append([go_term1, go_term0, train_scores_label1.mean()])
        # records_train.append([go_term1,go_term0,train_scores[1]])
        records_test.append([go_term0, go_term1, test_scores_label0.mean()])
        records_test.append([go_term1, go_term0, test_scores_label1.mean()])
    df_train = pd.DataFrame.from_records(
        records_train, columns=["pos_label", "neg_label", "f1_score"]
    ).pivot(index="pos_label", columns="neg_label", values="f1_score")
    # order = df_train.median().sort_values(ascending=False).index
    # df_train = df_train.loc[order, order]
    df_test = pd.DataFrame.from_records(
        records_test, columns=["pos_label", "neg_label", "f1_score"]
    ).pivot(index="pos_label", columns="neg_label", values="f1_score")
    # df_test = df_test.loc[df_test.mean(numeric_only=True).sort_values(ascending=False).index, df_test.mean(numeric_only=True,axis=1).sort_values(ascending=False).index]
    return df_train, df_test


def plot_results_as_heatmap(
    df_results_pairwise: pd.DataFrame,
    title: str = "",
    height: int = 14,
    width: int = 21,
    sort: bool = False,
):
    fig, ax = plt.subplots(figsize=(width, height))
    if sort:
        df_results_pairwise_plot = df_results_pairwise.loc[
            df_results_pairwise.min(skipna=True, axis=1)
            .sort_values(ascending=False)
            .index,
            df_results_pairwise.min(skipna=True, axis=0)
            .sort_values(ascending=False)
            .index,
        ]
    else:
        df_results_pairwise_plot = df_results_pairwise
    g = sns.heatmap(
        df_results_pairwise_plot,
        ax=ax,
        annot=True,
    )
    ax.set_title(title)
    return g


# Example:
# ----------------
# from subpred.go_prediction import (
#     get_model_evaluation_matrix_parallel,
#     process_pairwise_eval_results,
#     plot_results_as_heatmap,
# )

# pairwise_eval_results = get_model_evaluation_matrix_parallel(
#     df_sequences,
#     df_uniprot_goa,
#     exclude_iea=True,
#     standardize_samples=True,
#     multi_output=True,
#     min_samples_per_class=20,
#     min_unique_samples_per_class=20,
#     model_name="svc_multi",
#     n_jobs=-1,
# )

# df_train, df_test = process_pairwise_eval_results(
#     pairwise_eval_results=pairwise_eval_results, df_uniprot_goa=df_uniprot_goa
# )

# plot_results_as_heatmap(
#     df_train, title="average grid search test scores for each label across CV folds"
# )
