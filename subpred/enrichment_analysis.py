from scipy.stats import hypergeom, rankdata
import numpy as np
import pandas as pd


def enrichment_analysis(
    proteins_reference: set,
    proteins_subset: set,
    annotations: list,
    p_cutoff: float = None,
    min_lfc: float = None,
) -> pd.DataFrame:
    # Takes a protein set a sub-set and calculates p value and fold change.
    assert len(proteins_reference) > 0 and len(proteins_subset) > 0, "set was empty"
    assert all(
        protein in proteins_reference for protein in proteins_subset
    ), "a protein in the subset was not in the reference set."

    annotations_df = pd.DataFrame.from_records(
        annotations, columns=["identifier", "annotation"]
    )
    annotations_reference = annotations_df[
        annotations_df.identifier.isin(proteins_reference)
    ]
    annotations_subset = annotations_df[annotations_df.identifier.isin(proteins_subset)]

    annotation_scores = []
    for annotation in annotations_subset.annotation.unique():
        total_reference = len(proteins_reference)
        annotated_reference = (
            annotations_reference[annotations_reference.annotation == annotation]
            .identifier.unique()
            .shape[0]
        )
        total_subset = len(proteins_subset)
        annotated_subset = (
            annotations_subset[annotations_subset.annotation == annotation]
            .identifier.unique()
            .shape[0]
        )
        expected = (annotated_reference / total_reference) * total_subset
        fold_change = np.log2(annotated_subset / expected)
        percentage_of_annotated = round(annotated_subset / annotated_reference * 100, 2)

        # Scipy naming convention, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html
        # Normal convention: N, K, n, k
        # M proteins in the reference set (i.e. genome)
        M = total_reference
        # n of those are annotated with the annotation
        n = annotated_reference
        # we draw N proteins from M
        N = total_subset
        # k of those are annotated with the annotation
        k = annotated_subset
        # The cumulative distribution function (cdf) gives the probability of drawing k or fewer than k annotated proteins.
        # The cdf is the sum of the probability mass function for all values <=k.
        # The survival function sf is defined 1-cdf and is therefore the probability of getting more than k values.
        # Since we want the probability for k or more values, we have to calculate the sf at k-1.
        p_val = hypergeom.sf(k - 1, M, n, N)

        annotation_scores.append(
            [
                annotation,
                total_reference,
                annotated_reference,
                total_subset,
                annotated_subset,
                expected,
                percentage_of_annotated,
                fold_change,
                p_val,
            ]
        )
    annotation_scores = pd.DataFrame.from_records(
        annotation_scores,
        columns=[
            "annotation",
            "total_reference",
            "annotated_reference",
            "total_subset",
            "annotated_subset",
            "expected",
            "percentage_of_annotated",
            "lfc",
            "p",
        ],
    )

    annotation_scores["p_fdr"] = (
        annotation_scores.p * len(annotation_scores.p) / rankdata(annotation_scores.p)
    )
    annotation_scores.p_fdr = annotation_scores.p_fdr.clip(upper=1.0)

    annotation_scores["p_bonferroni"] = annotation_scores.p.transform(
        lambda p: min(p * annotation_scores.shape[0], 1.0)
    )
    annotation_scores = annotation_scores.sort_values(["p", "lfc"])

    if p_cutoff:
        annotation_scores = annotation_scores[annotation_scores.p_fdr <= p_cutoff]

    if min_lfc:
        annotation_scores = annotation_scores[annotation_scores.lfc >= min_lfc]

    return annotation_scores
