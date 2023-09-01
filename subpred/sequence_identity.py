from rpy2.robjects import r, packages, StrVector
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import pandas as pd
import numpy as np


def get_pairwise_sequence_identities(
    series_sequences: pd.Series,
    substitution_matrix: str = "BLOSUM62",
    gap_opening: int = 10,
    gap_extension: int = 4,
    identity_calculation_method: str = "PID1",
) -> pd.DataFrame:
    """Uses R packages Biostrings and future.apply for fast, parallel computation
    of pairwise sequence identity scores from optimal global sequence alignments
    calculated with the Needleman-Wunsch algorithm.

    Args:
        series_sequences (pd.Series): Series with identifiers as index and sequences as values
        substitution_matrix (str, optional): PAM or BLOSUM matrix.
            Options: "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
                "PAM30", "PAM40", "PAM70", "PAM120", "PAM250".
            Defaults to "BLOSUM62".
        gap_opening (int, optional): penalty for gap opening during alignment. Defaults to 10.
        gap_extension (int, optional): penalty for gap extension during alignment. Defaults to 4.
        identity_calculation_method (str, optional):
            "PID1":100 * (identical positions) / (aligned positions + internal gap positions)

            "PID2":100 * (identical positions) / (aligned positions)

            "PID3":100 * (identical positions) / (length shorter sequence)

            "PID4":100 * (identical positions) / (average length of the two sequences)

            Defaults to "PID1".
    Returns:
        pd.DataFrame: Sequence identity scores.
    """
    packages.importr("Biostrings")
    packages.importr("future.apply")
    r(
        """
        pairwiseAlignmentIdentityBatch <- function( 
            subject, 
            patterns, 
            ...
        ){
            # to parallelize outer loop
            alignments <- pairwiseAlignment(
                pattern=patterns, 
                subject=subject,
                ...
            )
            return (alignments)  
        }
        pairwiseAlignmentsList <- function(
            sequences,
            type="global", # Needleman-Wunsch. "local" for smith-waterman
            substitutionMatrix="BLOSUM62",
            gapOpening=10, 
            gapExtension=4
        ){
            # TODO would scale better if not re-calculating symmetrical scores
            plan(multicore)  # uses process fork, on windows use multisession (slower)
            alignments_list <- future_lapply(
                sequences, 
                pairwiseAlignmentIdentityBatch, 
                patterns=sequences, 
                type=type, 
                substitutionMatrix=substitutionMatrix, 
                gapOpening=gapOpening, 
                gapExtension=gapExtension
            )
            return (alignments_list)
        }
        calculateIdentitiesMatrix <- function(alignmentList, accessions, pidType="PID1"){
            identities <- lapply(alignmentList, pid, type=pidType)
            df_identities <- as.data.frame(identities, row.names=accessions, col.names=accessions)
            return (df_identities)
        }
        calculateScoresMatrix <- function(alignmentList, accessions){
            scores <- lapply(alignmentList, score)
            df_scores <- as.data.frame(scores, row.names=accessions, col.names=accessions)
            return (df_scores)
        }
        
    """
    )

    pairwise_alignments_list = r["pairwiseAlignmentsList"]
    calculate_identities_matrix = r["calculateIdentitiesMatrix"]
    calculate_scores_matrix = r["calculateScoresMatrix"]
    # parameter "local" does not make sense here, as local alignments usually have sequence identity of around 100%.
    # local would only make sense together with score() instead of pid()
    alignments_list = pairwise_alignments_list(
        sequences=StrVector(series_sequences.values),
        substitutionMatrix=substitution_matrix,
        gapOpening=gap_opening,
        gapExtension=gap_extension,
    )
    df_identity_r = calculate_identities_matrix(
        alignmentList=alignments_list,
        accessions=StrVector(series_sequences.index),
        pidType=identity_calculation_method,
    )
    df_scores_r = calculate_scores_matrix(
        alignmentList=alignments_list, accessions=StrVector(series_sequences.index)
    )

    with (ro.default_converter + pandas2ri.converter).context():
        df_protein_identity = ro.conversion.get_conversion().rpy2py(df_identity_r)
        df_protein_alignment_scores = ro.conversion.get_conversion().rpy2py(df_scores_r)

    return df_protein_identity, df_protein_alignment_scores


def get_pairwise_go_scores(
    go_terms: np.array,
    go_to_proteins_arr: dict,
    df_protein_scores: pd.DataFrame,
    aggr_method: str = "mean",
):
    func = {"mean": np.mean, "median": np.median, "min": np.min, "max": np.max}[
        aggr_method
    ]
    records = list()
    for i, go_term in enumerate(go_terms):
        # protein arrays have to be unique and sorted!
        protein_arr1 = go_to_proteins_arr[go_term]
        for j in range(i, len(go_terms)):
            go_term_other = go_terms[j]
            protein_arr2 = go_to_proteins_arr[go_term_other]
            df_scores_proteinsets = df_protein_scores.loc[protein_arr1, protein_arr2]
            score = func(df_scores_proteinsets.values)
            records.append([go_term, go_term_other, score])
            records.append([go_term_other, go_term, score])
    df_results = pd.DataFrame(data=records, columns=["go_id1", "go_id2", "score"])
    df_results = df_results.drop_duplicates()
    df_results = df_results.pivot(index="go_id1", columns="go_id2", values="score")
    return df_results


def get_pairwise_alignment_scores(
    df_sequences: pd.DataFrame, df_uniprot_goa: pd.DataFrame, exclude_iea: bool
):
    df_uniprot_goa_copy = df_uniprot_goa.copy(deep=True)
    if exclude_iea:
        df_uniprot_goa_copy = df_uniprot_goa_copy[
            df_uniprot_goa_copy.evidence_code != "IEA"
        ]
    # go_ids_unique = sorted(df_uniprot_goa_copy.go_id_ancestor.unique())

    series_sequences = (
        df_sequences.loc[df_uniprot_goa_copy.Uniprot.unique()]
        .sequence.drop_duplicates()
        .sort_index()
    )
    df_protein_identity, df_protein_alignment_scores = get_pairwise_sequence_identities(
        series_sequences=series_sequences
    )
    return df_protein_identity, df_protein_alignment_scores


def get_aggregated_sequence_alignments_go(
    df_uniprot_goa: pd.DataFrame,
    df_protein_scores: pd.DataFrame,
    exclude_iea: bool = True,
    aggr_method: str = "median",
):
    df_uniprot_goa_copy = df_uniprot_goa.copy(deep=True)
    if exclude_iea:
        df_uniprot_goa_copy = df_uniprot_goa_copy[
            df_uniprot_goa_copy.evidence_code != "IEA"
        ]
    go_terms_unique_sorted = np.sort(df_uniprot_goa_copy.go_id_ancestor.unique())
    go_to_proteins_arr = (
        df_uniprot_goa_copy[["Uniprot", "go_id_ancestor"]]
        .groupby("go_id_ancestor")
        .apply(lambda x: np.sort(x.Uniprot.unique().to_numpy()))
    )
    df_go_pairwise_score = get_pairwise_go_scores(
        go_terms=go_terms_unique_sorted,
        go_to_proteins_arr=go_to_proteins_arr,
        df_protein_scores=df_protein_scores,
        aggr_method=aggr_method,
    )

    return df_go_pairwise_score
