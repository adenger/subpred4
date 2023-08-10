from subpred.compositions import calculate_aac, calculate_paac
from subpred.pssm import calculate_pssm_feature
import pandas as pd


def calculate_features(series_sequences: pd.Series):
    df_aac = calculate_aac(series_sequences)
    df_paac = calculate_paac(series_sequences)
    df_pssm_50_1 = calculate_pssm_feature(
        series_sequences,
        tmp_folder="../data/intermediate/blast/pssm_uniref50_1it",
        blast_db="../data/raw/uniref/uniref50/uniref50.fasta",
        iterations=1,
        psiblast_threads=-1,
        verbose=False,
        feature_name="PSSM_50_1",
    )
    df_pssm_50_3 = calculate_pssm_feature(
        series_sequences,
        tmp_folder="../data/intermediate/blast/pssm_uniref50_3it",
        blast_db="../data/raw/uniref/uniref50/uniref50.fasta",
        iterations=3,
        psiblast_threads=-1,
        verbose=False,
        feature_name="PSSM_50_3",
    )
    df_pssm_90_1 = calculate_pssm_feature(
        series_sequences,
        tmp_folder="../data/intermediate/blast/pssm_uniref90_3it",
        blast_db="../data/raw/uniref/uniref90/uniref90.fasta",
        iterations=1,
        psiblast_threads=-1,
        verbose=False,
        feature_name="PSSM_90_1",
    )
    df_pssm_90_3 = calculate_pssm_feature(
        series_sequences,
        tmp_folder="../data/intermediate/blast/pssm_uniref90_3it",
        blast_db="../data/raw/uniref/uniref90/uniref90.fasta",
        iterations=3,
        psiblast_threads=-1,
        verbose=False,
        feature_name="PSSM_90_3",
    )
    df_features = pd.concat(
        [
            df_aac,
            df_paac,
            df_pssm_50_1,
            df_pssm_50_3,
            df_pssm_90_1,
            df_pssm_90_3,
        ],
        axis=1,
    )
    return df_features
