from .compositions import calculate_aac, calculate_paac
from .pssm import calculate_pssm_feature
import pandas as pd

class FeatureGenerator:
    def __init__(self):
    # TODO constants/paths here?
        pass

    def generate_features(
        self,
        sequences: pd.Series,
        feature_types: list = [
            "AAC",
            "PAAC",
            "PSSM_50_1",
            "PSSM_50_3",
            "PSSM_90_1",
            "PSSM_90_3",
        ],
        **kwargs,
    ):
        df_features = []
        for feature_type in set(feature_types):
            if feature_type == "AAC":
                df_features.append(self.__generate_aac(sequences))
            elif feature_type == "PAAC":
                df_features.append(self.__generate_paac(sequences))
            elif feature_type.startswith("PSSM"):
                _, ur_cluster_threshold, iterations = feature_type.split("_")
                df_features.append(
                    self.__generate_pssms(
                        sequences, int(ur_cluster_threshold), int(iterations), **kwargs
                    )
                )
            else:
                raise ValueError(f"Unknown feature type: {feature_type}")

        df_features = pd.concat(df_features, axis=1)
        return df_features

    def __generate_aac(self, sequences: pd.Series):
        return calculate_aac(sequences)

    def __generate_paac(self, sequences: pd.Series):
        return calculate_paac(sequences)

    def __generate_pssms(
        self,
        sequences: pd.Series,
        uniref_cluster_threshold: int,
        psiblast_iterations: int,
        # tmp_folder: str,
        psiblast_executable: str = "psiblast",
        psiblast_threads: int = 4,
    ):
        # TODO do something about hardcoding?
        df_pssm = calculate_pssm_feature(
            sequences,
            tmp_folder="../data/intermediate/blast/pssm_uniref{}_{}it".format(
                uniref_cluster_threshold, psiblast_iterations
            ),
            blast_db="../data/raw/uniref/uniref{}/uniref{}.fasta".format(
                uniref_cluster_threshold, uniref_cluster_threshold
            ),
            iterations=psiblast_iterations,
            psiblast_executable=psiblast_executable,
            psiblast_threads=psiblast_threads,
            verbose=False,
        )
        df_pssm = df_pssm.rename(
            columns=lambda c: f"PSSM_{uniref_cluster_threshold}_{psiblast_iterations}__{c}"
        )

        return df_pssm
