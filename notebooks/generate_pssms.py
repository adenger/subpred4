# TODO delete
from subpred.util import load_df
from subpred.graph import preprocess_data, get_substrate_matrix
from subpred.pssm import calculate_pssm_feature
from subpred.compositions import calculate_aac, calculate_paac
import pandas as pd
from subpred.cdhit import cd_hit

def get_classification_task(
    organism_ids: set,
    labels: set,
    clustering_threshold: int = None,
    dataset_folder_path: str = "../data/datasets",
) -> pd.DataFrame:
    # TODO handling for multi-substrate
    # TODO ability to use go terms or chebi terms (compare sample count, performance)

    (
        df_uniprot,
        df_uniprot_goa,
        graph_go_filtered,
        graph_chebi_filtered,
    ) = preprocess_data(
        organism_ids=organism_ids, datasets_folder_path=dataset_folder_path
    )
    # TODO go through method code
    df_substrate_overlaps, dict_chebi_to_uniprot = get_substrate_matrix(
        datasets_folder_path=dataset_folder_path,
        graph_chebi=graph_chebi_filtered,
        graph_go=graph_go_filtered,
        df_uniprot_goa=df_uniprot_goa,
        min_overlap=0,
        max_overlap=int(1e6),
    )
    assert df_substrate_overlaps.shape[0] == len(dict_chebi_to_uniprot.keys())
    chebi_name_to_term = {
        name: term for term, name in graph_chebi_filtered.nodes(data="name")
    }
    chebi_term_to_name = {
        term: name for term, name in graph_chebi_filtered.nodes(data="name")
    }
    molecule_counts = {
        chebi_term_to_name[term]: len(proteins)
        for term, proteins in dict_chebi_to_uniprot.items()
    }
    # print(sorted(molecule_counts.items(), key=lambda item: item[1], reverse=True))

    protein_to_label = list()
    for label in labels:
        label_proteins = dict_chebi_to_uniprot[chebi_name_to_term[label]]
        for protein in label_proteins:
            protein_to_label.append([protein, label])

    df_labels = pd.DataFrame.from_records(
        protein_to_label, columns=["Uniprot", "label"], index="Uniprot"
    )

    df_labels = df_labels[~df_labels.index.duplicated()]  # TODO series?
    # print(df_labels.label.value_counts())
    df_sequences = df_uniprot.loc[df_labels.index].sequence.to_frame()
    # print("number of sequences", df_sequences.shape[0])
    if clustering_threshold:
        cluster_representatives = cd_hit(
            df_sequences.sequence, identity_threshold=clustering_threshold
        )
        # print(cluster_representatives)
        df_sequences = df_sequences.loc[cluster_representatives]
        df_labels = df_labels.loc[cluster_representatives]
    return pd.concat([df_sequences, df_labels], axis=1)

def get_features(series_sequences:pd.Series):
    # df_aac = calculate_aac(series_sequences)
    # df_paac = calculate_paac(series_sequences)
    df_pssm_50_1 = calculate_pssm_feature(
        series_sequences,
        tmp_folder="../data/intermediate/blast/pssm_uniref50_1it",
        blast_db="../data/raw/uniref/uniref50/uniref50.fasta",
        iterations=1,
        psiblast_threads=-1,
        verbose=True,
        feature_name="PSSM_50_1"
    )
    df_pssm_50_3 = calculate_pssm_feature(
        series_sequences,
        tmp_folder="../data/intermediate/blast/pssm_uniref50_3it",
        blast_db="../data/raw/uniref/uniref50/uniref50.fasta",
        iterations=3,
        psiblast_threads=-1,
        verbose=True,
        feature_name="PSSM_50_3"
    )
    df_pssm_90_1 = calculate_pssm_feature(
        series_sequences,
        tmp_folder="../data/intermediate/blast/pssm_uniref90_3it",
        blast_db="../data/raw/uniref/uniref90/uniref90.fasta",
        iterations=1,
        psiblast_threads=-1,
        verbose=True,
        feature_name="PSSM_90_1"
    )
    df_pssm_90_3 = calculate_pssm_feature(
        series_sequences,
        tmp_folder="../data/intermediate/blast/pssm_uniref90_3it",
        blast_db="../data/raw/uniref/uniref90/uniref90.fasta",
        iterations=3,
        psiblast_threads=-1,
        verbose=True,
        feature_name="PSSM_90_3"
    )
    df_features = pd.concat(
        [
            # df_aac,
            # df_paac,
            df_pssm_50_1,
            df_pssm_50_3,
            df_pssm_90_1,
            df_pssm_90_3,
        ], axis=1
    )
    return df_features

dataset_name_to_organism_ids = {
    "human": {9606},
    "athaliana": {3702},
    "ecoli": {83333},
    "yeast": {559292},
}

test_cases = [
    ("athaliana", "potassium(1+)", "calcium(2+)"),
    ("athaliana", "inorganic anion", "inorganic cation"),
    ("athaliana", "carboxylic acid anion", "inorganic anion"),
    ("ecoli", "carbohydrate derivative", "monosaccharide"),
    ("ecoli", "monocarboxylic acid", "amino acid"),
    ("human", "calcium(2+)", "sodium(1+)"),
    ("human", "calcium(2+)", "potassium(1+)"),
    ("human", "sodium(1+)", "potassium(1+)"),
    ("human", "inorganic anion", "inorganic cation"),
    ("yeast", "amide", "amino acid derivative"),
]

dataset_name_to_organism_ids["all"] = {
    list(s)[0] for s in dataset_name_to_organism_ids.values() if len(s) == 1
}

for dataset_name, substrate1, substrate2 in test_cases:
    organism_ids = dataset_name_to_organism_ids[dataset_name]
    df_dataset = get_classification_task(
        organism_ids=organism_ids,
        labels={substrate1, substrate2},
        clustering_threshold=70,
    )

    df_features = get_features(df_dataset.sequence)