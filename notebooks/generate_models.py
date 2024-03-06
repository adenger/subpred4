
# %%
from subpred.transmembrane_transporters import get_transmembrane_transporter_dataset

# TODO look at screen session

df_sequences, df_uniprot_goa, df_go_chebi = get_transmembrane_transporter_dataset(
    organism_ids=[3702, 9606, 83333, 559292],
    swissprot_only=True,
    datasets_path="../data/datasets/",
    exclude_iea_go_terms=True,
    max_sequence_evidence_code=1,
    remove_proteins_without_gene_names=True,
)
# display(df_sequences)
# display(df_uniprot_goa)
# display(df_go_chebi)

# %%
from subpred.cdhit import cd_hit
cluster_representatives_70 = cd_hit(df_sequences.sequence, identity_threshold=70)

df_sequences = df_sequences[df_sequences.index.isin(cluster_representatives_70)]

assert (df_sequences.reviewed == True).all()
assert (df_sequences.protein_existence == True).all()

df_uniprot_goa = df_uniprot_goa[(df_uniprot_goa.Uniprot.isin(cluster_representatives_70))].reset_index(drop=True)

assert (df_uniprot_goa.evidence_code != "IEA").all()
assert (df_uniprot_goa.qualifier == "enables").all()

df_go_chebi = df_go_chebi[df_go_chebi.go_id.isin(df_uniprot_goa.go_id_ancestor)]

# %%
from subpred.go_redundancy import get_pairwise_test_scores
get_pairwise_test_scores(
    df_sequences, df_uniprot_goa, min_samples_unique=5, exclude_iea_go_terms=True, dataset_name="uniprot"
)

# %%
get_pairwise_test_scores(
    df_sequences, df_uniprot_goa, min_samples_unique=10, exclude_iea_go_terms=True, dataset_name="uniprot"
)

# %%
get_pairwise_test_scores(
    df_sequences, df_uniprot_goa, min_samples_unique=15, exclude_iea_go_terms=True, dataset_name="uniprot"
)

# %%
get_pairwise_test_scores(
    df_sequences, df_uniprot_goa, min_samples_unique=20, exclude_iea_go_terms=True, dataset_name="uniprot"
)
