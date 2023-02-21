import seaborn as sns
import pandas as pd

sns.set(rc={"figure.figsize": (15, 12)})
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from subpred.compositions import calculate_aac, calculate_paac
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

from subpred.dataset import create_dataset, SUBSTRATE_KEYWORDS
from subpred.dataset import get_go_df, get_keywords_df, get_tcdb_substrates
from subpred.go_utils import GeneOntology
from subpred.pssm import calculate_pssm_feature


df_uniprot = create_dataset(
    input_file="../data/raw/uniprot/uniprot_2022_05_evidence1-2_nofragments.tsv",
    # keywords_classes = None,
    # keywords_classes_all = SUBSTRATE_KEYWORDS,
    # keywords_filter = None,
    multi_substrate="keep",
    # outliers=outliers,
    verbose=True,
    # tax_ids_filter=[3702, 9606, 83333, 559292],
    # sequence_clustering=70,
    evidence_code=2,
    invalid_amino_acids="remove_amino_acids",
    # gene_names_only = True,
    # force_update=True,
    # remove_sequence_fragments = True,
    # force_update = False,
    tcdb_substrates_file="../data/raw/tcdb/tcdb_substrates.tsv",
    swissprot_only=True,
)
print(df_uniprot.shape)
df_uniprot.head()


go = GeneOntology("../data/raw/ontologies/go.owl")

transmembrane_transporter_go_terms = go.get_descendants(
    go.get_identifier("transmembrane transporter activity")
)

df_go_uniprot = get_go_df(df_uniprot, go)
df_go_uniprot.head()

df_go_uniprot_long = df_go_uniprot.melt(
    id_vars=["Uniprot"],
    value_vars=["go_id", "go_term"],
    var_name="dataset",
    value_name="annotation",
)
# df_go_uniprot_long = df_go_uniprot_long[df_go_uniprot_long.dataset == "go_term"]
df_go_uniprot_long.head()


transmembrane_transporters = set(
    df_go_uniprot_long[
        df_go_uniprot_long.annotation.isin(transmembrane_transporter_go_terms)
    ].Uniprot.unique()
)
sequences_transmembrane_transporters = df_uniprot.sequence[
    df_uniprot.sequence.index.isin(transmembrane_transporters)
]


# df_pssm_all = pd.DataFrame()
for uniref_cluster_threshold in [50,90]:
    for psiblast_iterations in [1,3]:
        df_pssm = calculate_pssm_feature(
            sequences_transmembrane_transporters,
            tmp_folder="../data/intermediate/blast/pssm_uniref{}_{}it".format(
                uniref_cluster_threshold, psiblast_iterations
            ),
            blast_db="../data/raw/uniref/uniref{}/uniref{}.fasta".format(
                uniref_cluster_threshold, uniref_cluster_threshold
            ),
            iterations=psiblast_iterations,
            psiblast_executable="psiblast",
            psiblast_threads=80,
            verbose=True,
        )
        # df_pssm = df_pssm.rename(
        #     columns=lambda c: f"PSSM_{uniref_cluster_threshold}_{psiblast_iterations}__{c}"
        # )

        # df_pssm_all = pd.concat([df_pssm_all, df_pssm], axis=1)
