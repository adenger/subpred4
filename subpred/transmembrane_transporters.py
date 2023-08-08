from subpred.protein_dataset import get_sequence_dataset
from subpred.go_annotations import get_go_annotations_subset
from subpred.chebi_annotations import get_go_chebi_annotations

# "human": 		9606
# "athaliana":	3702
# "ecoli": 		83333
# "yeast": 		559292

def get_transmembrane_transporter_dataset(
    organism_ids: set = None,
    swissprot_only: bool = False,
    datasets_path: str = "../data/datasets/",
    exclude_iea_go_terms: bool = False,
    max_sequence_evidence_code:int = 1
):
    # First, get all sequences with filtering criteriou:
    df_sequences = get_sequence_dataset(
        datasets_path=datasets_path,
        organism_ids=organism_ids,
        swissprot_only=swissprot_only,
        max_sequence_evidence_code=max_sequence_evidence_code,
    )
    # Get GO annotations from subset of transmembrane transporter go terms
    df_uniprot_goa = get_go_annotations_subset(
        datasets_path=datasets_path,
        root_go_term="transmembrane transporter activity",
        inner_go_relations={"is_a"},
        namespaces_keep={"molecular_function"},
        proteins_subset=set(df_sequences.index),
        go_protein_qualifiers_filter_set={"enables"},
        annotations_evidence_codes_remove={"IEA"} if exclude_iea_go_terms else None,
    )
    # Filter sequences for those with go annotations
    df_sequences = df_sequences[df_sequences.index.isin(df_uniprot_goa.Uniprot)]

    # Get chebi terms associated with go terms. Get them for ancestors, since go_id is subset of that.
    df_go_chebi = get_go_chebi_annotations(
        dataset_path=datasets_path,
        go_ids_subset=set(df_uniprot_goa.go_id_ancestor),
        go_chebi_relations_subset={"has_primary_input", "has_participant"},
        filter_by_3star=False  # TODO
    )
    return df_sequences, df_uniprot_goa, df_go_chebi
