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
    max_sequence_evidence_code: int = 1,
    additional_proteins: set = None,
    anatomical_entities_whitelist: set = None,
    remove_proteins_without_gene_names:bool=True,
):
    # First, get all sequences with filtering criteria:
    df_sequences = get_sequence_dataset(
        datasets_path=datasets_path,
        organism_ids=organism_ids,
        swissprot_only=swissprot_only,
        max_sequence_evidence_code=max_sequence_evidence_code,
        additional_proteins=additional_proteins,
        remove_proteins_without_gene_names=remove_proteins_without_gene_names
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
    # Filter sequences for those with transporter go annotations
    df_sequences = df_sequences[df_sequences.index.isin(df_uniprot_goa.Uniprot)]
    
    # Filter for cellular components
    if anatomical_entities_whitelist:
        df_uniprot_goa_anatomical_entity = get_go_annotations_subset(
            datasets_path=datasets_path,
            root_go_term="cellular anatomical entity",
            inner_go_relations={"is_a"},
            namespaces_keep={"cellular_component"},
            proteins_subset=set(df_sequences.index),
            go_protein_qualifiers_filter_set={"located_in", "is_active_in"},
            annotations_evidence_codes_remove={"IEA"} if exclude_iea_go_terms else None,
        )
        anatomical_entities_protein_subset = df_uniprot_goa_anatomical_entity[
            df_uniprot_goa_anatomical_entity.go_term_ancestor.isin(
                anatomical_entities_whitelist
            )
        ].Uniprot.unique()

        df_sequences = df_sequences[
            df_sequences.index.isin(anatomical_entities_protein_subset)
        ]

    # Get chebi terms associated with go terms. Get them for ancestors, since go_id is subset of that.
    df_go_chebi = get_go_chebi_annotations(
        dataset_path=datasets_path,
        go_ids_subset=set(df_uniprot_goa.go_id_ancestor),
        go_chebi_relations_subset={"has_primary_input", "has_participant"},
        filter_by_3star=False,  
        add_ancestors=True,
        molecules_only=True
    )
    return df_sequences, df_uniprot_goa, df_go_chebi
