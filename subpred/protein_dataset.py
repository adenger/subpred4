from subpred.util import load_df

def get_sequence_dataset(
    datasets_path:str,
    organism_ids: set = None,
    swissprot_only: bool = False,
    max_sequence_evidence_code: int = 2,
):
    df_uniprot = load_df("uniprot",folder_path=datasets_path)
    if swissprot_only:
        df_uniprot = df_uniprot[df_uniprot.reviewed]
    df_uniprot = df_uniprot[df_uniprot.protein_existence <= max_sequence_evidence_code]
    if organism_ids:
        df_uniprot = df_uniprot[df_uniprot.organism_id.isin(organism_ids)]
    df_uniprot = df_uniprot[["sequence", "reviewed", "protein_existence", "organism_id", "protein_names"]].drop_duplicates()
    return df_uniprot
