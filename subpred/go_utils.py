from .ontology import Ontology
import pandas as pd

EVIDENCE_CODE_TO_DESCRIPTION = {
    "IMP": "experimental_evidence",
    "IPI": "experimental_evidence",
    "IEP": "experimental_evidence",
    "IDA": "experimental_evidence",
    "EXP": "experimental_evidence",
    "IGI": "experimental_evidence",
    "HDA": "experimental_evidence_high_throughput",
    "HMP": "experimental_evidence_high_throughput",
    "HTP": "experimental_evidence_high_throughput",
    "HGI": "experimental_evidence_high_throughput",
    "HEP": "experimental_evidence_high_throughput",
    "IBA": "phylogenetically_inferred",
    "IBD": "phylogenetically_inferred",
    "IKR": "phylogenetically_inferred",
    "IRD": "phylogenetically_inferred",
    "ISS": "computational_analysis",
    "ISO": "computational_analysis",
    "ISA": "computational_analysis",
    "ISM": "computational_analysis",
    "IGC": "computational_analysis",
    "RCA": "computational_analysis",
    "NAS": "author_statement",
    "TAS": "author_statement",
    "IC": "curator_statement",
    "ND": "curator_statement",
    "IEA": "electronic_annotation",
}
GO_FILE = "../data/raw/ontologies/go.owl"


class GeneOntology(Ontology):
    def encode_identifier(self, identifier: str):
        return identifier.replace(":", "_", 1)

    def decode_identifier(self, identifier: str):
        return identifier.replace("_", ":", 1)

   



def read_go_uniprot(
    gaf_file_path: str, remove_not: bool = True, db: str = "UniProtKB"
) -> pd.DataFrame:
    """
    The GAF file format: http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
    """
    go_table = pd.read_table(
        gaf_file_path,
        header=None,
        names=[
            "db",
            "db_object_id",
            "db_object_symbol",
            "qualifier",
            "go_id",
            "db_reference",
            "evidence_code",
            "with_or_from",
            "aspect",
            "db_object_name",
            "db_object_synonym",
            "db_object_type",
            "taxon",
            "date",
            "assigned_by",
            "annotation_extension",
            "gene_product_form_id",
        ],
    )
    if db:
        # Subset gene products by dataset. Proteins, RNA, etc.
        go_table = go_table[go_table.db == db]
    if remove_not:
        # Remove annotations that are explicitly not
        go_table = go_table[~go_table.qualifier.str.startswith("NOT")]
    # evidence description for easier understanding of terms
    go_table = go_table.assign(
        evidence_description=go_table.evidence_code.map(EVIDENCE_CODE_TO_DESCRIPTION)
    )
    return go_table
