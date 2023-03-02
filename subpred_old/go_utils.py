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
