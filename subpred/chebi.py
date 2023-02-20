from .ontology import Ontology

class ChebiOntology(Ontology):
    def encode_identifier(self, identifier: str):
        return identifier.replace(":", "_", 1)

    def decode_identifier(self, identifier: str):
        return identifier.replace("_", ":", 1)