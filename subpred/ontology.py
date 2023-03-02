from owlready2 import get_ontology, Restriction
from abc import ABC, abstractmethod

# CHEBI_FILE = "../data/raw/ontologies/chebi.owl"

class Ontology(ABC):
    def __init__(
        self, owl_file_path: str, namespace_url: str = "http://purl.obolibrary.org/obo/"
    ):
        self.ontology = get_ontology(owl_file_path).load()
        self.namespace = self.ontology.get_namespace(namespace_url)

    @abstractmethod
    def encode_identifier(self, identifier):
        """
        Rename identifiers. GO identifiers commonly have the shape GO:[0-9]{7},
        except in the namespace, where they have the pattern GO_[0-9]{7}.
        This method can automatically rename the terms, when inherited by a subclass.
        """
        return identifier

    @abstractmethod
    def decode_identifier(self, identifier):
        """
        Rename identifiers. GO identifiers commonly have the shape GO:[0-9]{7},
        except in the namespace, where they have the pattern GO_[0-9]{7}.
        This method can automatically rename the terms, when inherited by a subclass.
        """
        return identifier

    def get_identifier(self, label: str, return_first=True) -> str | list:
        return (
            self.decode_identifier(self.ontology.search_one(label=label).name)
            if return_first
            else [
                self.decode_identifier(identifier.name)
                for identifier in self.ontology.search(label=label)
            ]
        )

    def update_identifer(self, identifier: str):
        # the problem was that deprecated identifiers did not have labels.
        # If a deprecated identifier does have a label then we keep it.
        # Use self.get_class(identifier).deprecated[0] to check for that
        if self.get_label(identifier) != "":
            return identifier
        else:
            cl = self.get_class(identifier)
            assert len(cl.IAO_0100001) == 1
            term_replaced_by = cl.IAO_0100001[0]
            if isinstance(term_replaced_by, str):
                print(term_replaced_by)
                return term_replaced_by
            else:
                return term_replaced_by.name
            # properties_dict = self.get_properties(identifier)
            # new_identifier = properties_dict["term replaced by"].name
            # return self.decode_identifier(new_identifier)

    def get_label(self, identifier: str) -> str:
        labels = self.namespace[self.encode_identifier(identifier)].label
        return "" if len(labels) == 0 else list(labels)[0]

    def __to_set(self, classes) -> set:
        # classes is a generator in case of subclasses and a set for ancestors/descendants
        return {self.decode_identifier(cl.name) for cl in classes if len(cl.label) > 0}

    def get_ancestors(self, identifier: str) -> set:
        return self.__to_set(self.namespace[self.encode_identifier(identifier)].ancestors())

    def get_descendants(self, identifier: str) -> set:
        return self.__to_set(self.namespace[self.encode_identifier(identifier)].descendants())

    def get_children(self, identifier: str) -> set:
        return self.__to_set(self.namespace[self.encode_identifier(identifier)].subclasses())

    def get_parents(
        self,
        identifier: str,
        include_restrictions: bool = False,
        default_relationship: str = "is_a",
    ) -> set | list:
        identifier_enc = self.encode_identifier(identifier)
        classes = {
            cl
            for cl in self.namespace[identifier_enc].is_a
            if not isinstance(cl, Restriction)
        }
        classes_set = self.__to_set(classes)
        if include_restrictions:
            supercl_with_restr = [
                (default_relationship, class_id) for class_id in classes_set
            ]
            restrictions = [
                cl
                for cl in self.namespace[identifier_enc].is_a
                if isinstance(cl, Restriction)
            ]
            assert all(
                [len(restriction.property.label) == 1 for restriction in restrictions]
            )
            supercl_with_restr.extend(
                [
                    (
                        list(restriction.property.label)[0],
                        self.decode_identifier(restriction.value.name),
                    )
                    for restriction in restrictions
                ]
            )
            return supercl_with_restr
        else:
            return classes_set

    def get_class(self, identifier: str):
        # get owlready class object
        return self.namespace[self.encode_identifier(identifier)]

    def get_properties(self, identifier: str) -> dict:
        # get all properties of class.
        properties = {}
        cl = self.get_class(self.encode_identifier(identifier))
        for property in cl.get_class_properties():
            labels = list(property.label)
            # should be checked once, to see if the additional labels are different.
            # not necessary for chebi, go
            # if len(labels) >= 2:
            #     print(f"Warning: more than one label for a property: {labels}")
            label = property.name if len(labels) == 0 else labels[0]
            values = list(getattr(cl, property.name))
            values = values[0] if len(values) == 1 else values
            if label in properties.keys():
                print(f"Warning: label {label} occurred more than once")
            properties[label] = values

        return properties

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


class ChebiOntology(Ontology):
    def encode_identifier(self, identifier: str):
        return identifier.replace(":", "_", 1)

    def decode_identifier(self, identifier: str):
        return identifier.replace("_", ":", 1)