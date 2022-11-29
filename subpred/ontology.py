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
            self.rename_identifier(self.ontology.search_one(label=label).name)
            if return_first
            else [
                self.rename_identifier(identifier.name)
                for identifier in self.ontology.search(label=label)
            ]
        )

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
