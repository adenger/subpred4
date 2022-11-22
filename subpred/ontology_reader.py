from owlready2 import get_namespace, get_ontology, Restriction

# from collections import frozenset
# TODO convert : to _?

CHEBI_FILE = "../data/raw/ontologies/chebi.owl"

GO_FILE = "../data/raw/ontologies/go.owl"


class Ontology:
    def __init__(
        self, owl_file_path: str, namespace_url: str = "http://purl.obolibrary.org/obo/"
    ):
        self.ontology = get_ontology(owl_file_path).load()
        self.namespace = self.ontology.get_namespace(namespace_url)

    def get_identifier(self, label: str, return_first=True) -> str | list:
        return (
            self.ontology.search_one(label=label).name
            if return_first
            else [identifier.name for identifier in self.ontology.search(label=label)]
        )

    def get_label(self, identifier: str) -> str:
        labels = self.namespace[identifier].label
        return "" if len(labels) == 0 else list(labels)[0]

    def __to_set(self, classes) -> set:
        # classes is a generator in case of subclasses and a set for ancestors/descendants
        return {cl.name for cl in classes if len(cl.label) > 0}

    def get_ancestors(self, identifier: str) -> set:
        return self.__to_set(self.namespace[identifier].ancestors())

    def get_descendants(self, identifier: str) -> set:
        return self.__to_set(self.namespace[identifier].descendants())

    def get_subclasses(self, identifier: str) -> set:
        return self.__to_set(self.namespace[identifier].subclasses())

    def get_superclasses(
        self,
        identifier: str,
        include_restrictions: bool = False,
        default_relationship: str = "is_a",
    ) -> set | list:
        classes = {
            cl
            for cl in self.namespace[identifier].is_a
            if not isinstance(cl, Restriction)
        }
        classes_set = self.__to_set(classes)
        if include_restrictions:
            supercl_with_restr = [
                (default_relationship, class_id) for class_id in classes_set
            ]
            restrictions = [
                cl
                for cl in self.namespace[identifier].is_a
                if isinstance(cl, Restriction)
            ]
            assert all(
                [len(restriction.property.label) == 1 for restriction in restrictions]
            )
            supercl_with_restr.extend(
                [
                    (list(restriction.property.label)[0], restriction.value.name)
                    for restriction in restrictions
                ]
            )
            return supercl_with_restr
        else:
            return classes_set

    def get_class(self, identifier: str):
        # get owlready class object
        return self.namespace[identifier]

    def get_properties(self, identifier: str) -> dict:
        # get all properties of class.
        properties = {}
        cl = self.get_class(identifier)
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
