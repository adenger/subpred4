import pandas as pd

header = [
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
]
outfile_header = ["Uniprot", "qualifier", "go_id", "evidence_code", "aspect", "date"]

current_uniprot_proteins = set(
    pd.read_table(
        "/home/ad/uniprot_2022_04/uniprot_evidence_go_kw.tsv", index_col=0
    ).index.str.strip().to_list()
)

with open("/home/ad/gene_ontology/goa_uniprot_all_filtered3.gaf", "w") as outfile:
    outfile.write("\t".join(outfile_header) + "\n")
    with open("/home/ad/gene_ontology/goa_uniprot_all.gaf") as infile:
        for line in infile:
            if not line.startswith("UniProtKB"):
                # also takes care of comments
                continue
            items = line.split("\t")
            if items[1].strip() in current_uniprot_proteins:
                outfile.write("\t".join([items[i].strip() for i in (1, 3, 4, 6, 8, 13)] + ["\n"]))

print("done.")
