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
with open("/home/ad/gene_ontology/goa_uniprot_all_filtered.gaf", "w") as outfile:
    outfile.write("\t".join(outfile_header) + "\n")
    with open("/home/ad/gene_ontology/goa_uniprot_all.gaf") as infile:
        for line in infile:
            if not line.startswith("UniProtKB"):
                # also takes care of comments
                continue
            items = line.split("\t")
            outfile.write(
                "\t".join([items[i] for i in (1,3,4,6,8,13)]+["\n"])
            )

print("done.")
