import pandas as pd
from typing import List
import re
from pathlib import Path
from subpred.fasta import write_fasta
from subpred.cdhit import cd_hit


def get_clustering_stats(
    df_dataset: pd.DataFrame,
    keyword_types: list = [
        "keywords_transport",
        "keywords_location",
        "keywords_transport_related",
        "tcdb_class",
    ],
    identity_thresholds: list = [40, 50, 60, 70, 80, 90, 100],
    explode: bool = True,
):
    columns = ["identity_threshold", "kw_type", "keyword", "count"]
    records = []
    for identity_threshold in identity_thresholds:
        cluster_representatives = cd_hit(
            df_dataset.sequence, identity_threshold=identity_threshold, verbose=False
        )
        df_clustered = df_dataset.loc[cluster_representatives]
        for keyword_type in keyword_types:
            for kw, count in df_clustered[keyword_type].value_counts().iteritems():
                records.append([identity_threshold, keyword_type, kw, count])
    df_stats_long = pd.DataFrame.from_records(records, columns=columns)
    if explode:
        df_stats_long.keyword = df_stats_long.keyword.str.split(";")
        df_stats_long = df_stats_long.explode("keyword", ignore_index=True)
        df_stats_long = df_stats_long.groupby(
            ["identity_threshold", "kw_type", "keyword"], as_index=False
        ).sum()
    return df_stats_long


def read_raw(input_file: str, force_update: bool = False):
    input_path = Path(input_file)
    pickle_path = Path(input_path.parent, input_path.name + ".pkl")
    print(pickle_path)
    if pickle_path.exists() and not force_update:
        print("Found pickle, reading...")
        df = pd.read_pickle(pickle_path)
    else:
        print("Reading text file...")
        if not force_update:
            print("Did not find pickle, creating new version...")
        else:
            print("Overwriting existing pickle...")
        df = pd.read_table(input_file, index_col=0, dtype=str)
        df.to_pickle(pickle_path)
    return df


def parse_columns(df: pd.DataFrame):
    df = df.rename(
        columns={
            # Columns in old version of REST API
            "Keyword ID": "keyword_ids",
            "Gene ontology IDs": "go_ids",
            "Gene ontology (GO)": "go_terms",
            "Cross-reference (TCDB)": "tcdb_id",
            # Columns in new version of REST API
            "Gene Ontology IDs": "go_ids",
            "Gene Ontology (GO)": "go_terms",
            "TCDB": "tcdb_id",
            "Organism (ID)": "organism_id",
        },
    )
    df.columns = df.columns.map(lambda c: c.lower().replace(" ", "_"))
    # for compatibility with older scripts
    # df["tax_id"] = df["organism_id"]
    if "entry_name" in df.columns:
        # In new version of Uniprot API, there is an additional column
        df = df.drop("entry_name", axis=1)
    df.index = df.index.rename("Uniprot")
    # Remove trailing ";" from single TCDB IDs, without affecting entries that have >1 TCBD IDs or are NaN
    df.tcdb_id = df.tcdb_id.str.split(";").apply(
        lambda x: ";".join([elem.strip() for elem in x if elem != ""])
        if type(x) != float
        else x
    )
    # Only column that is not string
    df.organism_id = df.organism_id.astype(int)
    df.keywords = df.keywords.astype(str)
    return df


def parse_sequences(df: pd.DataFrame, invalid_amino_acids: str):
    match (invalid_amino_acids):
        case "remove_protein":
            df = df[df.sequence.str.match(re.compile("[ACDEFGHIKLMNPQRSTVWY]+"))]
        case "remove_amino_acids":
            df = df[df.sequence.str.replace(re.compile("[^ACDEFGHIKLMNPQRSTVWY]+"), "")]
        case _:
            raise ValueError("Invalid value of invalid_amino_acids")
    df = df[df.fragment.isnull()]
    df = df.drop(["fragment"], axis=1)

    return df


def parse_rows(
    df: pd.DataFrame, evidence_code: int, tax_ids_filter: List[int], outliers: List[str]
):
    df = df[~df.keywords.isnull()]
    # Mostly peptides, apparently. Like Pollen
    df = df[~df.gene_names.isnull()]

    # df = df[
    #     df.protein_existence.isin(
    #         {"Evidence at protein level", "Evidence at transcript level"}
    #     )
    # ]
    assert evidence_code > 0 and evidence_code <= 5
    evidence_levels_filter = [
        "Evidence at protein level",
        "Evidence at transcript level",
        "Inferred from homology",
        "Predicted",
        "Uncertain",
    ][:evidence_code]
    df = df[df.protein_existence.isin(set(evidence_levels_filter))]

    if tax_ids_filter:
        df = df[~df.organism_id.isnull()]
        tax_ids_keep = set(tax_ids_filter)
        df = df[df.organism_id.isin(tax_ids_keep)]
        for tax_id in tax_ids_keep:
            if tax_id not in df.organism_id.values:
                raise RuntimeError(f"No proteins found for tax id {tax_id}")

    if outliers:
        df = df[~df.index.isin(outliers)]

    return df


# TODO refactoring
# TODO evidence code level
# TODO create seperate df for keywords?
# TODO "" instead of NaN from the start: check if there are still nan checks in script
# TODO docstring
# TODO revert imports? (subpred.)


def create_dataset(
    input_file: str,
    keywords_substrate_filter: List[str],
    keywords_component_filter: List[str],
    keywords_transport_filter: List[str],
    multi_substrate: str = "keep",
    outliers: List[str] = None,
    verbose: bool = False,
    tax_ids_filter: List[int] = None,
    output_tsv: str = None,
    output_fasta: str = None,
    sequence_clustering: int = None,
    invalid_amino_acids: str = "remove_protein",
    evidence_code: int = 2,
    force_update: bool = False,
):

    df = read_raw(input_file=input_file, force_update=force_update)
    # df = df.fillna("")

    df = parse_columns(df)

    df = parse_sequences(df, invalid_amino_acids=invalid_amino_acids)

    df = parse_rows(
        df,
        evidence_code=evidence_code,
        tax_ids_filter=tax_ids_filter,
        outliers=outliers,
    )

    # list of go terms and keywords to filter by
    # list of go terms and keywords to remove
    # list of possible classes
    # list of actual classes

    # TODO field for class labels?
    srs_keywords = df.keywords.str.split(";").explode().str.strip()
    srs_keyword_ids = df.keyword_ids.str.split(";").explode().str.strip()

    # keywords and keyword ids are both sorted alphanumerically, they do not match when axis=1.
    df_keywords = pd.concat([srs_keywords, srs_keyword_ids], axis=0)
    df_keywords.columns = ["keyword"]

    print(df_keywords.loc["P50402"])

    srs_go_long = df.go_terms.str.split(";").explode().str.strip()
    go_id_pattern = re.compile("\[(GO\:[0-9]{7})\]")
    srs_go_long_ids = srs_go_long.str.extract(go_id_pattern)
    srs_go_long_terms = srs_go_long.str.replace(go_id_pattern, "").str.strip()
    df_go = pd.concat([srs_go_long_ids, srs_go_long_terms], axis=1)
    df_go.columns = ["go_id", "go_term"]

    print(df_go)

    ######################
    # Keyword annotation #
    ######################

    keywords_transport = {
        "Ion transport",
        "Anion exchange",
        "Protein transport",
        "Sodium/potassium transport",
        "Polysaccharide transport",
        "Bacteriocin transport",
        "Peptide transport",
        "Translocation",
        "Bacterial flagellum protein export",
        "Amino-acid transport",
        "Electron transport",
        "Lipid transport",
        "mRNA transport",
        "Neurotransmitter transport",
        "Oxygen transport",
        "Phosphate transport",
        "Ammonia transport",
        "Phosphonate transport",
        "Viral movement protein",
        "Sulfate transport",
        "Sugar transport",
        "Calcium transport",
        "Cobalt transport",
        "Copper transport",
        "Hydrogen ion transport",
        "Iron transport",
        "Zinc transport",
        "Nickel transport",
        "Potassium transport",
        "Sodium transport",
        "Chloride",
    }
    df["keywords_transport"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in keywords_transport]
        )
    )
    keywords_transport_related = {
        "Antibiotic resistance",
        "Transport",
        "Symport",
        "Antiport",
        "ER-Golgi transport",
        "Ion channel",
        "Calcium channel",
        "Potassium channel",
        "Chloride channel",
        "Sodium channel",
        "Viral ion channel",
        "Voltage-gated channel",
        "Ligand-gated ion channel",
        "Porin",
        "Nuclear pore complex",
        "Respiratory chain",
    }
    df["keywords_transport_related"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in keywords_transport_related]
        )
    )
    keywords_location = {
        "Membrane",
        "Cell membrane",
        "Cell inner membrane",
        "Cell outer membrane",
        "Transmembrane",
        "Nucleus",
        "Mitochondrion",
        "Endoplasmic reticulum",
        "Plastid outer membrane",
        "Plastid inner membrane",
        "Mitochondrion outer membrane",
        "Mitochondrion inner membrane",
        "Postsynaptic cell membrane",
    }
    df["keywords_location"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in keywords_location]
        )
    )

    ################################################
    # Removing proteins without necessary keywords #
    ################################################

    df = df[df.keywords_transport != ""]

    ###############################
    # Transport keyword filtering #
    ###############################
    substrate_keywords = {kw.strip() for kw in keywords_transport_filter}
    df = df[
        df.keywords_transport_related.str.split(";").apply(
            lambda l: len(set(l) & substrate_keywords) > 0
        )
    ]

    #################################
    # Compartment keyword filtering #
    #################################
    df = df[
        df.keywords_location.str.split(";").apply(
            lambda l: len(set(l) & {kw.strip() for kw in keywords_component_filter}) > 0
        )
    ]

    ###############################
    # Substrate keyword filtering #
    ###############################

    keywords_filter_set = {kw.strip() for kw in keywords_substrate_filter}

    if multi_substrate == "keep":
        # Keep protein if at least one of the substrates is in parameter
        df = df[
            df.keywords_transport.str.split(";").apply(
                lambda l: len(set(l) & keywords_filter_set) > 0
            )
        ]
    elif multi_substrate in {"remove", "integrate"}:
        if multi_substrate == "integrate":
            # Remove all other keywords. Only those proteins where exactly one desired keywords is left are kept
            df.keywords_transport = df.keywords_transport.str.split(";").apply(
                lambda kw_list: ";".join(
                    [kw for kw in kw_list if kw in keywords_filter_set]
                )
            )
        # Only keep protein if it is annotated with one substrate, and that substrate is in parameter
        df = df[df.keywords_transport.apply(lambda s: s.strip() in keywords_filter_set)]
    else:
        # Should not happen, handled by argparse
        raise ValueError("Invalid parameter for multi_substrate")

    ########################
    # Sequence clustering  #
    ########################

    if sequence_clustering:
        cluster_repr = cd_hit(
            df.sequence, identity_threshold=sequence_clustering, verbose=verbose
        )
        df = df.loc[cluster_repr]

    ########################
    # TCDB class field     #
    ########################

    df = df.assign(tcdb_class=df.tcdb_id.fillna("0.0").apply(lambda x: x[:3]))

    ########################
    # Writing Output Files #
    ########################

    df_fasta = df[
        [
            "keywords_transport",
            "gene_names",
            "protein_names",
            "tcdb_id",
            "organism_id",
            "sequence",
        ]
    ]

    if output_fasta:
        df_fasta = df[
            [
                "keywords_transport",
                "gene_names",
                "protein_names",
                "tcdb_id",
                "organism_id",
                "sequence",
            ]
        ]
        fasta_data = (
            df_fasta.reset_index()
            .apply(
                lambda row: (
                    (">sp" + "|{}" * 6).format(
                        row.Uniprot,
                        ";".join(row.gene_names.split()),
                        row.organism_id,
                        row.tcdb_id,
                        row.keywords_transport,
                        row.protein_names,
                    ),
                    row.sequence,
                ),
                axis=1,
                result_type="reduce",
            )
            .tolist()
        )
        write_fasta(fasta_file_name=output_fasta, fasta_data=fasta_data)

    df_tsv = df[
        [
            "keywords_transport",
            "keywords_location",
            "keywords_transport_related",
            "gene_names",
            "protein_names",
            "tcdb_id",
            "tcdb_class",
            "organism_id",
            "sequence",
        ]
    ]
    if output_tsv:
        df_tsv.to_csv(output_tsv, sep="\t")
    return df_tsv


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument(
#         "--keywords-substrate",
#         type=str,
#         required=True,
#         help="List of Uniprot keywords, in Quotes and separated by semicolons",
#     )
#     parser.add_argument(
#         "--keywords-component",
#         type=str,
#         required=True,
#         help="List of Uniprot keywords, in Quotes and separated by semicolons",
#     )
#     parser.add_argument(
#         "--keywords-transport",
#         type=str,
#         required=True,
#         help="List of Uniprot keywords, in Quotes and separated by semicolons",
#     )
#     parser.add_argument("--input-file", type=str, required=True)
#     parser.add_argument("--output-tsv", type=str)
#     parser.add_argument("--output-fasta", type=str)
#     parser.add_argument("--output-log", type=str)
#     parser.add_argument(
#         "--tax-ids",
#         type=int,
#         nargs="+",
#         help="tax id(s) to filter for",
#         default=None,
#     )
#     parser.add_argument("--verbose", action="store_true")
#     parser.add_argument(
#         "--multi-substrate",
#         choices=["keep", "remove", "integrate"],
#         default="keep",
#         help="How to treat proteins with multiple substrates. \
#             Remove them, keep them, or integrate them into \
#             single-substrate classes if the other substrates are not in the dataset",
#     )
#     args = parser.parse_args()
#     create_dataset(
#         keywords_substrate_filter=args.keywords_substrate.split(";"),
#         keywords_component_filter=args.keywords_component.split(";"),
#         keywords_transport_filter=args.keywords_transport.split(";"),
#         input_file=args.input_file,
#         multi_substrate=args.multi_substrate,
#         verbose=args.verbose,
#         tax_ids_filter=args.tax_ids,
#         output_tsv=args.output_tsv,
#         output_fasta=args.output_fasta,
#         output_log=args.output_log,
# )


if __name__ == "__main__":
    # test
    outliers = (
        ["Q9HBR0", "Q07837"]
        + ["P76773", "Q47706", "P02943", "P75733", "P69856", "P64550"]
        + ["O81775", "Q9SW07", "Q9FHH5", "Q8S8A0", "Q3E965", "Q3EAV6", "Q3E8L0"]
    )
    df = create_dataset(
        keywords_substrate_filter=["Amino-acid transport", "Sugar transport"],
        keywords_component_filter=["Transmembrane"],
        keywords_transport_filter=["Transport"],
        input_file="data/raw/swissprot/uniprot-reviewed_yes.tab.gz",
        multi_substrate="integrate",
        verbose=True,
        tax_ids_filter=[3702, 9606, 83333, 559292],
        outliers=outliers,
        sequence_clustering=70,
        # force_update=True
    )
