import pandas as pd
import numpy as np
from typing import List
import re
from pathlib import Path
from subpred.cdhit import cd_hit


SUBSTRATE_KEYWORDS = {
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

KEYWORDS_TRANSPORT_RELATED = {
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

KEYWORDS_LOCATION = {
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
    df = df.assign(tcdb_class=df.tcdb_id.fillna("0.0").apply(lambda x: x[:3]))
    # Only column that is not string
    df.organism_id = df.organism_id.astype(int)
    df.keywords = df.keywords.astype(str)
    return df


def parse_sequences(df: pd.DataFrame, invalid_amino_acids: str):
    # "BJOUXZ"
    match (invalid_amino_acids):
        case "remove_protein":
            df = df[df.sequence.str.match(re.compile("[ACDEFGHIKLMNPQRSTVWY]+"))]
        case "remove_amino_acids":
            df.sequence = df.sequence.str.replace(
                re.compile("[^ACDEFGHIKLMNPQRSTVWY]+"), ""
            )
        case "keep":
            df = df
        case _:
            raise ValueError(
                "Invalid value of invalid_amino_acids:", invalid_amino_acids
            )
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


def annotate_keywords(df: pd.DataFrame):
    df["keywords_transport"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in SUBSTRATE_KEYWORDS]
        )
    )

    df["keywords_transport_related"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in KEYWORDS_TRANSPORT_RELATED]
        )
    )

    df["keywords_location"] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in KEYWORDS_LOCATION]
        )
    )
    return df


def filter_by_keywords(
    df: pd.DataFrame,
    keywords_transport_filter,
    keywords_component_filter,
    keywords_substrate_filter,
    multi_substrate,
):

    # Transport keyword filtering
    substrate_keywords = {kw.strip() for kw in keywords_transport_filter}
    df = df[
        df.keywords_transport_related.str.split(";").apply(
            lambda l: len(set(l) & substrate_keywords) > 0
        )
    ]

    # Compartment keyword filtering
    df = df[
        df.keywords_location.str.split(";").apply(
            lambda l: len(set(l) & {kw.strip() for kw in keywords_component_filter}) > 0
        )
    ]

    # Substrate keyword filtering
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
        raise ValueError("Invalid parameter for multi_substrate")

    return df


def get_keywords_df(df: pd.DataFrame, use_keyword_names: bool = True):
    srs_keywords = (
        df.keywords.str.split(";").explode().str.strip()
        if use_keyword_names
        else df.keyword_ids.str.split(";").explode().str.strip()
    )

    df_keywords = srs_keywords.to_frame(name="keyword").reset_index(drop=False)
    df_keywords = df_keywords[~df_keywords.keyword.isnull()]

    return df_keywords


def get_go_df(df: pd.DataFrame):
    df_go = df.go_terms.str.split(";").explode().str.strip().reset_index(drop=False)
    go_id_pattern = re.compile("\[(GO\:[0-9]{7})\]")
    df_go["go_id"] = df_go.go_terms.str.extract(go_id_pattern)
    df_go["go_term"] = df_go.go_terms.str.replace(go_id_pattern, "").str.strip()
    df_go = df_go.drop("go_terms", axis=1)
    return df_go


# TODO remove test
# TODO prints
# TODO more generalized version of annoatation/filtering?

# TODO idea: make interesting keyword cols independent of filtering
#
#
# accessions_or = df_keywords[df_keywords.keyword.isin(set(keywords_substrate_filter))].Entry.tolist()
# df_new.loc[accessions_or]

# filter for keywords, then check how many matches exist for each protein
# keyword_matches = df_keywords[df_keywords.keyword.isin(set(keywords_keep))].groupby("Entry").apply(len)
# accessions_and = keyword_matches[keyword_matches == len(keywords_keep)].index

# df_new.loc[accessions_and]


def create_dataset(
    input_file: str,
    keywords_classes: List[str] = None,
    keywords_classes_all: List[str] = list(SUBSTRATE_KEYWORDS),
    keywords_filter: List[str] = None,
    multi_substrate: str = "keep",
    outliers: List[str] = None,
    verbose: bool = False,
    tax_ids_filter: List[int] = None,
    sequence_clustering: int = None,
    invalid_amino_acids: str = "remove_protein",
    evidence_code: int = 2,
    force_update: bool = False,
) -> pd.DataFrame:
    """Creates machine learning dataset from Uniprot data

    Args:
        input_file (str): Uniprot custom download, see Makefile
        keywords_filter (List[str]):
            Uniprot keywords to filter for. Only proteins annotated with all keywords are kept
        keywords_classes (List[str]): The class labels to use for the classification task.
            For list of substrates, look at dataset.SUBSTRATE_KEYWORDS
        keywords_classes_all (List[str]): All possible class labels.
            Only used when for the multi_substrate="remove" or multi_substrate="keep".
            Defaults to dataset.SUBSTRATE_KEYWORDS
        multi_substrate (str, optional):
            How to deal with proteins that are annotated with multiple substrates.
            "keep": return all class labels in keywords_classes_all, separated by ";"
            "remove": remove all proteins annotated with more than one class in keywords_classes_all
            "integrate": only keep substrates that are class labels, remove proteins with more than one class label
            Defaults to "keep".
        outliers (List[str], optional):
            List of uniprot accessions to exclude from the dataset. Defaults to None.
        verbose (bool, optional): Print messages about progress. Defaults to False.
        tax_ids_filter (List[int], optional):
            List of organism identifiers to filter dataset with. Only proteins from organisms are kept.
            Defaults to None.
        sequence_clustering (int, optional): Sequence clustering threshold to use with cd-hit.
            Can be 40,50,60,70,80,90 or 100.
            Defaults to None.
        invalid_amino_acids (str, optional): What to do with protein sequences that contain amino acids like "X".
            "remove_protein": Remove entire protein from dataset
            "remove_amino_acids": Keep protein, but remove invalid amino acids from sequence
            "keep": Do nothing
            Defaults to "remove_protein".
        evidence_code (int, optional): Lower value leads to more accurate data, but less samples.
            1: Experimental evidence at protein level.
            2: Above + experimental evidence at transcript level.
            3: Above + inferred by homology
            4: Above + existence predicted
            5: Above + uncertain
            Defaults to 2.
        force_update (bool, optional):
            Read raw data and write pickle, even if pickle already exists. Defaults to False.

    Returns:
        pd.DataFrame: The finished dataset.
    """

    df = read_raw(input_file=input_file, force_update=force_update)

    df = parse_columns(df)

    df = parse_sequences(df, invalid_amino_acids=invalid_amino_acids)

    df = parse_rows(
        df,
        evidence_code=evidence_code,
        tax_ids_filter=tax_ids_filter,
        outliers=outliers,
    )

    # df_go = get_go_df(df)

    df = annotate_keywords(df)

    # ----------------------------------------
    ## Old version:
    # df = df[df.keywords_transport != ""]

    # df = filter_by_keywords(
    #     df,
    #     keywords_transport_filter=keywords_transport_filter,
    #     keywords_component_filter=keywords_component_filter,
    #     keywords_substrate_filter=keywords_substrate_filter,
    #     multi_substrate=multi_substrate,
    # )
    # ---------------------------------------

    df_keywords = get_keywords_df(df)

    if keywords_filter:
        # get set of proteins that are annotated with keywords_filter
        keyword_matches = (
            df_keywords[df_keywords.keyword.isin(keywords_filter)]
            .groupby("Uniprot")
            .apply(len)
        )
        proteins_all_keywords = (
            keyword_matches[keyword_matches == len(keywords_filter)].index.unique().values
        )

        # only keep those proteins
        df_keywords = df_keywords[
            df_keywords.Uniprot.isin(proteins_all_keywords)
        ].reset_index(drop=True)

    # if keywords_classes:
        # get proteins where at least one substrate is in keywords classes
    if multi_substrate != "integrate":
        df_classes = df_keywords[df_keywords.keyword.isin(keywords_classes_all)]
    else:
        df_classes = df_keywords

    if multi_substrate == "keep":
        df_classes = (
            df_classes.groupby("Uniprot").keyword.apply(list).str.join(";").to_frame()
        )
    elif multi_substrate == "integrate":
        df_classes = df_classes[df_classes.keyword.isin(keywords_classes)]
        df_classes = df_classes[~df_classes.Uniprot.duplicated(keep=False)]
        df_classes = df_classes.set_index("Uniprot")
    elif multi_substrate == "remove":
        df_classes = df_classes[~df_classes.Uniprot.duplicated(keep=False)]
        df_classes = df_classes[df_classes.keyword.isin(keywords_classes)]
        df_classes = df_classes.set_index("Uniprot")

    df_classes.rename(columns={"keyword": "label"}, inplace=True)

    df = df.join(df_classes, how="inner")

    if sequence_clustering:
        cluster_repr = cd_hit(
            df.sequence, identity_threshold=sequence_clustering, verbose=verbose
        )
        df = df.loc[cluster_repr]

    return df


if __name__ == "__main__":
    # print entire df (for debugging)
    pd.set_option("expand_frame_repr", False)
    # test
    outliers = (
        ["Q9HBR0", "Q07837"]
        + ["P76773", "Q47706", "P02943", "P75733", "P69856", "P64550"]
        + ["O81775", "Q9SW07", "Q9FHH5", "Q8S8A0", "Q3E965", "Q3EAV6", "Q3E8L0"]
    )
    from time import time

    t = time()
    df = create_dataset(
        keywords_classes=[
            "Amino-acid transport",
            "Sugar transport",
            "Ion transport",
            "Potassium transport",
        ],
        keywords_filter=["Transmembrane", "Transport"],
        input_file="data/raw/swissprot/uniprot_data_2022_04.tab.gz",
        multi_substrate="integrate",
        verbose=True,
        tax_ids_filter=[3702, 9606, 83333, 559292],
        outliers=outliers,
        sequence_clustering=70,
        evidence_code=2,
        invalid_amino_acids="remove_protein",
        force_update=False
    )
    print(df.shape)
    print(time() - t)
