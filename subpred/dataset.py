import pandas as pd
import re
from pathlib import Path
from subpred.cdhit import cd_hit
from urllib.parse import urlencode, quote
from subpred.go_utils import GeneOntology

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

# columns to be added to dataframe
KEYWORD_SETS = {
    "keywords_substrates": SUBSTRATE_KEYWORDS,
    "keywords_transport_related": KEYWORDS_TRANSPORT_RELATED,
    "keywords_location": KEYWORDS_TRANSPORT_RELATED,
}


def __read_raw(input_file: str, force_update: bool = False):
    # does not work if file paths contain "~"
    input_path = Path(input_file)
    pickle_path = Path(input_path.parent, input_path.name + ".pkl")
    if pickle_path.exists() and not force_update:
        print("Found pickle, reading...")
        df = pd.read_pickle(pickle_path)
    else:
        print("Reading text file...")
        if not force_update:
            print("Did not find pickle, creating new version...")
        else:
            print("Overwriting existing pickle...")
        df = pd.read_table(input_file, index_col=0, low_memory=False)
        df.to_pickle(pickle_path)
    return df


def parse_columns(df: pd.DataFrame):
    """Basic cleanup of raw data"""
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


def __parse_sequences(
    df: pd.DataFrame, invalid_amino_acids: str, remove_sequence_fragments: bool
):
    """Handle non-standard amino acids, and sequence fragments"""
    match (invalid_amino_acids):
        case "remove_protein":
            df = df[df.sequence.str.fullmatch(re.compile("[ACDEFGHIKLMNPQRSTVWY]+"))]
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

    if remove_sequence_fragments and "fragment" in df.columns:
        df = df[df.fragment.isnull()]
        df = df.drop(["fragment"], axis=1)

    return df


def __add_keyword_column(df: pd.DataFrame, keyword_set: set, colname: str):
    df[colname] = df.keywords.str.split(";").apply(
        lambda keywords: ";".join(
            [keyword for keyword in keywords if keyword in keyword_set]
        )
    )

    return df


def __filter_protein_evidence(df: pd.DataFrame, evidence_code: int):
    assert evidence_code > 0 and evidence_code <= 5
    evidence_levels_filter = [
        "Evidence at protein level",
        "Evidence at transcript level",
        "Inferred from homology",
        "Predicted",
        "Uncertain",
    ][:evidence_code]
    df = df[df.protein_existence.isin(set(evidence_levels_filter))]
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


def get_tcdb_substrates(df: pd.DataFrame):
    df_tcdb = df.reset_index(drop=False)[["Uniprot", "tcdb_id", "tcdb_substrates"]]
    df_tcdb = df_tcdb[(~df_tcdb.tcdb_id.isnull()) & (~df_tcdb.tcdb_substrates.isnull())]
    df_substrates = (
        df_tcdb.tcdb_substrates.str.split("|")
        .explode()
        .str.split(";", expand=True)
        .rename(columns={0: "chebi_id", 1: "chebi_term"})
    )
    df_substrates = (
        df_tcdb.drop("tcdb_substrates", axis=1)
        .join(df_substrates)
        .drop_duplicates()
        .reset_index(drop=True)
    )
    return df_substrates


def get_go_df(df: pd.DataFrame, go: GeneOntology):
    # df_go = df.go_terms.str.split(";").explode().str.strip().reset_index(drop=False)
    # go_id_pattern = re.compile("\[(GO\:[0-9]{7})\]")
    # df_go["go_id"] = df_go.go_terms.str.extract(go_id_pattern)
    # df_go["go_term"] = df_go.go_terms.str.replace(go_id_pattern, "").str.strip()
    # df_go = df_go.drop("go_terms", axis=1)

    # transform to long df of uniprot/go_id
    df_go = df.go_ids.str.split(";").explode().str.strip().reset_index(drop=False)
    df_go = df_go.rename(columns={"go_ids": "go_id"})
    df_go = df_go[~df_go.go_id.isnull()]

    go_ids_unique = df_go.go_id.unique()
    go_ids_update_dict = {go_id: go.update_identifer(go_id) for go_id in go_ids_unique}

    # faster than replace, same result:
    df_go = df_go.assign(go_id=df_go.go_id.map(go_ids_update_dict))
    assert df_go.go_id[df_go.go_id.isnull()].shape[0] == 0
    # df_go = df_go.replace({"go_id": deprecated_go_map})

    id_to_term = {
        identifier: go.get_label(identifier) for identifier in df_go.go_id.unique()
    }

    # df_go = df_go.assign(go_term = df_go.go_id.apply(go.get_label))
    df_go = df_go.assign(go_term=df_go.go_id.map(id_to_term))
    df_go.head()
    return df_go


def filter_by_keywords(df: pd.DataFrame, keywords_filter: set):
    df_keywords = get_keywords_df(df)
    keyword_matches = (
        df_keywords[df_keywords.keyword.isin(keywords_filter)]
        .groupby("Uniprot")
        .apply(len)
    )
    proteins_all_keywords = set(
        keyword_matches[keyword_matches == len(keywords_filter)].index.unique().values
    )

    df = df[df.index.isin(proteins_all_keywords)]
    return df


def __add_class_labels(
    df: pd.DataFrame,
    keywords_classes: set,
    keywords_classes_all: set,
    multi_substrate: str,
):
    df_classes = get_keywords_df(df)

    if multi_substrate == "keep":
        df_classes = df_classes[df_classes.keyword.isin(keywords_classes_all)]
        df_classes = (
            df_classes.groupby("Uniprot").keyword.apply(list).str.join(";").to_frame()
        )
    elif multi_substrate == "integrate":
        df_classes = df_classes[df_classes.keyword.isin(keywords_classes)]
        df_classes = df_classes[~df_classes.Uniprot.duplicated(keep=False)]
        df_classes = df_classes.set_index("Uniprot")
    elif multi_substrate == "remove":
        df_classes = df_classes[df_classes.keyword.isin(keywords_classes_all)]
        df_classes = df_classes[~df_classes.Uniprot.duplicated(keep=False)]
        df_classes = df_classes[df_classes.keyword.isin(keywords_classes)]
        df_classes = df_classes.set_index("Uniprot")

    df_classes.rename(columns={"keyword": "label"}, inplace=True)

    df = df.join(df_classes, how="inner")
    return df


def create_dataset(
    input_file: str,
    keywords_classes: set = None,
    keywords_classes_all: set = SUBSTRATE_KEYWORDS,
    keywords_filter: set = None,
    multi_substrate: str = "keep",
    outliers: set = None,
    verbose: bool = False,
    tax_ids_filter: set = None,
    sequence_clustering: int = None,
    invalid_amino_acids: str = "remove_protein",
    evidence_code: int = 2,
    gene_names_only: bool = True,
    remove_sequence_fragments: bool = True,
    force_update: bool = False,
    tcdb_substrates_file: str = None,
    swissprot_only: str = True,
) -> pd.DataFrame:
    """Creates machine learning dataset from Uniprot data

    ## Args:
        input_file (str): Uniprot custom download, see Makefile
        keywords_classes (List[str]): The class labels to use for the classification task.
            For list of substrates, look at dataset.SUBSTRATE_KEYWORDS
            Defaults to None
        keywords_classes_all (List[str]): All possible class labels.
            Only used when for the multi_substrate="remove" or multi_substrate="keep".
            Defaults to dataset.SUBSTRATE_KEYWORDS
        keywords_filter (List[str]):
            Uniprot keywords to filter for. Only proteins annotated with all keywords are kept
            Defaults to None
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
        gene_names_only (bool, optional):
            If True, remove all proteins with no associated gene name.
            Useful when using gene annotation data.
            Defaults to True
        remove_sequence_fragments (bool, optional):
            Remove fragmented sequences from dataset.
            Defaults to True
        force_update (bool, optional):
            Read raw data and write pickle, even if pickle already exists. Defaults to False.
        tcdb_substrates_file (str, optional):
            Tab-seperated file with TCDB IDs and substrates, from TCDB website. See makefile for download link.
        swissprot_only (bool, optional):
            Removes proteins from the dataset that have not been manually reviewed

    ## Returns:
        pd.DataFrame: The finished dataset.

    ## Examples:

    ### Download raw data::

        from subprocess import run
        url = __create_uniprot_url()
        filename = "uniprot_data_2022_04.tab.gz"
        run(f'curl "{url}" > "{filename}"', shell=True)

    ### Create transporter dataset::

        from subpred.dataset import create_dataset
        create_dataset(
            input_file="uniprot_data_2022_04.tab.gz",
            keywords_classes={
                "Amino-acid transport",
                "Sugar transport",
                "Ion transport",
                "Potassium transport",
            },
            keywords_filter={"Transmembrane", "Transport"},
            multi_substrate="integrate",
            verbose=True,
            tax_ids_filter={3702, 9606, 83333, 559292},
            outliers={
                "Q9HBR0",
                "Q07837",
                "P76773",
                "Q47706",
                "P02943",
                "P75733",
                "P69856",
                "P64550",
                "O81775",
                "Q9SW07",
                "Q9FHH5",
                "Q8S8A0",
                "Q3E965",
                "Q3EAV6",
                "Q3E8L0",
            },
            sequence_clustering=70,
            evidence_code=2,
            invalid_amino_acids="remove_protein",
            force_update=False,
        )

    ### Add keyword column to dataset::

        from subpred import dataset
        col_name = "new_col"
        keyword_set = {"Transmembrane", "Cell membrane"}
        dataset.KEYWORD_SETS.update({col_name: keyword_set})

    """

    df = __read_raw(input_file=input_file, force_update=force_update)

    df = parse_columns(df)

    if gene_names_only:
        # Mostly peptides, apparently. Like Pollen
        df = df[~df.gene_names.isnull()]

    if swissprot_only and "reviewed" in df.columns:
        df = df[df.reviewed == "reviewed"]

    df = __parse_sequences(
        df,
        invalid_amino_acids=invalid_amino_acids,
        remove_sequence_fragments=remove_sequence_fragments,
    )

    df = __filter_protein_evidence(df, evidence_code=evidence_code)

    if tax_ids_filter:
        df = df[df.organism_id.isin(tax_ids_filter)]

    if outliers:
        df = df[~df.index.isin(outliers)]

    # df_go = get_go_df(df)

    if keywords_filter:
        df = filter_by_keywords(df, keywords_filter=keywords_filter)

    if keywords_classes:
        df = __add_class_labels(
            df=df,
            keywords_classes=keywords_classes,
            keywords_classes_all=keywords_classes_all,
            multi_substrate=multi_substrate,
        )

    if tcdb_substrates_file:
        df_substrates = pd.read_table(tcdb_substrates_file, header=None)
        df_substrates.columns = ["tcdb_id", "tcdb_substrates"]
        # add uniprot accessions
        df_substrates = df_substrates.merge(
            df.tcdb_id.to_frame().reset_index(drop=False), how="left", on="tcdb_id"
        )
        df_substrates = df_substrates.drop("tcdb_id", axis=1)
        df_substrates = df_substrates[~df_substrates.Uniprot.isnull()]
        df_substrates = df_substrates.set_index("Uniprot", drop=True)
        df = df.join(df_substrates, how="left")

    for keyword_name, keyword_set in KEYWORD_SETS.items():
        df = __add_keyword_column(df, keyword_set=keyword_set, colname=keyword_name)

    if sequence_clustering:
        cluster_repr = cd_hit(
            df.sequence, identity_threshold=sequence_clustering, verbose=verbose
        )
        df = df.loc[cluster_repr]

    return df


if __name__ == "__main__":
    # print entire df (for debugging)
    pd.set_option("expand_frame_repr", False)
    # test input
    from time import time

    t = time()
    df = create_dataset(
        keywords_classes={
            "Amino-acid transport",
            "Sugar transport",
            "Ion transport",
            "Potassium transport",
        },
        keywords_filter={"Transmembrane", "Transport"},
        input_file="data/raw/swissprot/uniprot_data_2022_04.tab.gz",
        multi_substrate="integrate",
        verbose=True,
        tax_ids_filter={3702, 9606, 83333, 559292},
        outliers={
            "Q9HBR0",
            "Q07837",
            "P76773",
            "Q47706",
            "P02943",
            "P75733",
            "P69856",
            "P64550",
            "O81775",
            "Q9SW07",
            "Q9FHH5",
            "Q8S8A0",
            "Q3E965",
            "Q3EAV6",
            "Q3E8L0",
        },
        sequence_clustering=70,
        evidence_code=2,
        invalid_amino_acids="remove_protein",
        force_update=False,
    )
    print(df.shape)
    print(time() - t)
