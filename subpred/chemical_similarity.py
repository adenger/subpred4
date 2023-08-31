import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# from rdkit.Chem import AllChem
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.rdMolDescriptors import (
    GetMorganFingerprintAsBitVect,
    GetHashedAtomPairFingerprintAsBitVect,
    GetHashedTopologicalTorsionFingerprintAsBitVect,
    GetMACCSKeysFingerprint,
)
from subpred.util import load_df

# from rdkit.Chem.Draw import IPythonConsole
# from rdkit.Chem import Descriptors


def get_chebi_smiles_dict(graph_chebi):
    chebi_to_smiles = dict()

    pattern_smiles = re.compile(
        '^http://purl.obolibrary.org/obo/chebi/smiles "(.*?)" xsd:string$'
    )

    for chebi_id, data in graph_chebi.nodes(data="property_value"):
        if not data:
            continue
        for data_point in data:
            if match_obj := re.search(pattern_smiles, data_point):
                chebi_to_smiles[chebi_id] = match_obj.group(1)
                break

    return chebi_to_smiles


# def draw_morgan(
#     smiles_str: str,
#     radius: int = 2,
#     bits: int = 2048,
# ):
#     bitinfo = dict()
#     mol = Chem.MolFromSmiles(smiles_str)
#     fingerprint = GetMorganFingerprintAsBitVect(
#         mol, radius=radius, nBits=bits, bitInfo=bitinfo
#     )
#     return Draw.DrawMorganBits(
#         [(mol, on_bit, bitinfo) for on_bit in fingerprint.GetOnBits()],
#         molsPerRow=4,
#         legends=[str(on_bit) for on_bit in fingerprint.GetOnBits()],
#     )


def get_fingerprint(
    smiles_str: str,
    method="morgan",
    radius: int = 2,
    n_bits: int = 2048,
):
    # radius determines the number of bonds the algo with look around an atom to find substructures
    mol = Chem.MolFromSmiles(smiles_str)

    # ideas: dice similarity, int vectors, hashed int vectors, compare?

    # https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html
    # https://scikit-chem.readthedocs.io/en/stable/_modules/skchem/descriptors/fingerprints.html
    # https://www.rdkit.org/UGM/2012/Landrum_RDKit_UGM.Fingerprints.Final.pptx.pdf
    match (method):
        case "morgan":
            fingerprint = GetMorganFingerprintAsBitVect(
                mol, radius=radius, nBits=n_bits
            )
        case "atompairs":
            fingerprint = GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits)
        case "torsions":
            fingerprint = GetHashedTopologicalTorsionFingerprintAsBitVect(
                mol, nBits=n_bits
            )
        case "maccs":
            fingerprint = GetMACCSKeysFingerprint(mol)
        # case "morgan_hash":
        #     fingerprint = GetHashedMorganFingerprint(
        #         mol, radius=radius, nBits=n_bits
        #     )
        # case "atompairs_hash":
        #     fingerprint = GetHashedAtomPairFingerprint(mol, nBits=n_bits)
        # case "torsions_hash":
        #     fingerprint = GetHashedTopologicalTorsionFingerprint(
        #         mol, nBits=n_bits
        #     )
        case _:
            raise ValueError(f"invalid fingerprint method: {method}")

    return fingerprint


# def tanimoto(set1, set2):
#     tanimoto = len(set1 & set2) / len(set1 | set2)
#     return tanimoto


def get_pairwise_similarity(chebi_ids: list, fingerprint_method: str = "morgan"):
    """Calculate pairwise tanimoto similarities

    Args:
        chebi_ids (list): list of chebi identifiers
        fingerprint_method (str, optional):
            Options: "morgan", "atompairs", "torsions", "maccs".
            Defaults to "morgan".

    Returns:
        pd.DataFrame: DataFrame with pairwise tanimoto scores.
            Only keeps ids for which a smiles string can be found
            and for which a fingerprint can be calculated
    """
    graph_chebi = load_df("chebi_obo")
    chebi_smiles_dict = get_chebi_smiles_dict(graph_chebi=graph_chebi)
    smiles = [chebi_smiles_dict.get(chebi_id) for chebi_id in chebi_ids]
    df_smiles = (
        pd.DataFrame({"chebi_id": chebi_ids, "smiles": smiles})
        .drop_duplicates()
        .set_index("chebi_id")
    )
    df_smiles = df_smiles[~df_smiles.smiles.isnull()]

    df_smiles["fingerprint"] = df_smiles.smiles.apply(
        get_fingerprint, method=fingerprint_method
    )
    # FingerprintSimilarity seems to give same results as tanimoto. throws error for int vectors
    df_chem_similarity = pd.DataFrame(
        [
            [
                FingerprintSimilarity(fingerprint1, finterprint2)
                for fingerprint1 in df_smiles.fingerprint
            ]
            for finterprint2 in df_smiles.fingerprint
        ],
        index=df_smiles.index,
        columns=df_smiles.index.rename("chebi_id2"),
    )
    return df_chem_similarity


def tanimoto_chebi_to_go(
    df_tanimoto_chebi: pd.DataFrame,
    df_go_chebi: pd.DataFrame,
    agg_function="mean",
    primary_input_only: bool = True,
) -> pd.DataFrame:
    """Converts tanimoto scores between molecules into aggregated tanimoto scores between go terms annotated with those molecules.

    Args:
        df_tanimoto_chebi (pd.DataFrame): The dataframe to convert. Created with get_pairwise_similarity
        df_go_chebi (pd.DataFrame): GO-Chebi mapping df from subpred.transmembrane_transporters
        agg_function (str, optional): aggregation function, if two GO terms have multiple tanimoto scores.
            Can be any aggr. function, for example min, max, median, mean, etc. Defaults to "mean".
        primary_input_only (bool, optional): Whether to only look at transported substrates (True),
            Or also at interacting molecules such as ATP, H2O, etc. Defaults to True.

    Returns:
        pd.DataFrame: Aggretated Tanimoto scores between go terms, based on the substrates of their proteins.
    """
    df_go_chebi_local = (
        df_go_chebi[df_go_chebi.chebi_go_relation == "has_primary_input"].copy()
        if primary_input_only
        else df_go_chebi.copy()
    )

    df_tanimoto_go = df_tanimoto_chebi.unstack().reset_index(name="tanimoto")
    chebi_id_to_go_ids = (
        df_go_chebi_local.groupby("chebi_id")
        .apply(lambda x: x.go_id.sort_values().unique().tolist())
        .to_dict()
    )
    df_tanimoto_go["go_id1"] = df_tanimoto_go.chebi_id.map(chebi_id_to_go_ids)
    df_tanimoto_go["go_id2"] = df_tanimoto_go.chebi_id2.map(chebi_id_to_go_ids)
    df_tanimoto_go = (
        df_tanimoto_go.explode("go_id1").explode("go_id2").reset_index(drop=True)
    )

    df_tanimoto_go = (
        df_tanimoto_go.drop(["chebi_id", "chebi_id2"], axis=1)
        .groupby(["go_id1", "go_id2"], as_index=False)
        .agg(agg_function)
    )
    return df_tanimoto_go.pivot(index="go_id1", columns="go_id2", values="tanimoto")
