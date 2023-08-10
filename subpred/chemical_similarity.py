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
        columns=df_smiles.index,
    )
    return df_chem_similarity
