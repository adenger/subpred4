import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit import DataStructs
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


def get_fingerprint(smiles_str, radius: int = 2, bits: int = 1024, draw: bool = False):
    # radius determines the number of bonds the algo with look around an atom to find substructures
    bit_info = {}
    mol = Chem.MolFromSmiles(smiles_str)
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=radius, nBits=bits, bitInfo=bit_info
    )
    if draw:
        return Draw.DrawMorganBits(
            [(mol, on_bit, bit_info) for on_bit in fingerprint.GetOnBits()],
            molsPerRow=4,
            legends=[str(on_bit) for on_bit in fingerprint.GetOnBits()],
        )
    return fingerprint


def tanimoto(set1, set2):
    tanimoto = len(set1 & set2) / len(set1 | set2)
    return tanimoto


def get_pairwise_tanimoto(chebi_ids: list):
    graph_chebi = load_df("chebi_obo")
    chebi_smiles_dict = get_chebi_smiles_dict(graph_chebi=graph_chebi)
    smiles = [chebi_smiles_dict.get(chebi_id) for chebi_id in chebi_ids]
    df_smiles = (
        pd.DataFrame({"chebi_id": chebi_ids, "smiles": smiles})
        .drop_duplicates()
        .set_index("chebi_id")
    )
    df_smiles = df_smiles[~df_smiles.smiles.isnull()]

    df_smiles["fingerprint"] = df_smiles.smiles.map(get_fingerprint)

    df_chem_similarity = pd.DataFrame(
        [
            [DataStructs.TanimotoSimilarity(fingerprint1, finterprint2) for fingerprint1 in df_smiles.fingerprint]
            for finterprint2 in df_smiles.fingerprint
        ],
        index=df_smiles.index,
        columns=df_smiles.index,
    )
    return df_chem_similarity