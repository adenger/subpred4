from joblib.parallel import Parallel, delayed
from collections import Counter
import numpy as np

AMINO_ACIDS = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]
DIPEPTIDES = [a + b for a in AMINO_ACIDS for b in AMINO_ACIDS]
DIPEPTIDES_DICT = {dp: 0 for dp in DIPEPTIDES}


def paac1(seq):
    return [seq.count(aa) / (len(seq) - 1) for aa in DIPEPTIDES]


def paac2(seq):
    counter = Counter(DIPEPTIDES_DICT)
    counter.update([seq[i : i + 2] for i in range(len(seq) - 1)])
    return [x / (len(seq) - 1) for x in counter.values()]


def paac3(seq):
    d = DIPEPTIDES_DICT
    for i in range(len(seq) - 1):
        peptide = seq[i : i + 2]
        d[peptide] += 1
    return [v / (len(seq) - 1) for v in d.values()]


def paac_batch1(seqs):
    ret = []
    counter = Counter()
    for seq in seqs:
        counter.update(DIPEPTIDES_DICT)
        counter.update([seq[i : i + 2] for i in range(len(seq) - 1)])
        ret.append([x / (len(seq) - 1) for _, x in counter.items()])
        counter.clear()
    return ret


def paac_batch2(seqs):
    return [[seq.count(aa) / (len(seq) - 1) for aa in DIPEPTIDES] for seq in seqs]


def unpack(list3d):
    for paac_lists in list3d:
        for paac_list in paac_lists:
            yield paac_list


def paac(sequences, protein_names=None, n_chunks=80):

    df_uniprot = df_uniprot[["sequence"]]

    # sequences = df_uniprot.sequence.to_numpy()[:500302]
    if not isinstance(sequences, np.ndarray):
        sequences = np.array(sequences)

    sequences_chunks = np.array_split(sequences, n_chunks)
    paac_function = paac_batch1

    results = Parallel(n_jobs=-1)(
        delayed(paac_function)(chunk) for chunk in sequences_chunks
    )
    results = list(unpack(results))

    # TODO dataframe

    # x = Parallel(n_jobs=-1)(delayed(paac)(seq) for seq in sequences)
    # x = Parallel(n_jobs=-1)(delayed(paac2)(seq) for seq in sequences)
    # x = Parallel(n_jobs=-1)(delayed(paac3)(seq) for seq in sequences)


# TODO actually get dipeptides
