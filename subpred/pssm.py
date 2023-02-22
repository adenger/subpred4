import pandas as pd
import argparse
import os
from sklearn.preprocessing import minmax_scale
import subprocess
import platform
from .fasta import read_fasta, write_fasta

PSSM_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
PSSM_AA_LIST = list(PSSM_AA_ORDER)
PSSM_AA_SET = set(PSSM_AA_ORDER)

def __process_pssm_file(pssm_file_name, sequence:str):
    # print(pssm_file_name, tmp_folder_name)
    with open(pssm_file_name) as pssm_file:
        next(pssm_file)
        next(pssm_file)

        amino_acids = pssm_file.readline().strip().split()[:20]
        assert amino_acids == PSSM_AA_LIST, f"Unexpexted amino acid order: {amino_acids}"

        # dict keeps insertion order in python 3.7+, therefore is same order as PSSM_AA_ORDER
        amino_acid_to_sum_vector = {
            amino_acid: [0.0] * 20 for amino_acid in amino_acids
        }

        sequence_pssm_file = ""

        for line in pssm_file:
            if line == "\n":  # end of file, before overall scores
                break

            values = line.strip().split()
            amino_acid = values[1]
            assert amino_acid in PSSM_AA_SET, f"unexpected amino acid in pssm file {pssm_file_name}: {amino_acid}"
            sequence_pssm_file += amino_acid

            scores = [float(score) for score in values[2:22]]
            sum_vector = amino_acid_to_sum_vector.get(amino_acid)

            assert len(scores) == 20, f"incomplete PSSM file: {pssm_file_name}. delete and rerun program."

            for pos in range(20):
                sum_vector[pos] += scores[pos]
            amino_acid_to_sum_vector[amino_acid] = sum_vector


        # Can happen for sequence conflicts, like in Q91Y77 position 5
        assert sequence_pssm_file == sequence, f"Sequence from PSSM file {pssm_file_name} did not match input sequence:\n{sequence_pssm_file}\n{sequence}"

        sum_amino_acids = ""
        pssm = []
        for sum_aa, sum_vector in amino_acid_to_sum_vector.items():
            sum_amino_acids += sum_aa
            pssm.extend(sum_vector)

        assert sum_amino_acids == PSSM_AA_ORDER, f"unexpected amino acid in pssm file {pssm_file_name}: {amino_acid}"

        pssm = minmax_scale(pssm).tolist()  # scale to [0,1]

        return pssm


def __create_pssm_file(
    psiblast_location: str,
    fasta_file_name: str,
    pssm_file_name: str,
    blastdb_location: str,
    iterations: int,
    evalue: float,
    threads: int,
) -> None:
    
    # TODO create tmp files, then rename after program finished
    log_file_name = f"{pssm_file_name}.log"
    subprocess.run(
        "{} -query {} -db {} -num_iterations {} -inclusion_ethresh {} -num_threads {} -save_pssm_after_last_round\
             -out_ascii_pssm {} -out {} -comp_based_stats {}".format(
            psiblast_location,
            fasta_file_name,
            blastdb_location,
            iterations,
            evalue,
            threads,
            pssm_file_name,
            log_file_name,
            2
            if iterations == 1
            else 1,  # default is 2, but not supported when matrix is PSSM instead of BLOSUM
        ),
        check=True,
        shell=False if platform.system() == "Windows" else True,
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.DEVNULL,
    )


def __get_pssm_feature(
    accession: str,
    sequence: str,
    blastdb_fasta_file: str,
    pssm_folder_path: str,
    psiblast_location: str,
    iterations: int = 1,
    evalue: float = 0.002,
    threads: int = 1,
    verbose: bool = False,
) -> list:

    if not os.path.exists(pssm_folder_path):
        os.makedirs(pssm_folder_path)

    fasta_file_name = f"{pssm_folder_path}/{accession}.fasta"
    pssm_file_name = f"{pssm_folder_path}/{accession}.pssm"

    pssm = []
    if os.path.isfile(pssm_file_name):
        pssm = __process_pssm_file(pssm_file_name, sequence)
        # if verbose:
            # print(f"PSSM for accession {accession} was found in tmp folder {pssm_folder_path}")
    else:
        if verbose:
            print(
                f"PSSM for accession {accession} was not found in tmp folder {pssm_folder_path}, calling psiblast"
            )
        write_fasta(
            fasta_file_name=fasta_file_name, fasta_data=[(">" + accession, sequence)]
        )
        __create_pssm_file(
            psiblast_location=psiblast_location,
            fasta_file_name=fasta_file_name,
            pssm_file_name=pssm_file_name,
            blastdb_location=blastdb_fasta_file,
            iterations=iterations,
            evalue=evalue,
            threads=threads,
        )
        pssm = __process_pssm_file(pssm_file_name, sequence)
        if verbose:
            print(f"PSSM for accession {accession} was generated")

    return pssm


def calculate_pssm_feature(
    sequences: pd.Series,
    tmp_folder: str,
    blast_db: str,
    iterations: int,
    psiblast_executable: str = "psiblast",
    psiblast_threads: int = 4,
    verbose: bool = False,
):
    accessions = list()
    features = list()
    errors = list()

    for i in range(len(sequences)):
        try:
            pssm = __get_pssm_feature(
                accession=sequences.index[i],
                sequence=sequences.values[i],
                blastdb_fasta_file=blast_db,
                pssm_folder_path=tmp_folder,
                iterations=iterations,
                psiblast_location=psiblast_executable,
                threads=psiblast_threads,
                verbose=verbose,
            )
            accessions.append(sequences.index[i])
            features.append(pssm)
        except StopIteration:
            errors.append(f"Error: Stopiteration occurred for {f'{tmp_folder}/{sequences.index[i]}.pssm'}. File might be empty")
            continue
        except AssertionError as e:
            errors.append(f"Error: AssertionError occurred for {f'{tmp_folder}/{sequences.index[i]}.pssm'}. Message:{e}")
            continue
    
    for error in errors: 
        print(error)

    # pssm_aa_order = "ARNDCQEGHILKMFPSTWYV"
    pssm_aa_substitutions = [
        aa1 + aa2 for aa1 in PSSM_AA_ORDER for aa2 in PSSM_AA_ORDER
    ]

    pssm_df = pd.DataFrame(
        data=features, index=accessions, columns=pssm_aa_substitutions
    )
    return pssm_df


# special hardcoded function for the notebooks
# def calculate_pssms_notebook(
#     sequences: pd.Series,
#     uniref_cluster_thresholds: list = [50, 90],
#     psiblast_iterations_list: list = [1, 3],
#     n_threads: int = 4,
# ):
#     df_pssm_all = pd.DataFrame()
#     for uniref_cluster_threshold in uniref_cluster_thresholds:
#         for psiblast_iterations in psiblast_iterations_list:
#             df_pssm = calculate_pssm_feature(
#                 sequences,
#                 tmp_folder="../data/intermediate/blast/pssm_uniref{}_{}it".format(
#                     uniref_cluster_threshold, psiblast_iterations
#                 ),
#                 blast_db="../data/raw/uniref/uniref{}/uniref{}.fasta".format(
#                     uniref_cluster_threshold, uniref_cluster_threshold
#                 ),
#                 iterations=psiblast_iterations,
#                 psiblast_executable="psiblast",
#                 psiblast_threads=n_threads,
#                 verbose=False,
#             )
#             df_pssm = df_pssm.rename(
#                 columns=lambda c: f"PSSM_{uniref_cluster_threshold}_{psiblast_iterations}__{c}"
#             )

#             df_pssm_all = pd.concat([df_pssm_all, df_pssm], axis=1)

#     return df_pssm_all
