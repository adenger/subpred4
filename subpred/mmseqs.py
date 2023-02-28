import tempfile
import subprocess
import os
from subpred.fasta import read_fasta, write_fasta
import pandas as pd


def __flatten_kwargs(**kwargs):
    kwargs_list = list()
    for k, v in kwargs.items():
        kwargs_list.append(k)
        kwargs_list.append(v)
    return kwargs_list


def mmseq(
    sequences: pd.Series,
    fast_algo: bool = False,
    cluster_mode: int = 0,
    min_seq_id: float = 0.0,
    min_coverage: float = 0.8,
    cov_mode: int = 0,
    max_evalue: float = 1e-3,
    verbose: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """Python wrapper for MMSeqs2

    Args:
        sequences (pd.Series):
            Series with identifiers as index and amino acid sequences as values
        fast_algo (bool, optional):
            Use linclust algorithm (for huge datasets). Defaults to False.
        cluster_mode (int, optional):
            0: Set-Cover (greedy)
            1: Connected component (BLASTclust)
            2,3: Greedy clustering by sequence length (CDHIT).
            Defaults to 0.
        min_seq_id (float, optional):
            List matches above this sequence identity (for clustering) (range 0.0-1.0). Defaults to 0.0.
        min_coverage (float, optional):
            -c. List matches above this fraction of aligned (covered) residues (see cov-mode). Defaults to 0.8.
        cov_mode (int, optional):
            0: coverage of query and target
            1: coverage of target
            2: coverage of query
            3-5 also exist
            Defaults to 0.
        max_evalue (float, optional):
            List matches below this E-value (range 0.0-inf). Defaults to 1e-3.
        verbose (bool, optional):
            Print command, and output of MMSeq2. Defaults to False.
    """

    with (
        tempfile.NamedTemporaryFile(suffix=".fasta") as tmp_fasta_in,
        tempfile.NamedTemporaryFile() as tmp_out,
        tempfile.TemporaryDirectory() as tmp_folder,
    ):
        fasta_data = list(
            zip(
                [">" + ac for ac in sequences.index.tolist()], sequences.values.tolist()
            )
        )
        write_fasta(fasta_file_name=tmp_fasta_in.name, fasta_data=fasta_data)

        seq_db_path = tmp_fasta_in.name.replace(".fasta", ".seqdb")

        cluster_fasta = tmp_out.name + "_all_seqs.fasta"
        cluster_tsv = tmp_out.name + "_cluster.tsv"
        cluster_representatives_fasta = tmp_out.name + "_rep_seq.fasta"

        execution = [
            "mmseqs",
            "easy-linclust" if fast_algo else "easy-cluster",
            # "easy-linclust",
            tmp_fasta_in.name,
            tmp_out.name,
            tmp_folder,
            "--cluster-mode",
            cluster_mode,
            "--min-seq-id",
            min_seq_id,
            "-c",
            min_coverage,
            "--cov-mode",
            cov_mode,
            "-e",
            max_evalue,
        ] + __flatten_kwargs(**kwargs)
        execution = [
            str(argument) if not isinstance(argument, str) else argument
            for argument in execution
        ]
        result = subprocess.run(
            execution, check=True, stdout=subprocess.PIPE, universal_newlines=True
        )

        if verbose:
            print("==COMMAND==")
            print(" ".join(execution))
            print("==STDOUT==")
            print(result.stdout)
            print("==STDERR==")
            print(result.stderr)

        cluster_tsv = pd.read_table(cluster_tsv, header=None, names=["Rep", "Cluster"])

    return cluster_tsv
