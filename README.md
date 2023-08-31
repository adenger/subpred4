# Subpred4

## Setup

1. Clone repository

2. Download subpred4_data.tar.gz and place in repository folder

3. Extract raw data

    ```bash
    make data_import
    ```

4. Install Mambaforge

    ```bash
    wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
    bash Mambaforge-$(uname)-$(uname -m).sh
    source ~/.bashrc
    ```

5. Recreate exact conda environment (use the other yml file if there is an error)

    ```bash
    mamba env create --file environment.yml
    ```

6. Activate conda environment

    ```bash
    conda activate subpred4
    ```

7. Install code as python package into environment

    ```bash
    make package
    ```

8. Create BLAST databases (Needs >100GB of space and takes several hours, pre-computed pssms are availabe in data/intermediate)

    ```bash
    make blast_databases
    ```

## Data concept

All raw data is left untouched in data/raw. The download commands can be found in the Makefile, and in the preprocessing notebook. We made sure that all files are based on the same version of Uniprot (2022_05). Re-downloading the raw data will upgrade then to the latest version, although that can lead to incompatibilities, since not all databases based on a particular Uniprot version are released at the same time.

Preprocessing if performed on the raw data, then the processed data is saved as pickles in data/datasets for fast i/o. The method subpred.util.load_df can be used to read these pickles.

A transporter dataset can be created manually with all parameters using the methods in *subpred.protein_dataset*, *subpred.go_annotations* and *subpred.chebi_annotations*. This process is simplified by the function *subpred.transmembrane_transporters.get_transmembrane_transporter_dataset*, which sets most of the parameters for you.

The function *get_transmembrane_transporter_dataset* returns three dataframes: One with sequences, one with GO annotations, and one with ChEBI annotations. All of the remaining methods in the package take one or multiple of these dataframes as input to carry out their calculations, since they contain all the necessary data.
