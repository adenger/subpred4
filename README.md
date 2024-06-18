# Subpred4

## Setup (tested on Ubuntu 22.04 LTS)

1. Clone repository

2. Download subpred4_data.tar.gz and place in repository folder

    [OneDrive download link (~50GB)](https://unisaarlandde-my.sharepoint.com/:u:/g/personal/ande010_uni-saarland_de/EdtikTFsnuJGoUhtmvnM1PkBXGHGBB15ipbmWZco3ZrQag?e=6kcVpd)

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

5. Recreate exact conda environment (use environment_history.yml instead if there is an error, which can happen on different OS)

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

9. Run the notebooks in order, according to their filenames

## API concept

All raw data is left untouched in data/raw. The download commands and versions can be found in the preprocessing notebook. All files are based on the same version of Uniprot (2022_05). Re-downloading the raw data using the same commands can upgrade them to the latest version, but that can lead to incompatibilities, since not all databases based on a particular Uniprot version are released at the same time.

Preprocessing is performed on the raw data, then the processed data is saved as pickles in data/datasets for fast i/o. The method subpred.util.load_df can be used to read these pickles.

A transporter dataset can be created manually with all parameters using the methods in *subpred.protein_dataset*, *subpred.go_annotations* and *subpred.chebi_annotations*. This process is simplified through the function *subpred.transmembrane_transporters.get_transmembrane_transporter_dataset*, which sets most of the parameters.

The function *get_transmembrane_transporter_dataset* returns three dataframes: One with sequences, one with GO annotations, and one with ChEBI annotations. These three dataframes essentially act like data classes. All of the remaining methods in the package take one or multiple of these dataframes as input to carry out their calculations, and the data should ideally not be changed before using the methods on them.
