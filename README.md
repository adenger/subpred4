The Makefiles included here download the most recent versions of the raw datasets.

TODO readme text for manuscript 2

## Setup:

1. Clone repository
2. Download data_backup.tar.gz and place in repository folder
3. Extract raw data:
```
make data_import
```
4. Install Mambaforge
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
source ~/.bashrc
```
5. Recreate exact conda environment (use the other yml file if there is an error):
```
mamba env create --file environment.yml
```
6. Activate conda environment: 
```
conda activate subpred4
```
7. Install code as python package into environment: 
```
make package
```
8. Create BLAST databases (Needs >100GB of space and takes several hours, pre-computed pssms are availabe in data/intermediate):
```
make blast_databases
```
