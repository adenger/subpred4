The Makefiles included here download the most recent versions of the raw datasets.

TODO readme text for manuscript 2

## Setup:

1. Install Mambaforge
```
https://github.com/conda-forge/miniforge#mambaforge
```
2. Recreate conda environment:
```
mamba env create --file environment.yml
```
3. Activate conda environment: 
```
conda activate subpred4
```
4. Install code as python package into environment: 
```
pip install -e .
```
5. Download raw data: 
```
make raw_data
```
6. Create BLAST databases (Needs >100GB of space and takes several hours): 
```
make blast_databases
```

<!-- ## Reproduce results from manuscript:

 1. Install miniconda
2. Recreate conda environment:
```
conda env create --file environment.yml
```
3. Activate conda environment: 
```
conda activate subpred
```
4. Install code as python package into environment: 
```
pip install -e .
```
5. Download data_full.tar from https://cloud.hiz-saarland.de/s/sGTyGApAqdgAQiB and place it in repository
6. Rename existing data folder:
```
mv data data_bak
```
7. Extract tar archive:
```
make raw_data_manuscript
```
8. Create BLAST databases (Needs >100GB of space and takes several hours):
    - This step is optional, as the previous step extracts pre-computed PSSMs for all proteins to *data/intermediate/blast*
  
```
make blast_databases
``` -->
