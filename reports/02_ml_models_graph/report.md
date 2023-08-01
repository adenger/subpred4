# Training ML models on Chebi graph

## SVM models, comparisons

TODO implement pipeline(s)

## TODOs

- test cases for organisms and substrates from chebi graph picture
  - wait for PSSMs from tera to finish
  - why are there non-transporters in the dataset?
    - Quickfix: just filter them again in the notebook, see how many are left
      - class for creating go-slims would make things easier
    - Actual fix: Create flowchart of dataset creation to check how hard it is to redo the dataset creation pipeline in a simpler way.
- extract eval scores for train and test, put them into table
  - together with sample counts, dataset name, substrates
- pipeline is not deterministic
  - only changes sometimes
  - tested: saving sequence dataset and comparing to new one, no difference. 

## Test cases

- athaliana
  - Ca2+ und K+
  - inorganic anion/cation
  - carboxylic acid anion/inorganic anion
- ecoli
  - carbohydrate derivate / monosaccharide
  - minicarboxylic acid / amino acid
- human
  - Ca2+ / Na1+
  - Ca2+ / K+
  - Na+ / K+
  - inorganic anion/cation
- hefe
  - amid / amino acid derivative

## ML-Model

I implemented a simple SVM-RBF model with basic hyperparameter optimization. There are two kinds of (optional) features selection:

- Anova-based feature selection that optimized the percentile of best features
- FeatureCombinator from manuscript 1 tries all combinations of features and feature generation parameters (such as psiblast iterations or blast database) and selects the optimal one based on training data.

The Chebi-terms are used as labels. 

The model is optimized using F1 score, and evaluated with F1/precision/recall for the individual classes.

## Comparison of SVM results with other measures

je mehr trennung desto weniger accuracy?
Macht auftrennung chemisch sinn?
go semanitc smilarity, maximale sequenzidentit, overlap

## One organism vs different meta organisms

## Automate process of property annotation

BFS

## Transfer learning from E coli to similar bacteria, compare to BLAST
