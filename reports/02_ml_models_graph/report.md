# Training ML models on Chebi graph

## SVM models, comparisons

TODO implement pipeline(s)

## TODOs (monday here)

- test cases for organisms and substrates from chebi graph picture
  - Include preprocessing in flowchart
  - why are there non-transporters in the dataset? also: pipeline not deterministic
    - Quickfix: just filter them again in the notebook, see how many are left
      - class for creating go-slims would make things easier
    - Actual fix: Create flowchart of dataset creation to check how hard it is to redo the dataset creation pipeline in a simpler way.
- extract scores for train and test, put them into table
  - together with sample counts, dataset name, substrates and other measures (see below)
  - create protein/go/chebi dataset, with ability to subset by substrate combination and organism
- For adjacency matrix:
  - If substrates are 1 in a column (or row?) then they have same parent
  - Automation for all pairs in the matrix with enough samples, not just the test cases
- Then: calculate correlations, statistical tests, etc.
  - For the matrices

## Data matrices for test cases

- For each (organism, substrate1, substrate2) multiple matrices/heatmaps:
  - is_a connections (nx.from_pandas_adjacency, nx.adjacency_matrix)
  - SVM scores for one model (select one)
    - F1 train score for each substrate (or average)
    - F1 test score for each substrate (or average)
  - Overlap matrix
  - average/max sequence similarity
  - blast/mmseqs2 performance im Vergleich mit SVM
  - protein embeddings cosine similarity
  - chemical similarity (smiles etc.)

- For each (organism, go_term1, go_term2) multiple matrices/heatmaps:
  - is_a connections (nx.from_pandas_adjacency, nx.adjacency_matrix)
  - SVM scores for one model (select one)
    - F1 train score for each substrate (or average)
    - F1 test score for each substrate (or average)
  - Overlap matrix
  - average/max sequence similarity
  - blast/mmseqs2 performance im Vergleich mit SVM
  - protein embeddings cosine similarity
  - semantic similarity

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

## Done

I created a data flowchart to get an overview of the processing and went through all the code. Adding ancestor chebi terms did not make difference.

I implemented a SVM-RBF model with basic hyperparameter optimization. There are three kinds of features selection that are compared:

- Anova-based feature selection that optimized the percentile of best features
- FeatureCombinator from manuscript 1 tries all combinations of features and feature generation parameters (such as psiblast iterations or blast database) and selects the optimal one based on training data.
- FeatureCombinator followed by Anova

Precision-Recall curves are created for each test case and each model. ROC curves have problems with imbalanced data. Also there is the F1 score for the training dataset, and a classification_report for the testing dataset. ROC curves showed misleading results. This could be due to the two separate evaluations, or because of PR curves in general. Link with explanation: <https://medium.com/@douglaspsteen/precision-recall-curves-d32e5b290248>

The model is optimized using F1 score, and evaluated with F1/precision/recall for the individual classes.

## Comparison of SVM results with other measures


## One organism vs different meta organisms

## Automate process of property annotation

BFS

## Transfer learning from E coli to similar bacteria, compare to BLAST
