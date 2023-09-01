# What determines a good classification performance of a substrate class?


## Overview

On one hand, we have a sub-graph of the Gene Ontology, which contains all GO terms that are descendants of *transmembrane transporter activity*, a molecular function term. On the other hand, we have a set of proteins that are annotated with at least one term in that sub-graph. The task is to predict the transport function term or terms with which a particular protein is annotated.

There are two issues: There are too many terms for multiclass classification, and there are too few samples per term. There are approaches that train a classifier for each term in GO, but they don't achieve accuracies above 0.65.

This is called *hierarchical classification*, and there are multiple apporaches:

- Multi-Class classification between all leafs. There are many leafs, and they often have very few samples.
- Find a set of abstract terms with little overlap in terms of descendants and annotated proteins, and put all of their respective descendant terms into a single class
- Train a classifier for each junction in the tree. Since we have to multiply the accuracies, this will automatically become less accurate the more we move down the tree.

## Dataset creation

The dataset creation is abstracted to a point where it only needs a set of organism identifiers, a flag that includes or excludes proteins that have not been manually reviewed (TrEMBL), a flag that can exclude electronically inferred GO terms (useful for organisms that are not model organisms), and the evidence code for the existence of the sequence (proven at transcript level, or proven at protein level).

After filtering the proteins from UniprotKB, a subset of membrane transporter substrate predictions is created:

- The root term is *transmembrane transporter activity*
- We only keep edges that are annotated with *is_a* (direct logical relation)
- Only keep *molecular_function* terms
- Only terms that annotate proteins from our filtered protein set
- Only terms that have the *enables* relation to their protein, which means that the function is directly carried out by the protein. The *part_of* relation only added outliers, like parts of transport protein complexes that don't carry out the transport directly.
- Finally, the each term is annotated with its ancestors, all the way up to the root node, but only using *is_a* edges. This increases the number of samples drastically. For example, a protein annotated with *sugar transmembrane transporter activity* is also annotated with *carbohydrate transmembrane transporter activity*, because there is an *is_a* relation between the two terms.

Lastly, the chebi terms for the GO terms are retrieved from a dataset of cross-ontology relations. ChEBI is a database of biologically relevant molecules, and also an ontology. These can be filtered by 3-star terms, which are manually curated. For our yeast dataset, this did not make a difference. There are two types of cross-ontology relations that are relevant: *has_primary_input* typically denotes the transported substrate or substrates, and *has_participant* typically includes other molecules that contribute to the transport, such as ATP or H2O.

In contrast to GO, ChEBI is not acyclid. This means that adding ancestors is not as trivial as with GO, and traversing the graph to find ancestors can often add thousands of terms that have some path to the term we are looking at. Therefore, we decided against using ChEBI terms and the ChEBI ontology directly, and instead see them as annotations for our GO graph.

Another point is that not every GO term has its appropriate ChEBI term in the database. For example, there is a ChEBI term for *Ion*, but the *monoatomic ion transmembrane transport* GO term is not annotated with that.

## Stats

For our test notebook, we created a dataset of *S. cerevisiae* transporters, with the following parameters:

```python
    ORGANISM_IDS={559292}
    SWISSPROT_ONLY=False
    MAX_SEQUENCE_EVIDENCE_CODE = 1  # Only evidence at protein level
    EXCLUDE_IEA_GO_TERMS=False
```

The yeast dataset contains:

- 322 unique proteins before filtering out similar sequences with cd-hit.
- 288 unique GO terms related to *transmembrane transporter activity*, after adding ancestors
- 211 unique GO terms related to *transmembrane transporter activity*, *before* adding ancestors
- 226 GO term are annotated with ChEBI terms
- 181 unique ChEBI terms (some of them abstract, not directly usable)
- 89 unique ChEBI terms that are usable for Tanimoto coefficient (which usually means they are concrete molecules, not abstract terms)
- 131 unique GO terms annotated to the 89 ChEBI terms that are usable for Tanimoto coefficient
- 36 GO terms with more than 20 samples per class and at least 20 unique proteins per class (i.e. usable for ML). We could try going down to at least 15 unique proteins per class.


## Matrices

### GO pairwise score matrices

- Binary adjacency matrix of GO terms, only using *is_a* relations. **$288 \times 288$**
- Overlap matrix of GO terms, i.e. how many proteins two terms have in common. **$288 \times 288$**
- Semantic similarity matrices. **$288 \times 288$**
  - Wang similarity algorithm: Fast, parallel
  - TODO: other five algorithms (only one core, takes very long)
- Sequence scores. **$288 \times 288$**, 8 matrices
  - Sequence identity between proteins annotated with GO terms. Identity function: (identical positions) / (aligned positions + internal gap positions)
  - Needleman-Wunsch alignment score between sequences
  - Aggr. functions: mean, median, max, min
- Tanimoto coefficients between ChEBI terms related to GO terms. **$131 \times 131$** 16 matrices
  - Fingerprint algorithms: Morgan, atom pairs, topological torsions, MACCS
  - Aggr. functions: mean, median, max, min
  - Problem: GO terms with enough proteins for ML are too abstract to have concrete ChEBI terms, only two molecules available for a total of 12/36 terms

### ChEBI matrices

- ChEBI adjacency matrix of ChEBI terms directly annotated to ancestor GO terms **$159 \times 159$**
- ChEBI overlap matrix (of proteins related to go terms related to ChEBI terms) **$159 \times 159$**
- Tanimoto coefficients **$89 \times 89$**, four matrices

### ML results matrices

- GO terms Machine learning F1 score matrices **$36 \times 36$**
  - SVM pipeline with no feature selection, all 1600 features
    - training scores
    - test scores
    - Matrices are asymetrical: score (i,j) is classification if GO term i is positive class, and GO term j is negative class. The average of i,j and j,i is the macro-averaged F1 score.
  - SVM pipeline with feature selection TODO
  - SVM pipeline with PCA TODO

## Analysis

TODO method that exports pandas to markdown file, then include file here

### Stats

|       |   mean_train_score |   mean_test_score |   overlap |   semantic_sim_wang |
|:------|-------------------:|------------------:|----------:|--------------------:|
| count |             404    |            404    |    404    |              404    |
| mean  |               0.86 |              0.87 |      4.86 |                0.37 |
| std   |               0.09 |              0.09 |     11.82 |                0.13 |
| min   |               0.52 |              0.54 |      0    |                0.13 |
| 25%   |               0.82 |              0.83 |      0    |                0.27 |
| 50%   |               0.87 |              0.89 |      0    |                0.35 |
| 75%   |               0.93 |              0.94 |      3.25 |                0.45 |
| max   |               0.98 |              0.99 |     89    |                0.74 |

### Pearson correlation between matrices and train/test score

|                                    |   mean_train_score |   mean_test_score |
|:-----------------------------------|-------------------:|------------------:|
| mean_train_score                   |          1         |         0.971266  |
| mean_test_score                    |          0.971266  |         1         |
| overlap                            |         -0.0288888 |        -0.0469669 |
| semantic_sim_wang                  |         -0.281933  |        -0.257315  |
| go_median_sequence_identity        |         -0.0563852 |        -0.0619873 |
| go_median_sequence_alignment_score |         -0.245457  |        -0.250395  |
| go_mean_sequence_identity          |          0.0426981 |         0.022596  |
| go_mean_sequence_alignment_score   |         -0.148291  |        -0.156709  |
| go_max_sequence_identity           |         -0.317217  |        -0.285438  |
| go_max_sequence_alignment_score    |         -0.214941  |        -0.204546  |
| go_min_sequence_identity           |          0.304803  |         0.275133  |
| go_min_sequence_alignment_score    |          0.450618  |         0.425358  |
| tanimoto_morgan_go_mean            |        nan         |       nan         |
| tanimoto_atompairs_go_mean         |        nan         |       nan         |
| tanimoto_torsions_go_mean          |          0.160065  |         0.277435  |
| tanimoto_maccs_go_mean             |          0.160065  |         0.277435  |
| tanimoto_morgan_go_median          |        nan         |       nan         |
| tanimoto_atompairs_go_median       |        nan         |       nan         |
| tanimoto_torsions_go_median        |          0.160065  |         0.277435  |
| tanimoto_maccs_go_median           |          0.160065  |         0.277435  |
| tanimoto_morgan_go_min             |        nan         |       nan         |
| tanimoto_atompairs_go_min          |        nan         |       nan         |
| tanimoto_torsions_go_min           |          0.160065  |         0.277435  |
| tanimoto_maccs_go_min              |          0.160065  |         0.277435  |
| tanimoto_morgan_go_max             |        nan         |       nan         |
| tanimoto_atompairs_go_max          |        nan         |       nan         |
| tanimoto_torsions_go_max           |          0.160065  |         0.277435  |
| tanimoto_maccs_go_max              |          0.160065  |         0.277435  |

### Spearman correlation between matrices and train/test score

|                                    |   mean_train_score |   mean_test_score |
|:-----------------------------------|-------------------:|------------------:|
| mean_train_score                   |          1         |         0.9683    |
| mean_test_score                    |          0.9683    |         1         |
| overlap                            |         -0.390259  |        -0.387189  |
| semantic_sim_wang                  |         -0.341199  |        -0.317716  |
| go_median_sequence_identity        |         -0.0515734 |        -0.0645598 |
| go_median_sequence_alignment_score |         -0.11866   |        -0.120075  |
| go_mean_sequence_identity          |         -0.0232387 |        -0.0385864 |
| go_mean_sequence_alignment_score   |         -0.0239223 |        -0.035882  |
| go_max_sequence_identity           |         -0.550163  |        -0.520922  |
| go_max_sequence_alignment_score    |         -0.535236  |        -0.533938  |
| go_min_sequence_identity           |          0.398639  |         0.378385  |
| go_min_sequence_alignment_score    |          0.503448  |         0.477742  |
| tanimoto_morgan_go_mean            |        nan         |       nan         |
| tanimoto_atompairs_go_mean         |        nan         |       nan         |
| tanimoto_torsions_go_mean          |          0.09759   |         0.29277   |
| tanimoto_maccs_go_mean             |          0.09759   |         0.29277   |
| tanimoto_morgan_go_median          |        nan         |       nan         |
| tanimoto_atompairs_go_median       |        nan         |       nan         |
| tanimoto_torsions_go_median        |          0.09759   |         0.29277   |
| tanimoto_maccs_go_median           |          0.09759   |         0.29277   |
| tanimoto_morgan_go_min             |        nan         |       nan         |
| tanimoto_atompairs_go_min          |        nan         |       nan         |
| tanimoto_torsions_go_min           |          0.09759   |         0.29277   |
| tanimoto_maccs_go_min              |          0.09759   |         0.29277   |
| tanimoto_morgan_go_max             |        nan         |       nan         |
| tanimoto_atompairs_go_max          |        nan         |       nan         |
| tanimoto_torsions_go_max           |          0.09759   |         0.29277   |
| tanimoto_maccs_go_max              |          0.09759   |         0.29277   |

### Scatter plots between pairs of scores

train/test scores, and sequence identity scores

### Score distribution

TODO Max, min, hist

## Next steps/TODOs

- graph with sample count depending on min sample count and max overlap count
- ML models with feature selection and/or PCA
- is_a is not working, since no two go terms in the ML dataset are direct descendants. 
  - Other graph similarity score, like distance of common ancestor, or whether they have the same parent
- Find way to implement other semantic similarity algorithms more efficiently (parallel)
