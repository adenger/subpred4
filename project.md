# Manuscript 2

## General ideas

Problems with methods from manuscript 1:

- Not that many substrate classes
    - Increase sample count for training/testing data
        - Explore new data sources
    - Find more optimal substrate classes through clustering/unsupervised learning
        - If it works well: Use unsupervised learning for prediction, compare to SVM
    - Improve accuracy for classes with few samples
        - More sophisticated feature generation algorithms: Other PSSM encodings, protein embeddings, etc.
- Complicated setup
    - Provide pre-trained model that is easily applicable
    - Or snakemake pipeline, python package, webserver

## Dataset improvements

**Uniprot & TrEMBL**

We updated the Uniprot version from version 2021_04 to 2022_05. This led to a considerable increase in membrane proteins and transmembrane transporters. 

In manuscript 1, we only used manually curated SwissProt data. Now we want to include proteins where the sequence existence has been verified at protein level (or at transcript level?).

We also added more columns to the dataset. They might be useful when trying to find annotations that match the clusters in the data. The new annotations include ChEBI terms, for which we implemented a parser for the OWL file. Another addition are Interpro domains. The Interpro database now includes Pfam.

*TODO stats on changes between versions.*

**GO terms**

We will start using GO terms instead of keywords. GO-MF annotations are more precise, since they relate to the function of a specific protein, not its upstream or downstream function, and since more terms exist overall. 

Keywords relate more to the biological process that a protein is involved in, and we often need to create intersections between Keyword sets in order to get a specific set of proteins, which can lead to edge cases. For example, the GO dataset contains a term "transmembrane transporter activity", while we have to combine the "Transport" and "Transmembrane" keywords, which can include peripheral membrane proteins that are involved in transport, for example.

The number of GO annotations can be significantly increased by using the OWL file to retrieve all child terms of a particular term. While there might not be enough transmembrane transporters annotated with "transmembrane transporter activity", there are many more that are annotated with a descendant of that term.

We have thoroughly compared the GO annotations from muliple sources:

- The goa_uniprot_all file from gene-ontology.org
- The goa_uniprot_all file from the EBI FTP server
- The GO column from a Uniprot custom download

Strangely, the gene-ontology file had little overlap with the Uniprot custom download annotations. 

The problem with the custom download annotations is that the qualifiers and the evidence codes are missing. To get accurate result, we should stick to the "enables" qualifier and the non-IEA evidence codes.  

Overall, the EBI file was by far the most complete. Linking the custom download file with the qualifiers and evidence codes from the EBI file yielded fewer proteins and annotations than simply filtering the EBI file for Uniprot identifiers in the Uniprot dataset. Although there might be a reason why Uniprot removed some of the annotations?

Link: https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz

**PSSM dataset**

We also updated the Uniref datasets used for local BLAST searches to version 2022_05. Some sequences have changed, and this leads to sequences in the PSSM files being different from the sequences in the Uniprot dataset.

*PSSMs were generated again from scratch.*


## New methods

### Clustering/unsupervised learning


First, we tried special tools for sequence clustering, such as MMSeq2 and CD-Hit. They both encode and cluster the proteins. This approach did not work, since they create too many clusters. Even with the least strict clustering parameters, more than half the proteins were in their own cluster. We want methods that create around 5-20 clusters, not hundreds.

Clustering a protein sequence dataset follows two steps: Encoding the sequence into a vector of fixed length *n*, and the clustering of those vectors. Finally, the clusters are compared to training labels, which are given by annotations such as GO terms.

**Encodings:**

- AAC, PAAC
- PSSM
    - Improve algorithms? *TODO*
- Biovec
- Word2Vec custom algorithm
- ProtNLM

**Clustering**

Algorithms:

- Kmeans

Finding the optimal number of clusters k: Cluster quality metrics.

- Silhouette Coefficient
- Calinski-Harabasz Index
- Davis-Bouldin Index
- Elbow plot

**Assessment of clusters**

Gene annotation enrichment analysis: 

- Hypergeometric test
- Rand metric?
