# Manuscript 2

## General ideas & goals

- More classes
  - Requires more samples per class
  - How much should the dataset be divided? Evaluate.
- More samples
  - Explore new methods of finding substrate annotations
  - Include less accurate or predicted data into the training
- More automation
  - Pre-trained model or container
  - Input: Fasta. Output: Transporters and their substrates.

## Datasets & Preprocessing

- Protein data
  - Uniprot accession
  - Gene Names
  - Protein names
  - Reviewed status
  - Protein existence
  - Organism ID
  - Sequence
- Ontologies
  - GO
  - ChEBI
- Annotations
  - Uniprot to GO
  - Uniprot to Interpro
  - Uniprot to TCDB
  - Uniprot to Keywords
  - GO to ChEBI (from QuickGO)
  - TCDB to ChEBI
- Blast databases (Uniref50/90)

## Graph analysis GO & ChEBI

Goal:

Finding substrate classes with enough samples, maybe automatically. Clustering based on features did not work, so we will try ontologies and protein annotations to divide the dataset into clusters.

Results:

- Transmembrane transporter GO terms in the dataset
- Set of primary substrate ChEBI terms associated with the active transmembrane transport GO terms
- Overlap of Chebi terms, in terms of proteins that are annotated with them
  - Co-transport and logical relations can be reasons

Next:

Aufspaltungen im Graphen analysieren. Wie gut funktionieren die? Wie fein können wir die Proteine unterteilen? Gibt es korrelationen?

## ML models on levels of tree

Idee:

- Auftrennen der Kategorien in den Graphen mit ML Modellen
  - Einfache SVM Pipeline
  - Viele verschiedene Test cases
- Wann funktioniert das? Vergleichen mit verschiedenen Maßen
  - Overlap der Samples
  - GO term semantic similarity
  - maximale Sequenzidentität zwischen Proteinen beider Gruppen
- Vergleich Meta-Datensatz mit einzelnen Organismen
  - Und Gruppen ähnlicher Spezies wie Hefen

Data results:

- For each (organism, substrate1, substrate2) multiple matrices/heatmaps:
  - is_a connections (nx.from_pandas_adjacency, nx.adjacency_matrix)
  - SVM scores for one model (select one)
    - F1 train score for each substrate (or average)
    - F1 test score for each substrate (or average)
  - Overlap matrix
  - average sequence similarity
  - chemical similarity (smiles etc.)

- For each (organism, go_term1, go_term2) multiple matrices/heatmaps:
  - is_a connections (nx.from_pandas_adjacency, nx.adjacency_matrix)
  - SVM scores for one model (select one)
    - F1 train score for each substrate (or average)
    - F1 test score for each substrate (or average)
  - Overlap matrix
  - average sequence similarity
  - semantic similarity

- Then: calculate correlations, statistical tests, etc.
- For adjacency matrix: If substrates are 1 in a column (or row?) then they have same parent

### SVM Pipeline results
