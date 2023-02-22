## Download links:

ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz

## Unzipping files:

conda activate subpred4

cd manuscript2/subpred4/data/raw/uniref/uniref50/
xz -T0 -d uniref50.fasta.xz && makeblastdb -in uniref50.fasta -parse_seqids -dbtype prot

cd manuscript2/subpred4/data/raw/uniref/uniref90/
xz -T0 -d uniref90.fasta.xz && makeblastdb -in uniref90.fasta -parse_seqids -dbtype prot
