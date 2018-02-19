# serum_readfilter

serum_readfilter is a program designed to filter whole genome sequence data to only obtain reads that may touch your gene's of interest. This is accomplished by creating a database off a set of sequences you want to filter by (ie MLST, resistance genes, etc) and then filtering (via Kraken https://github.com/DerrickWood/kraken or Kaiju https://github.com/bioinformatics-centre/kaiju) the raw reads against this database. By default with Kraken this will filter all reads without 1 k-mer match to the database which can reduce your data set substantially while still retaining almost all options that would normally be a potential match to the database. This is very useful when you want to map reads against a gene or have a large set of data to work with and want to make it more managable.

To use:

serum_readfilter makedb kraken -db <database> -i <reference.fasta/or directory of fasta>

serum_readfilter runfilter kraken -db <database> -R1 <R1_reads.fastq.gz> -R2 <R2_reads.fastq.gz>


bioRxiv link:
https://www.biorxiv.org/content/early/2018/02/15/266080