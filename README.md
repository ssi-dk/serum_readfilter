# serum_readfilter-private

serum_readfilter is a program designed to filter whole genome seequence data to only obtain reads that may touch your gene's of interest. By creating a database of your gene's of interest then havign any reads that map 1 kmer entry to them you can reduce your data set substantially while still retaining almost all options that would normally be a hit. This is very useful when you want to map reads against a gene such as with mlst, resistance finder, viruelence finder, etc. or have a large set of data to work with.

To use:

serum_readfilter makedb kraken -db <database> -i <reference.fasta>
serum_readfilter runfilter kraken -db <database> -R1 <R1_reads.fastq.gz> -R2 <R2_reads.fastq.gz>
