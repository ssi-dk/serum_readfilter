# FastFilter-private

FastFilter is a program designed to filter whole genome seequence data to only obtain reads that may touch your gene's of interest. By creating a database of your gene's of interest then havign any reads that map 1 kmer entry to them you can reduce your data set substantially. This is very useful when you want to map reads against a gene such as with mlst, resistance finder, viruelence finder, etc.

Future additions
Add in BWT on AA from kaiju as alternative filtering approach