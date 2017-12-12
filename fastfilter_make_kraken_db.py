#!/usr/bin/env python3
import Bio.SeqIO
import os
import subprocess
import argparse
import shutil

# for f in ../../cge_dbs/resfinder_db/*.fsa; do (cat "${f}"; echo) >> resfinder.fasta; done


def program_initialization():
    parser = argparse.ArgumentParser(description='Fastfilter - make kraken db, generate a DB to be used for filtering of reads')
    parser.add_argument("-db", "--database_to_create",
                        help="Create a kraken db at location",
                        required=True)
    parser.add_argument("-i", "--input_fasta",
                        help="Fasta file containing sequences to filter on",
                        required=True)
    parser.add_argument("-f", "--force_clean",
                        help="Remove DB folder if present before",
                        action="store_true",
                        default=False)
    parser.add_argument("-t", "--threads",
                        help="Number of threads",
                        default=1)
    parser.add_argument("-k", "--kmer_size",
                        help="Kmer size for DB creation default 21",
                        default=21)
    args = parser.parse_args()

    return args


def make_kraken_db_from_fasta(fasta_file, db_location, threads=1, kmer_size=21):
    if shutil.which("kraken-build") is None:
        print("Error finding kraken-build (from kraken) in PATH")
        exit()

    with open(fasta_file, "r") as fasta_input:
        records = list(Bio.SeqIO.parse(fasta_input, "fasta"))

    if not os.path.isdir(db_location):
        os.mkdir(db_location)
    os.chdir(db_location)

    for i, record in enumerate(records):
        records[i].description = ""
        records[i].name = ""
        records[i].id = "kraken:taxid|1"

    os.mkdir("taxonomy")
    with open("taxonomy/names.dmp", "w") as output:
        output.write("1\t|\tPlaceholder\t|\t\t|\tscientific name\t|\n")

    with open("taxonomy/nodes.dmp", "w") as output:
        output.write("1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|\n")

    with open("taxonomy/gi_taxid_nucl.dmp", "w") as output:
        output.write("1\t1\n")

    with open("seqid2taxid.map", "w") as output:
        output.write("kraken:taxid|1\t1\n")

    with open("kraken.fasta", "w") as output:
        Bio.SeqIO.write(records, "kraken.fasta", "fasta")

    subprocess.call("kraken-build --threads {} --add-to-library kraken.fasta --db .".format(threads), shell=True)
    subprocess.call("kraken-build --threads {} --build --kmer-len {} --minimizer-len 1 --db .".format(threads, kmer_size), shell=True)


if __name__ == "__main__":
    args = program_initialization()
    print("Starting")
    if args.force_clean and os.path.isdir(args.database_to_create):
        shutil.rmtree(args.database_to_create)
    make_kraken_db_from_fasta(args.input_fasta, args.database_to_create, args.threads, args.kmer_size)
    print("Complete")
