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
                        help="Fasta file/directory containing sequences to filter on",
                        required=True)
    parser.add_argument("-f", "--force_clean",
                        help="Remove DB folder if present before",
                        action="store_true",
                        default=False)
    parser.add_argument("-t", "--threads",
                        help="Number of threads",
                        default=1)
    parser.add_argument("-k", "--kmer_size",
                        help="Kmer size for DB creation",
                        default=31)
    parser.add_argument("-x", "--file_extensions",
                        help="Acceptable fasta file extenstions for reference",
                        default=".fna,.fa,.fasta")
    args = parser.parse_args()

    args.file_extensions = args.file_extensions.split(",")

    return args


def make_kraken_db_from_fasta(fasta, db_location, threads, kmer_size, file_extensions):
    if shutil.which("kraken-build") is None:
        print("Error finding kraken-build (from kraken) in PATH")
        return 1

    records = []
    if os.path.isdir(fasta):
        for file in os.listdir(fasta):
            if os.path.isfile(os.path.join(fasta, file)) and file.split["."][1] in file_extensions:
                with open(os.path.join(fasta, file), "r") as fasta_input:
                    records.extend(list(Bio.SeqIO.parse(fasta_input), "fasta"))
    elif os.path.isfile(fasta):
        with open(fasta, "r") as fasta_input:
            records = list(Bio.SeqIO.parse(fasta_input, "fasta"))
    else:
        print("No references found")
        return 1

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

    return 0


if __name__ == "__main__":
    args = program_initialization()
    print("Starting")
    if args.force_clean and os.path.isdir(args.database_to_create):
        shutil.rmtree(args.database_to_create)
    make_kraken_db_from_fasta(args.input_fasta, args.database_to_create, args.threads, args.kmer_size, args.file_extensions)
    print("Complete")
