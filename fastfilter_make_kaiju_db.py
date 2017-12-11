#!/usr/bin/env python3
import subprocess
import argparse

# for f in ../../cge_dbs/resfinder_db/*.fsa; do (cat "${f}"; echo) >> resfinder.fasta; done


def program_initialization():
    parser = argparse.ArgumentParser(description='Fastfilter - make kaiju db, generate a DB to be used for filtering of reads')
    parser.add_argument("-db", "--database_to_create",
                        help="Create a kaiju.fmi db at location",
                        required=True)
    parser.add_argument("-i", "--input_fasta",
                        help="Fasta file containing sequences to filter on",
                        required=True)
    parser.add_argument("-t", "--threads",
                        help="Number of threads",
                        default=1)
    args = parser.parse_args()

    return args


def make_kaiju_db_from_fasta(fasta_file, db_location, threads=1):
    subprocess.call("mkbwt -n {} -a protein -o {}".format(threads, db_location), shell=True)
    subprocess.call("mkfmi {}".format(db_location), shell=True)


if __name__ == "__main__":
    args = program_initialization()
    print("Starting")
    make_kaiju_db_from_fasta(args.input_fasta, args.database_to_create, args.threads)
    print("Complete")
