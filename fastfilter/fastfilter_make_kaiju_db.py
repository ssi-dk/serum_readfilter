#!/usr/bin/env python3
import subprocess
import argparse
import os
import shutil

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
    if shutil.which("mkbwt") is None:
        print("Error finding mkbwt (from kaiju) in PATH")
        return 1
    if shutil.which("mkfmi") is None:
        print("Error finding mkfmi (from kaiju) in PATH")
        return 1
    subprocess.call("mkbwt -n {} -a protein -o {} {}".format(threads, db_location, fasta_file), shell=True)
    subprocess.call("mkfmi {}".format(db_location), shell=True)
    os.remove(db_location + ".bwt")
    os.remove(db_location + ".sa")

    return 0


if __name__ == "__main__":
    args = program_initialization()
    print("Starting")
    make_kaiju_db_from_fasta(args.input_fasta, args.database_to_create, args.threads)
    print("Complete")
