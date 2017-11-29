#!/usr/bin/env python3
import os
import subprocess
import argparse

# for f in ../../cge_dbs/resfinder_db/*.fsa; do (cat "${f}"; echo) >> resfinder.fasta; done


def program_initialization():
    parser = argparse.ArgumentParser(description='Fastfilter - makedb, generate a DB to be used for filtering of reads')
    parser.add_argument("-db", "--database_to_use",
                        help="Create a kraken db at location",
                        type=str,
                        required=True)
    parser.add_argument("-R1", "--R1_reads",
                        help="Fasta file containing sequences to filter on",
                        type=str,
                        required=True)
    parser.add_argument("-R2", "--R2_reads",
                        help="Fasta file containing sequences to filter on",
                        type=str)
    parser.add_argument("-o", "--output_name",
                        help="Name to add to filtered samples",
                        default="filtered")
    parser.add_argument("-t", "--threads",
                        help="threads",
                        default=8)
    parser.add_argument("-inv", "--inverse",
                        help="Returns reads that aren't in this set instead of ones that are",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    return args


def filter_reads_on_kraken(R1_reads, R2_reads, outfile, db_location, threads, inverse):
    if R2_reads is None:
        R2_reads = ""
    if not inverse:
        subprocess.call("kraken --db {} --threads {} --quick --min-hits 1 --classified-out {} {} {} 1> /dev/null".format(db_location, threads, outfile, R1_reads, R2_reads), shell=True)
    else:
        subprocess.call("kraken --db {} --threads {} --quick --min-hits 1 --unclassified-out {} {} {} 1> /dev/null".format(db_location, threads, outfile, R1_reads, R2_reads), shell=True)


if __name__ == "__main__":
    args = program_initialization()
    print("Starting")
    filter_reads_on_kraken(args.R1_reads, args.R2_reads, args.output_name, args.database_to_use, args.threads, args.inverse)
    print("Complete")
