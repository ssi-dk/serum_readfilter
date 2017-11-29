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


def filter_reads_on_kraken(reads, outfile, db_location, threads, inverse):
    with open(os.devnull, 'w') as devnull:
        if not inverse:
            subprocess.Popen(["kraken", "--db", db_location, "--threads", threads, "--quick", "--min-hits", "1", "--classified-out", outfile, reads], stdout=devnull, stderr=subprocess.STDOUT)
        else:
            subprocess.Popen(["kraken", "--db", db_location, "--threads", threads, "--quick", "--min-hits", "1", "--unclassified-out", outfile, reads], stdout=devnull, stderr=subprocess.STDOUT)


if __name__ == "__main__":
    args = program_initialization()
    filter_reads_on_kraken(args.R1_reads, args.output_name + "_" + os.path.split(args.R1_reads)[1], args.database_to_use, args.threads, args.inverse)
    if args.R2_reads is not None:
        filter_reads_on_kraken(args.R2_reads, args.output_name + "_" + os.path.split(args.R2_reads)[1], args.database_to_use, args.threads, args.inverse)