#!/usr/bin/env python3
import argparse
from serum_readfilter import makedb
from serum_readfilter import runfilter


def program_initialization():
    parser = argparse.ArgumentParser(
        prog='serum_readfilter',
        usage='serum_readfilter <mode> <options>',
        description='serum_readfilter quickly filters reads against a database of interest'
    )
    subparsers = parser.add_subparsers(
        title='Available modes',
        dest='mode'
    )

    add_subparser__makedb(subparsers)
    add_subparser__runfilter(subparsers)

    args = parser.parse_args()

    return args


def add_subparser__makedb(subparsers):
    makedb_parser = subparsers.add_parser(
        'makedb',
        description='Make database for run'
    )
    makedb_subparsers = makedb_parser.add_subparsers(
        title='Avaiable methods',
        dest='method'
    )
    add_subparser__makekrakendb(makedb_subparsers)
    add_subparser__makekaijudb(makedb_subparsers)


def add_subparser__makekrakendb(subparsers):
    makekrakendb_parser = subparsers.add_parser(
        'kraken',
        help='',
        usage='',
        description='serum_readfilter - make kraken db, generate a DB to be used for filtering of reads'
    )
    makekrakendb_parser.add_argument(
        "-db",
        "--db_location",
        help="Create a kraken db at location",
        required=True
    )
    makekrakendb_parser.add_argument(
        "-ref",
        "--fasta_location",
        help="Reference fasta file/directory containing sequences to filter on",
        required=True
    )
    makekrakendb_parser.add_argument(
        "-f",
        "--force_clean",
        help="Remove DB folder if present before",
        action="store_true",
        default=False
    )
    makekrakendb_parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads",
        default=1
    )
    makekrakendb_parser.add_argument(
        "-k",
        "--kmer_size",
        help="Kmer size for DB creation",
        default=31
    )
    makekrakendb_parser.add_argument(
        "-ext",
        "--file_extensions",
        nargs="+",
        type=str,
        help="Acceptable fasta file extenstions for reference",
        default=["fna", "fa", "fasta"]
    )


def add_subparser__makekaijudb(subparsers):
    makekaijudb_parser = subparsers.add_parser(
        'kaiju',
        help='',
        usage='',
        description='serum_readfilter - make kraken db, generate a DB to be used for filtering of reads'
    )
    makekaijudb_parser.add_argument(
        "-db",
        "--db_location",
        help="Create file <db_location>.fmi at location",
        required=True)
    makekaijudb_parser.add_argument(
        "-ref",
        "--fasta_location",
        help="Reference fasta file/directory containing sequences to filter on",
        required=True)
    makekaijudb_parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads",
        default=1)
    makekaijudb_parser.add_argument(
        "-ext",
        "--file_extensions",
        nargs="+",
        type=str,
        help="Acceptable fasta file extenstions for reference",
        default=["fna", "fa", "fasta"]
    )


def add_subparser__runfilter(subparsers):
    runfilter_parser = subparsers.add_parser(
        'runfilter',
        description='Run filter against a DB'
    )
    runfilter_subparsers = runfilter_parser.add_subparsers(
        title='Avaiable methods',
        dest='method'
    )
    add_subparser__runkrakenfilter(runfilter_subparsers)
    add_subparser__runkaijufilter(runfilter_subparsers)


def add_subparser__runkrakenfilter(subparsers):
    runfilter_parser = subparsers.add_parser(
        'kraken',
        help='',
        usage='',
        description='serum_readfilter - filter a set of reads against a database (default: kraken)'
    )
    runfilter_parser.add_argument(
        "-db", "--database_to_use",
        help="db at location",
        type=str,
        required=True
    )
    runfilter_parser.add_argument(
        "-R1",
        "--R1_reads",
        help="Fasta file containing sequences to filter on",
        type=str,
        required=True
    )
    runfilter_parser.add_argument(
        "-R2",
        "--R2_reads",
        help="Fasta file containing sequences to filter on",
        type=str,
        required=True
    )
    runfilter_parser.add_argument(
        "-o",
        "--output_name",
        help="Name to add to filtered samples",
        default="filtered"
    )
    runfilter_parser.add_argument(
        "-t",
        "--threads",
        help="threads",
        default=8
    )
    runfilter_parser.add_argument(
        "-inv",
        "--inverse",
        help="Returns reads that aren't in this set instead of ones that are",
        action="store_true",
        default=False
    )
    runfilter_parser.add_argument(
        "--norm",
        help="Uses bbnorm to normalize reads down to ensure maximum coverage on kmers",
        action="store_true",
        default=False
    )
    runfilter_parser.add_argument(
        "--config",
        help="Config file to overwrite",
        type=str,
        default=""
    )


def add_subparser__runkaijufilter(subparsers):
    runfilter_parser = subparsers.add_parser(
        'kaiju',
        help='',
        usage='',
        description='serum_readfilter - filter a set of reads against a database (default: kraken)'
    )
    runfilter_parser.add_argument(
        "-db", "--database_to_use",
        help="db at location",
        type=str,
        required=True
    )
    runfilter_parser.add_argument(
        "-R1",
        "--R1_reads",
        help="Fasta file containing sequences to filter on",
        type=str,
        required=True
    )
    runfilter_parser.add_argument(
        "-R2",
        "--R2_reads",
        help="Fasta file containing sequences to filter on",
        type=str
    )
    runfilter_parser.add_argument(
        "-o",
        "--output_name",
        help="Name to add to filtered samples",
        default="filtered"
    )
    runfilter_parser.add_argument(
        "-t",
        "--threads",
        help="threads",
        default=8
    )
    runfilter_parser.add_argument(
        "-inv",
        "--inverse",
        help="Returns reads that aren't in this set instead of ones that are",
        action="store_true",
        default=False
    )
    runfilter_parser.add_argument(
        "--norm",
        help="Uses bbnorm to normalize reads down to ensure maximum coverage on kmers",
        action="store_true",
        default=False
    )
    runfilter_parser.add_argument(
        "--config",
        help="Config file to overwrite",
        type=str,
        default=""
    )


def main():
    args = program_initialization()

    if args.mode == "makedb":
        makedb.make_db(args)
    if args.mode == "runfilter":
        runfilter.run_filter(args)


if __name__ == "__main__":
    main()
