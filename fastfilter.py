#!/usr/bin/env python3
import subprocess
import argparse
import ruamel.yaml
import pkg_resources

# for f in ../../cge_dbs/resfinder_db/*.fsa; do (cat "${f}"; echo) >> resfinder.fasta; done
config_file = pkg_resources.resource_filename(__name__, "config.yaml")
yaml = ruamel.yaml.YAML(typ='safe')
yaml.default_flow_style = False
with open(config_file, "r") as yaml_stream:
    config = yaml.load(yaml_stream)


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
    parser.add_argument("--kaiju",
                        help="Uses Kaiju a BWT AA approach over Kraken",
                        action="store_true",
                        default=False)
    parser.add_argument("--norm",
                        help="Uses bbnorm to normalize reads down to ensure maximum coverage on kmers",
                        action="store_true",
                        default=False)
    parser.add_argument("--config",
                        help="Config file to overwrite",
                        type=str,
                        default="")
    args = parser.parse_args()

    return args


def filter_reads_on_kraken(R1_reads, R2_reads, outfile, db_location, threads, inverse):
    paired = "--paired --out-fmt paired "
    if R2_reads is None:
        R2_reads = ""
        paired = ""

    if not inverse:
        subprocess.call("kraken --db {} --threads {} {} {} --classified-out {} {} {} 1> /dev/null".format(db_location, threads, config["kraken"]["options"], paired, outfile, R1_reads, R2_reads), shell=True)
    else:
        subprocess.call("kraken --db {} --threads {} {} {} --unclassified-out {} {} {} 1> /dev/null".format(db_location, threads, config["kraken"]["options"], paired, outfile, R1_reads, R2_reads), shell=True)

    if R2_reads is None:
        return (outfile)
    return 0


def filter_reads_on_kaiju(R1_reads, R2_reads, outfile, db_location, threads, inverse):
    # subprocess.call("kaiju -t {} -f {} -i {} -j {} -z {} {}")
    return 0


def bbnorm_results(R1_reads, R2_reads, threads):
    R1_normalized = "norm" + R1_reads
    R2_normalized = "norm" + R2_reads
    if R2_reads is None:
        subprocess.call("bbnorm.sh threads={} in={} out={} {} 1> /dev/null".format(threads, R1_reads, R1_normalized, config["bbnorm"]["options"]), shell=True)
    else:
        subprocess.call("bbnorm.sh threads={} in={} in2={} out={} out2={} {} 1> /dev/null".format(threads, R1_reads, R2_reads, R1_normalized, R2_normalized, config["bbnorm"]["options"]), shell=True)


if __name__ == "__main__":
    args = program_initialization()
    print("Starting")

    if args.config != "":
        with open(args.config, "r") as yaml_stream:
            config = yaml.load(yaml_stream)

    if args.kaiju:
        filter_reads_on_kaiju()
    else:
        filter_reads_on_kraken(args.R1_reads, args.R2_reads, args.output_name, args.database_to_use, args.threads, args.inverse)
    if args.norm:
        if args.R2_reads is None:
            bbnorm_results(args.output_name + "_R1.fastq", args.R2_reads, args.threads)
        else:
            bbnorm_results(args.output_name + "_R1.fastq", args.output_name + "_R2.fastq", args.threads)
    print("Complete")
