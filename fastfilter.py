#!/usr/bin/env python3
import subprocess
import argparse
import ruamel.yaml
import pkg_resources
import skbio
import gzip

# for f in ../../cge_dbs/resfinder_db/*.fsa; do (cat "${f}"; echo) >> resfinder.fasta; done
config_file = pkg_resources.resource_filename(__name__, "config.yaml")
yaml = ruamel.yaml.YAML(typ='safe')
yaml.default_flow_style = False
with open(config_file, "r") as yaml_stream:
    config = yaml.load(yaml_stream)


def program_initialization():
    parser = argparse.ArgumentParser(description='Fastfilter - filter a set of reads against a database (default: kraken)')
    parser.add_argument("-db", "--database_to_use",
                        help="db at location",
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
    # subprocess.call("kaijux -t {} -f {} -i {} -j {} -z {} {}")
    tempfile = "kaiju.out"
    if R2_reads is None:
        subprocess.call("kaijux -z {} -f {} -i {} -o {} {}".format(threads, db_location, R1_reads, tempfile, config["kaiju"]["options"]), shell=True)
    else:
        subprocess.call("kaijux -z {} -f {} -i {} -j {} -o {} {}".format(threads, db_location, R1_reads, R2_reads, tempfile, config["kaiju"]["options"]), shell=True)

    if not inverse:
        extract_reads(tempfile, R1_reads, R2_reads, "C", outfile)
    else:
        extract_reads(tempfile, R1_reads, R2_reads, "U", outfile)
    return 0


def extract_reads(classifier_file, R1_reads, R2_reads, classification_symbol, outfile):
    read_dict = {}
    with open(classifier_file, "r") as classifier:
        for line in classifier:
            if line.startswith(classification_symbol):
                read_dict[line.split("\t")[1]] = ""
    with open(outfile + "_R1.fastq", "w") as R1_out:
        for seq in skbio.io.read(R1_reads, format="fastq"):
            if seq.metadata['id'] in read_dict:
                seq.write(R1_out, format="fastq")
    if R2_reads is not None:
        with open(outfile + "_R2.fastq", "w") as R2_out:
            for seq in skbio.io.read(R2_reads, format="fastq"):
                if seq.metadata['id'] in read_dict:
                    seq.write(R2_out, format="fastq")
    return 0


def bbnorm_results(R1_reads, R2_reads, threads):
    R1_normalized = "norm" + R1_reads
    R2_normalized = "norm" + R2_reads
    if R2_reads is None:
        subprocess.call("bbnorm.sh threads={} in={} out={} {} 1> /dev/null".format(threads, R1_reads, R1_normalized, config["bbnorm"]["options"]), shell=True)
    else:
        subprocess.call("bbnorm.sh threads={} in={} in2={} out={} out2={} {} 1> /dev/null".format(threads, R1_reads, R2_reads, R1_normalized, R2_normalized, config["bbnorm"]["options"]), shell=True)
    return 0


if __name__ == "__main__":
    args = program_initialization()
    print("Starting")

    if args.config != "":
        with open(args.config, "r") as yaml_stream:
            config = yaml.load(yaml_stream)

    if args.kaiju:
        filter_reads_on_kaiju(args.R1_reads, args.R2_reads, args.output_name, args.database_to_use, args.threads, args.inverse)
    else:
        filter_reads_on_kraken(args.R1_reads, args.R2_reads, args.output_name, args.database_to_use, args.threads, args.inverse)
    if args.norm:
        if args.R2_reads is None:
            bbnorm_results(args.output_name + "_R1.fastq", args.R2_reads, args.threads)
        else:
            bbnorm_results(args.output_name + "_R1.fastq", args.output_name + "_R2.fastq", args.threads)

    print("Complete")
