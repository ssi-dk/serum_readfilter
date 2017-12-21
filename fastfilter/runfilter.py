import subprocess
import ruamel.yaml
import pkg_resources
import Bio.SeqIO
import gzip
import shutil
import os
import tempfile


# for f in ../../cge_dbs/resfinder_db/*.fsa; do (cat "${f}"; echo) >> resfinder.fasta; done
config_file = pkg_resources.resource_filename(__name__, "config.yaml")
yaml = ruamel.yaml.YAML(typ='safe')
yaml.default_flow_style = False
with open(config_file, "r") as yaml_stream:
    config = yaml.load(yaml_stream)


def set_config(config_file):
    global config
    with open(config_file, "r") as yaml_stream:
        config = yaml.load(yaml_stream)

    return 0


def filter_reads_on_kraken(R1_reads, R2_reads, outfile, db_location, threads, inverse):
    if shutil.which("kraken") is None:
        print("Error finding kraken in PATH")
        exit()

    R1_params = R1_reads
    R2_params = R2_reads + " --paired --out-fmt paired" if R2_reads is not None else ""
    out_param = "--classified-out " + outfile if not inverse else "--unclassified-out" + outfile

    subprocess.call("kraken --db {} --threads {} {} {} {} {} 1> /dev/null".format(db_location, threads, config["kraken"]["options"], out_param, R1_params, R2_params), shell=True)

    return 0


def filter_reads_on_kaiju(R1_reads, R2_reads, outfile, db_location, threads, inverse):
    if shutil.which("kaijux") is None:
        print("Error finding kaijux in PATH")
        exit()

    R1_params = "-i " + R1_reads
    R2_params = "-j " + R2_reads if R2_reads is not None else ""

    temp = tempfile.TemporaryFile()

    subprocess.call("kaijux -z {} -f {} {} -o {} {}".format(threads, db_location, R1_params, R2_params, temp.name, config["kaiju"]["options"]), shell=True)

    classified_letter = "C"
    if inverse:
        classified_letter = "U"
    reads = [R1_reads, R2_reads] if R2_reads is not None else [R1_reads]
    extract_reads(temp.name, reads, classified_letter, outfile)
    os.remove(temp.name)

    return 0


def extract_reads(classifier_file, reads, classification_symbol, outfile):
    read_dict = {}
    with open(classifier_file, "r") as classifier:
        for line in classifier:
            if line.startswith(classification_symbol):
                read_dict[line.split("\t")[1]] = ""

    for i, read in enumerate(reads):
        filtered_records = []
        for record in Bio.SeqIO.parse(gzip.open(read, "rt"), "fastq"):
            if record.id in read_dict:
                filtered_records.append(record)
        with open(outfile + "_R" + str(i + 1) + ".fastq", "w") as read_out:
            Bio.SeqIO.write(filtered_records, read_out, "fastq")

    return 0


def bbnorm_results(R1_reads, R2_reads, threads):
    if shutil.which("bbnorm.sh") is None:
        print("Error finding bbnorm.sh in PATH")
        exit()

    R1_params = "in={} out={}".format(R1_reads, "norm" + R1_reads)
    R2_params = "in2={} out2={}".format(R2_reads, "norm" + R2_reads) if R2_reads is not None else ""
    subprocess.call("bbnorm.sh threads={} {} {} {} 1> /dev/null".format(threads, R1_params, R2_params, config["bbnorm"]["options"]), shell=True)
    return 0
