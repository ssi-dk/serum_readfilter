#!/usr/bin/env python3
import Bio.SeqIO
import os
import subprocess
import shutil


def get_fasta_records(fasta, file_extensions):
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
    return records


def make_kraken_db_from_fasta(fasta, db_location, threads, kmer_size, file_extensions):
    if shutil.which("kraken-build") is None:
        print("Error finding kraken-build (from kraken) in PATH")
        return 1

    records = get_fasta_records(fasta, file_extensions)
    if records == 1:
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


def make_kaiju_db_from_fasta(fasta_file, db_location, threads):
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
