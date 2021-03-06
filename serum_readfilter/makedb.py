import Bio.SeqIO
import os
import subprocess
import shutil
import gzip
import tempfile


def make_and_set_destination_folder(db_location):
    if not os.path.isdir(db_location):
        os.mkdir(db_location)
    os.chdir(db_location)


def get_fasta_records(fasta, file_extensions):
    records = []
    if os.path.isdir(fasta):
        for file in os.listdir(fasta):
            if os.path.isfile(os.path.join(fasta, file)) and "." in file and file.split(".")[1].strip() in file_extensions:
                print("{} being added".format(file))
                if file.endswith(".gz"):
                    with gzip.open(os.path.join(fasta, file), "rt") as fasta_input:
                        records.extend(list(Bio.SeqIO.parse(fasta_input, "fasta")))
                else:
                    with open(os.path.join(fasta, file), "r") as fasta_input:
                        records.extend(list(Bio.SeqIO.parse(fasta_input, "fasta")))
    elif os.path.isfile(fasta):
        with open(fasta, "r") as fasta_input:
            records = list(Bio.SeqIO.parse(fasta_input, "fasta"))
    else:
        print("No references found")
        return 1
    return records


def make_db(args):
    print(args.file_extensions)
    if args.method == "kraken":
        make_kraken_db_from_fasta(args.fasta_location, args.db_location, args.threads, args.kmer_size, args.file_extensions)
    elif args.method == "kaiju":
        make_kaiju_db_from_fasta(args.fasta_location, args.db_location, args.threads, args.file_extensions)
    return 0


def make_kraken_db_from_fasta(fasta_location, db_location, threads, kmer_size, file_extensions):
    if shutil.which("kraken-build") is None:
        print("Error finding kraken-build (from kraken) in PATH")
        return 1

    records = get_fasta_records(fasta_location, file_extensions)
    if records == 1:
        return 1

    previous_path = os.getcwd()
    make_and_set_destination_folder(db_location)

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
        Bio.SeqIO.write(records, output, "fasta")

    subprocess.call("kraken-build --threads {} --add-to-library kraken.fasta --db .".format(threads), shell=True)
    subprocess.call("kraken-build --threads {} --build --kmer-len {} --minimizer-len 1 --db .".format(threads, kmer_size), shell=True)

    os.chdir(previous_path)
    return 0


def make_kaiju_db_from_fasta(fasta_location, db_location, threads, file_extensions):
    if shutil.which("mkbwt") is None:
        print("Error finding mkbwt (from kaiju) in PATH")
        return 1
    if shutil.which("mkfmi") is None:
        print("Error finding mkfmi (from kaiju) in PATH")
        return 1

    temp_file = tempfile.NamedTemporaryFile(mode='w+t')
    records = get_fasta_records(fasta_location, file_extensions)

    Bio.SeqIO.write(records, temp_file, "fasta")

    subprocess.call("mkbwt -n {} -a protein -o {} {}".format(threads, db_location, temp_file.name), shell=True)
    subprocess.call("mkfmi {}".format(db_location), shell=True)
    os.remove(db_location + ".bwt")
    print("removed " + db_location + ".bwt")
    os.remove(db_location + ".sa")
    print("removed " + db_location + ".sa")
    temp_file.close()

    return 0
