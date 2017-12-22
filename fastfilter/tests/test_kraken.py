import unittest
import shutil
import os

from fastfilter import makedb
from fastfilter import runfilter


class fastfilterTests(unittest.TestCase):

    def setUp(self):
        pass

    def test_make_db(self):
        makedb.make_kraken_db_from_fasta("data/cdifficile_mlst", "cdifficile_db", 1, 31, ["tfa"])
        # self.assertTrue(os.path.isdir("cdifficile_db"))

    def test_runfilter(self):
        runfilter.filter_reads_on_kraken("data/reads/cdifficile_R1.fastq.gz", "data/reads/cdifficile_R1.fastq.gz", "filtered", "cdifficile_db", 1, False)
        # self.assertTrue(os.path.isfile("filtered_R1.fastq"))
        # self.assertTrue(os.path.isfile("filtered_R2.fastq"))

    def tearDown(self):
        # shutil.rmtree("cdifficile_db")
        # os.remove("filtered_R1.fastq")
        # os.remove("filtered_R2.fastq")
        pass


if __name__ == '__main__':
    unittest.main()
