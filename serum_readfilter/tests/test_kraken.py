#!/usr/bin/env python3

import unittest
import shutil
import os
import tempfile

from serum_readfilter import makedb
from serum_readfilter import runfilter


class serum_readfilterTests(unittest.TestCase):

    def setUp(self):
        self.here = os.path.dirname(os.path.abspath(__file__))
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)

    def test_make_db_and_filter(self):
        makedb.make_kraken_db_from_fasta(os.path.join(self.here, "data/cdifficile_mlst"), "cdifficile_db", 1, 31, ["tfa"])
        self.assertTrue(os.path.isdir("cdifficile_db"))
        self.assertTrue(os.path.isfile("cdifficile_db/lca.complete"))
        runfilter.filter_reads_on_kraken(os.path.join(self.here, "data/reads/cdifficile_R1.fastq.gz"), os.path.join(self.here, "data/reads/cdifficile_R1.fastq.gz"), "filtered", "cdifficile_db", 1, False)
        self.assertTrue(os.path.isfile("filtered_R1.fastq"))
        self.assertTrue(os.path.isfile("filtered_R2.fastq"))

    def tearDown(self):
        shutil.rmtree(self.test_dir)


if __name__ == '__main__':
    unittest.main()
