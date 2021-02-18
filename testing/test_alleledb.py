import unittest
import unittest.mock as mock
import sys,os,json, requests
import pandas as pd
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
from Scripts.data_access import alleledb
from Scripts.data_access.db import Location, VariantData

class TestAlleleDB(unittest.TestCase):
    def test_init(self):
        """Test initialization

        Cases:
            No file
            Invalid file
            Proper file
        """
        nonexistent_file = ""
        proper_file = "testing/db_resources/allele_db.gz"
        invalid_file = "testing/db_resources/invalid_file.gz"
        with self.assertRaises(FileNotFoundError) as ass:
            alleledb.FGAlleleDB(nonexistent_file)
        with self.assertRaises(Exception) as ass:
            alleledb.FGAlleleDB(invalid_file)
        alleledb.FGAlleleDB(proper_file)

    def test_get_alleles(self):
        """Test getting alleles for positions

        Cases:
            No alleles for position
            biallelic position
            Multiallelic position
        """
        filename = "testing/db_resources/allele_db.gz"
        db = alleledb.FGAlleleDB(filename)
        #magic knowledge!
        no_alleles_loc = [Location("1",1)]
        biallelic_loc = [Location("1",12)]
        multiallelic_loc = [Location("2",15)]
        manylocs = [
            Location("1",1),
            Location("1",12),
            Location("2",15)
        ]

        no_alleles = db.get_alleles(no_alleles_loc)
        biallelic = db.get_alleles(biallelic_loc)
        multiallelic = db.get_alleles(multiallelic_loc)
        manyvars = db.get_alleles(manylocs)
        no_all_ver = []
        biallelic_ver = [VariantData(
            "1",
            12,
            "A",
            ["C"],
            True
        )]
        multiallelic_ver = [
            VariantData(
            "2",
            15,
            "A",
            ["T","CT"],
            False
        )
        ]
        manyvar_ver = biallelic_ver + multiallelic_ver
        self.assertEqual(no_alleles,no_all_ver)
        self.assertEqual(biallelic,biallelic_ver)
        self.assertEqual(multiallelic,multiallelic_ver)
        self.assertEqual(manyvars,manyvar_ver)