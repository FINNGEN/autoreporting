import unittest
import unittest.mock as mock
import sys,os,json, requests
import pandas as pd
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
from Scripts.data_access import alleledb
from Scripts.data_access.db import Location, VariantData, Variant

class TestAlleleDB(unittest.TestCase):
    def test_vcf_init(self):
        """Test initialization

        Cases:
            No file
            Invalid file
            Proper file
        """
        nonexistent_file = ""
        proper_file = "testing/db_resources/test_vcf.vcf.gz"
        invalid_file = "testing/db_resources/invalid_file.gz"
        with self.assertRaises(FileNotFoundError) as ass:
            alleledb.VCFAlleleDB(nonexistent_file)
        with self.assertRaises(Exception) as ass:
            alleledb.VCFAlleleDB(invalid_file)
        alleledb.VCFAlleleDB(proper_file)

    def test_vcf_get_alleles(self):
        """Test getting alleles for positions

        Cases:
            No alleles for position
            biallelic position
            Multiallelic position
        """
        filename = "testing/db_resources/test_vcf.vcf.gz"
        db = alleledb.VCFAlleleDB(filename)
        #magic knowledge!
        no_alleles_loc = [Location("1",10)]
        biallelic_loc = [Location("1",10039)]
        multiallelic_loc = [Location("1",10169)]
        manylocs = [
            Location("1",10),
            Location("1",10039),
            Location("1",10169)
        ]

        no_alleles = db.get_alleles(no_alleles_loc)
        biallelic = db.get_alleles(biallelic_loc)
        multiallelic = db.get_alleles(multiallelic_loc)
        manyvars = db.get_alleles(manylocs)
        no_all_ver = []
        biallelic_ver = [VariantData(
            Variant("1",
            10039,
            "A",
            "C"),
            [],
            978760828
        )]
        multiallelic_ver = [
            VariantData(
            Variant("1",
            10169,
            "T",
            "C"),
            ["G"],
            1456517851
        )
        ]
        manyvar_ver = biallelic_ver + multiallelic_ver
        self.assertEqual(no_alleles,no_all_ver)
        self.assertEqual(biallelic,biallelic_ver)
        self.assertEqual(multiallelic,multiallelic_ver)
        self.assertEqual(manyvars,manyvar_ver)