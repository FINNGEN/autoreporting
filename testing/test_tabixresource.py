import unittest
import sys
import pysam
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
from Scripts.load_tabix import tb_resource_manager,TabixResource
from Scripts.data_access.db import Variant

class TestTBResource(unittest.TestCase):
    def test_creating_tabixresource(self):
        with tb_resource_manager("testing/db_resources/test_vcf.vcf.gz","#CHROM","POS","REF","ALT") as tbres:
            pass

    def test_read_region(self):
        with tb_resource_manager("testing/db_resources/test_vcf.vcf.gz","#CHROM","POS","REF","ALT") as tbres:
            region = tbres.load_region("1",10018,10039,["ID"])
            validation = {
                Variant("1",10019,"TA","T"):["rs775809821"],
                Variant("1",10039,"A","C"):["rs978760828"]
            }
            self.assertEqual(region,validation)

    def test_header_skip_file(self):
        with tb_resource_manager("testing/db_resources/skip_header.tsv.gz","chrom","pos","ref","alt") as tbres:
            #validate header
            self.assertEqual(
                tbres.header,
                ['chrom', 'pos', 'ref', 'alt', 'data']
            )
            
            region = tbres.load_region("1",2,3,["data"])
            validation = {
                Variant("1",2,"T","CAA"):["20"],
                Variant("1",3,"T","CAAA"):["30"]
            }
            self.assertEqual(region,validation)