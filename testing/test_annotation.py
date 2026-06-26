import unittest
import sys
import pysam
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
#from Scripts.load_tabix import tb_resource_manager,TabixResource
from Scripts.group_annotation import generate_chrom_ranges
from Scripts.data_access.db import Variant
from Scripts.autoreporting_utils import Region

class TestAnnotations(unittest.TestCase):

    def test_generate_ranges(self):
        test_data = [Variant(str(i//15+1),i,"A","C") for i in range(1,31)]
        out = generate_chrom_ranges(set(test_data))
        validation = {'1': [Region(chrom='1', start=0, end=14)],
            '2': [Region(chrom='2', start=14, end=29)],
            '3': [Region(chrom='3', start=29, end=30)]}
        self.assertEqual(out,validation)
        out2 = generate_chrom_ranges(set(test_data),maximum_range_length=4)
        validation2 = {'2': [
                Region(chrom='2', start=14, end=17),
                Region(chrom='2', start=17, end=20),
                Region(chrom='2', start=20, end=23),
                Region(chrom='2', start=23, end=26),
                Region(chrom='2', start=26, end=29)],
            '1': [
                Region(chrom='1', start=0, end=3),
                Region(chrom='1', start=3, end=6),
                Region(chrom='1', start=6, end=9),
                Region(chrom='1', start=9, end=12),
                Region(chrom='1', start=12, end=14)],
            '3' :[Region(chrom='3', start=29, end=30)]}
        self.assertEqual(out2,validation2)




if __name__=="__main__":
    unittest.main()