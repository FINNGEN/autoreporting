import unittest
import sys
sys.path.append("../")
sys.path.append("./")
import pandas as pd
from Scripts import gws_fetch


class TestGws(unittest.TestCase):

    def test_create_variant_column(self):
        df_data={"#chrom":"1","pos":1234,"ref":"CGC","alt":"TTC"}
        df=pd.DataFrame(df_data,index=[0])
        variant=pd.Series("chr{}_{}_{}_{}".format("1",1234,"CGC","TTC"))
        variant_2=gws_fetch.create_variant_column(df)
        value=(variant==variant_2).all()
        self.assertTrue(value)

    def test_region_pruning(self):
        #case no regions
        reg=pd.DataFrame([])
        df=pd.DataFrame(columns=["#chrom", "pos_rmax","pos_rmin"])
        reg2=gws_fetch.prune_regions(df)
        #value=(reg==reg2).all().all()
        self.assertTrue(reg.equals(reg2))
        #case one region
        reg=pd.DataFrame({"#chrom":"1","max":20,"min":10},index=[0])
        df=pd.DataFrame({"#chrom":"1", "pos_rmin":10,"pos_rmax":20},index=[0])
        reg2=gws_fetch.prune_regions(df)
        self.assertTrue(reg.equals(reg2))
        #case two regions: no overlap by position
        reg=pd.DataFrame({"#chrom":["1","1"],"max":[20,40],"min":[10,30]})
        df=pd.DataFrame({"#chrom":["1","1"],"pos_rmax":[20,40],"pos_rmin":[10,30]})
        reg2=gws_fetch.prune_regions(df)
        self.assertTrue(reg.equals(reg2))
        #case two regions: overlapping positions, different chromosomes
        reg=pd.DataFrame({"#chrom":["1","2"],"max":[40,30],"min":[10,20]})
        df=pd.DataFrame({"#chrom":["1","2"],"pos_rmax":[40,30],"pos_rmin":[10,20]})
        reg2=gws_fetch.prune_regions(df)
        self.assertTrue(reg.equals(reg2))
        #case two regions: one in other, i.e. 1_min 2_min 2_max 1_max
        reg=pd.DataFrame({"#chrom":["1"],"max":[40],"min":[10]})
        df=pd.DataFrame({"#chrom":["1","1"],"pos_rmax":[40,20],"pos_rmin":[10,30]})
        reg2=gws_fetch.prune_regions(df)
        self.assertTrue(reg.equals(reg2))
        #case two regions: 1_min 2_min 1_max 2_max
        reg=pd.DataFrame({"#chrom":["1"],"max":[40],"min":[10]})
        df=pd.DataFrame({"#chrom":["1","1"],"pos_rmax":[30,40],"pos_rmin":[10,15]})
        reg2=gws_fetch.prune_regions(df)
        self.assertTrue(reg.equals(reg2))
        #case two regions: 2_min 1_min 2_max 1_max
        reg=pd.DataFrame({"#chrom":["1"],"max":[40],"min":[10]})
        df=pd.DataFrame({"#chrom":["1","1"],"pos_rmax":[30,40],"pos_rmin":[10,15]})
        reg2=gws_fetch.prune_regions(df)
        self.assertTrue(reg.equals(reg2))
        #case two regions: 2_min 1_min 1_max 2_max
        reg=pd.DataFrame({"#chrom":["1"],"max":[40],"min":[10]})
        df=pd.DataFrame({"#chrom":["1","1"],"pos_rmax":[30,40],"pos_rmin":[15,10]})
        reg2=gws_fetch.prune_regions(df)
        self.assertTrue(reg.equals(reg2))

    def test_plink_parsing(self):
        pass

    def test_solve_groups(self):
        pass


if __name__=="__main__":
    unittest.main()