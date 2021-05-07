import unittest
import sys,os
sys.path.append("../")
sys.path.append("./")
import pandas as pd,numpy as np
from Scripts import autoreporting_utils
from Scripts.autoreporting_utils import Region

class TestUtils(unittest.TestCase):
    def test_create_variant_column(self):
        df_data={"#chrom":"1","pos":1234,"ref":"CGC","alt":"TTC"}
        df=pd.DataFrame(df_data,index=[0])
        variant=pd.Series("chr{}_{}_{}_{}".format("1",1234,"CGC","TTC"))
        variant_2=autoreporting_utils.create_variant_column(df)
        value=(variant==variant_2).all()
        self.assertTrue(variant.equals(variant_2))

    def test_region_pruning(self):
        #case no regions
        reg=[]
        out=[]
        reg2=autoreporting_utils.prune_regions(reg)
        self.assertEqual(reg,reg2)
        #case one region
        reg=[Region("1",10,20)]
        reg2=autoreporting_utils.prune_regions(reg)
        self.assertEqual(reg,reg2)
        #case two regions: no overlap by position
        regs = [Region("1",10,20),Region("1",30,40)]
        reg2=autoreporting_utils.prune_regions(regs)
        self.assertEqual(regs,reg2)
        #case two regions: overlapping positions, different chromosomes
        regs = [Region("1",10,40),Region("2",20,30)]
        reg2=autoreporting_utils.prune_regions(regs)
        self.assertEqual(regs,reg2)
        #case two regions: one in other, i.e. 1_min 2_min 2_max 1_max
        validate=[Region("1",10,40)]
        regs=[Region("1",10,40),Region("1",20,30)]
        reg2=autoreporting_utils.prune_regions(regs)
        self.assertEqual(validate,reg2)
        #case two regions: 1_min 2_min 1_max 2_max
        validate = [Region("1",10,40)]
        regs=[Region("1",10,30), Region("1",15,40)]
        reg2=autoreporting_utils.prune_regions(regs)
        self.assertEqual(validate,reg2)
        #case two regions: 2_min 1_min 2_max 1_max
        validate = [Region("1",10,40)]
        regs = [Region("1",15,40), Region("1",10,30)]
        reg2=autoreporting_utils.prune_regions(regs)
        self.assertEqual(reg2, validate)
        #case two regions: 2_min 1_min 1_max 2_max
        validate = [Region("1",10,40)]
        regs=[Region("1",15,30), Region("1",10,40)]
        df=pd.DataFrame({"#chrom":["1","1"],"pos_rmax":[30,40],"pos_rmin":[15,10]})
        reg2=autoreporting_utils.prune_regions(regs)
        self.assertEqual(reg2, validate)
        #6 regions, with 3 output regions
        regions=[
            Region("1",1,100),
            Region("1",75,150),
            Region("1",100,200),
            Region("1",250,251),
            Region("1",251,260),
            Region("1",300,305),
        ]
        validation=[
            Region("1",1,200),
            Region("1",250,260),
            Region("1",300,305)
        ]
        out=autoreporting_utils.prune_regions(regions)
        self.assertEqual(set(validation),set(out))

    #TODO: test get gzip header

    def test_column_labels(self):
        collabs = ["chrom","pos","ref","alt","pval"]
        collabs2 = ["1","2","3","4","5"]
        columns= autoreporting_utils.columns_from_arguments(collabs)
        columns2= autoreporting_utils.columns_from_arguments(collabs2)
        validate= {"chrom":"chrom","pos":"pos", "ref":"ref", "alt":"alt", "pval":"pval"}
        validate2= {"chrom":"1","pos":"2", "ref":"3", "alt":"4", "pval":"5"}
        self.assertEqual(columns, validate)
        self.assertEqual(columns2, validate2)

    def test_pysam_df(self):
        import gzip
        annotate_fpath = "testing/annotate_resources/gnomad_genomes.tsv.gz"
        data = {"chr":["1","1"],"position":[1,10],"otherdata":[1,2]}
        data=pd.DataFrame(data)
        columns = {
            "chrom":"chr",
            "pos":"position",
            "ref":"ref",
            "alt":"alt",
            "pval":"pval"
        }
        #validation data
        with gzip.open(annotate_fpath,"rt") as f:
            vlines= f.readlines()
            vheader = vlines[0].strip('\n').split("\t")
            vdata = [a.strip('\n').split("\t") for a in vlines[1:]]
            vdata = pd.DataFrame(vdata,columns=vheader)
            vdata[vdata.columns]=vdata[vdata.columns].apply(pd.to_numeric,errors="ignore")
        annotation_data = autoreporting_utils.load_pysam_df(data,annotate_fpath,columns)
        self.assertTrue(vdata.equals(annotation_data))
if __name__=="__main__":
    unittest.main()