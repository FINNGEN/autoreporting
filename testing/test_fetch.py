import unittest
import sys,os
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
import pandas as pd,numpy as np
from Scripts import gws_fetch


class TestGws(unittest.TestCase):

    def test_plink_parsing(self):
        #load data
        df=pd.read_csv("fetch_resources/test_plink.clumped",sep="\s+")
        df=df.loc[:,["SNP","TOTAL","SP2"]]
        #create end result
        df2={"#chrom":["1","1","X","X","X"],
            "pos":["111111111","111111112","111111113","111111114","111111115"],
            "ref":["A","C","A","G","C"],
            "alt":["T"]*5}
        df2=pd.DataFrame(df2)
        df2.loc[:,"#variant"]=gws_fetch.create_variant_column(df2)
        #test parsing
        df=gws_fetch.parse_plink_output(df)
        self.assertTrue(df2.equals(df))

    def test_solve_groups(self):
        #create parameters for function
        result_header=["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "beta", "sebeta", "maf", "maf_cases", "maf_controls","#variant","locus_id"]
        df=pd.DataFrame(columns=result_header)
        group_df=pd.read_csv("fetch_resources/test_plink.clumped",sep="\s+")
        group_df=group_df.loc[:,["SNP","TOTAL","SP2"]]
        tabixdf=pd.read_csv("fetch_resources/test_tsv",sep="\t")
        tabixdf.loc[:,"#variant"]=gws_fetch.create_variant_column(tabixdf)
        df2=pd.read_csv("fetch_resources/solve_groups_result.tsv",sep="\t")
        df=gws_fetch.solve_groups(df,group_df,tabixdf)
        df.loc[:,"pos"]=df.loc[:,"pos"].astype(np.int64)
        #test for equality
        self.assertTrue(df.equals(df2))
        pass


if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()