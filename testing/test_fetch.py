import unittest
import sys,os
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
import pandas as pd,numpy as np
from Scripts import gws_fetch
from Scripts import autoreporting_utils as autils
from io import StringIO

class Arg():
    def __init__(self):
        pass

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

    def test_group_range(self):
        #input possibilities: df with no entries, df with one entry, df with multiple entries of one group
        #column names: since they are specified in the function call, no need to test them
        # Case 1: No entries in df
        # Function should return None
        df=pd.DataFrame(columns=["#chrom","pos","ref","alt","locus_id"])
        gr_var="chr1_100_A_T"
        retval=gws_fetch.get_group_range(df,gr_var,{"pos":"pos"})
        self.assertEqual(retval, None)
        # Case 2: df with one entry, matching gr_var
        # Function should return {"min":pos,"max":pos}
        pos=12345
        df=pd.DataFrame([["1",pos,"A","T"]],columns=["#chrom","pos","ref","alt"])
        df.loc[:,"locus_id"]=autils.create_variant_column(df)
        gr_var="chr1_12345_A_T"
        retval=gws_fetch.get_group_range(df,gr_var,{"pos":"pos"})
        self.assertTrue(retval["min"]==pos & retval["max"] == pos)
        # Case 3: df with one entry, no matching gr_var
        #function should return None
        pos=12345
        df=pd.DataFrame([["1",pos,"A","T"]],columns=["#chrom","pos","ref","alt"])
        df.loc[:,"locus_id"]=autils.create_variant_column(df)
        gr_var="chr1_11111_G_C"
        retval=gws_fetch.get_group_range(df,gr_var,{"pos":"pos"})
        self.assertEqual(retval, None)
        # Case 4: df with multiple entries, matching gr_var
        #function should return group range min and max
        df=pd.read_csv("fetch_resources/group_range_df.tsv",sep="\t",dtype={"pos":np.int32,"#chrom":str})
        gr_var="chr1_1111_A_C"
        min_=1111
        max_=4321
        retval=gws_fetch.get_group_range(df, gr_var, columns={"pos":"pos"})
        self.assertEqual(min_,retval["min"])
        self.assertEqual(max_,retval["max"])
        pass

    def test_simple_filtering(self):
        #test simple filtering
        input_="fetch_resources/filter_test.tsv.gz"
        args=Arg()
        args.gws_fpath=input_
        args.sig_treshold=0.20
        args.sig_treshold_2=0.20
        args.loc_width=1
        args.grouping=False
        args.column_labels=["#chrom","pos","ref","alt","pval"]
        args.fetch_out=StringIO()
        gws_fetch.fetch_gws(args)
        args.fetch_out.seek(0)
        output=pd.read_csv(args.fetch_out,sep="\t")
        validation=pd.read_csv("fetch_resources/filter_test.tsv.gz",compression="gzip",sep="\t")
        validation=validation.loc[validation["pval"]<=args.sig_treshold,:]
        validation["pos_rmin"]=validation["pos"]
        validation["pos_rmax"]=validation["pos"]
        validation["#variant"]=autils.create_variant_column(validation)
        validation["locus_id"]=autils.create_variant_column(validation)
        validation=validation.reset_index(drop=True).sort_values(by="#variant")
        output=output.reset_index(drop=True)
        for col in output.columns:
            self.assertEqual( list(output[col]) , list(validation[col]) )

    def test_simple_grouping(self):
        #test simple grouping 
        input_="fetch_resources/test_grouping.tsv.gz"
        args=Arg()
        args.gws_fpath=input_
        args.sig_treshold=0.05
        args.sig_treshold_2=0.05
        args.column_labels=["#chrom","pos","ref","alt","pval"]
        args.loc_width=1
        args.grouping=True
        args.grouping_method="simple"
        args.overlap=False
        args.fetch_out=StringIO()

        gws_fetch.fetch_gws(args)
        args.fetch_out.seek(0)
        output=pd.read_csv(args.fetch_out,sep="\t")
        val_path="fetch_resources/group_validate_width_1k.tsv"
        validate=pd.read_csv(val_path,sep="\t")
        for col in output.columns:
            self.assertEqual(list(output[col]),list(validate[col]))

if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()