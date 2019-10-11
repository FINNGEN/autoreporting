import unittest
import unittest.mock as mock
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

def return_sp_mock(arg,stdout, stderr=None,encoding=None,*args):
    retval=mock.Mock()
    retval.returncode=0
    retval.wait=lambda:True
    retval.stdout.readlines=lambda:[""]
    return retval

class TestGws(unittest.TestCase):

    def test_solve_groups(self):
        #create parameters for function
        result_header=["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "beta", "sebeta", "maf", "maf_cases", "maf_controls","#variant","locus_id"]
        df=pd.DataFrame(columns=result_header)
        group_df=pd.read_csv("fetch_resources/test_plink.clumped",sep="\s+")
        group_df=group_df.loc[:,["SNP","TOTAL","SP2"]]
        df=pd.read_csv("fetch_resources/test_tsv",sep="\t")
        df.loc[:,"#variant"]=gws_fetch.create_variant_column(df)
        df2=pd.read_csv("fetch_resources/solve_groups_result.tsv",sep="\t")
        df=gws_fetch.solve_groups(df,group_df)
        df.loc[:,"pos"]=df.loc[:,"pos"].astype(np.int64)
        #test for equality
        for col in df:
            self.assertTrue(df[col].equals(df2[col]))
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
        args.ignore_region=""
        args.prefix=""
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
        #test simple grouping function
        #set up data and parameters
        input_="fetch_resources/test_grouping.tsv.gz"
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        data=pd.read_csv(input_,compression="gzip",sep="\t")
        data["#variant"]=autils.create_variant_column(data)
        data["locus_id"]=data["#variant"]
        data["pos_rmax"]=data["pos"]
        data["pos_rmin"]=data["pos"]
        sig_treshold=0.05
        sig_treshold_2=0.05
        loc_width=1000
        overlap=False
        df_p1=data[data["pval"] <=sig_treshold]
        df_p2=data[data["pval"] <=sig_treshold_2]
        #run function
        output=gws_fetch.simple_grouping(df_p1,df_p2,loc_width,overlap,columns)
        output=output.sort_values(["locus_id","#variant"])
        val_path="fetch_resources/group_validate_width_1k.tsv"
        validate=pd.read_csv(val_path,sep="\t")
        for col in output.columns:
            self.assertEqual(list(output[col]),list(validate[col]))

    #@mock.patch('')
    @mock.patch('Scripts.gws_fetch.subprocess.Popen',side_effect=return_sp_mock)
    @mock.patch('Scripts.gws_fetch.subprocess.call')
    def test_ld_grouping(self,mocked_subprocess,mocked_call):
        mocked_subprocess.return_value=0
        mocked_subprocess.wait=lambda:True
        input_="fetch_resources/test_grouping.tsv.gz"
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        data=pd.read_csv(input_,compression="gzip",sep="\t")
        data["#variant"]=autils.create_variant_column(data)
        data["locus_id"]=data["#variant"]
        data["pos_rmax"]=data["pos"]
        data["pos_rmin"]=data["pos"]
        sig_treshold=0.05
        sig_treshold_2=0.05
        loc_width=1000
        overlap=False
        df_p1=data[data["pval"] <=sig_treshold]
        df_p2=data[data["pval"] <=sig_treshold_2]
        ld_treshold=0.2
        ld_panel_path=""
        plink_mem=12000
        prefix=""
        ld_out=pd.read_csv("fetch_resources/plink_out.csv",sep="\t")
        with mock.patch("Scripts.gws_fetch.pd.DataFrame.to_csv") as mock_write:
            with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=ld_out) as mock_read:
                with mock.patch("Scripts.gws_fetch.os.path.exists") as mock_exists:
                    retval=gws_fetch.ld_grouping(df_p1,df_p2,sig_treshold,sig_treshold_2,loc_width,ld_treshold,ld_panel_path,plink_mem,overlap,prefix,columns)
        validate=pd.read_csv("fetch_resources/ld_group_validate.csv",sep="\t")
        retval=retval.astype(dtype={"pos":np.int64,"pos_rmax":np.int64,"pos_rmin":np.int64})
        retval=retval.reset_index(drop=True)
        for col in validate.columns:
            self.assertTrue(validate[col].equals(retval[col]))
if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()