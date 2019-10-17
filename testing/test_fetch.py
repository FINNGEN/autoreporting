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
        args.cred_set_file=""
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
        for col in validation.columns:
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

    def test_get_gws_vars(self):
        #test the get_gws_variants function
        #case 1: empty data, should return empty dataframe.
        empty_read=[pd.DataFrame(columns=["#chrom", "pos", "ref", "alt", "pval"])]
        with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=empty_read):
            retval=gws_fetch.get_gws_variants("")
        self.assertTrue(retval.empty)
        #case 2: valid data and hits
        valid_read = [pd.read_csv("fetch_resources/test_grouping.tsv.gz",sep="\t",compression="gzip")]
        with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=valid_read):
            retval=gws_fetch.get_gws_variants("",sign_treshold=0.4)
        validate=valid_read[0]
        validate=validate[validate["pval"]<0.4].reset_index(drop=True)
        for col in validate.columns:
            self.assertTrue(validate[col].equals(retval[col]))
        #case 3: valid data, but no hits 
        with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=valid_read):
            retval=gws_fetch.get_gws_variants("",sign_treshold=5e-8)
        validate=valid_read[0]
        validate=validate[validate["pval"]<5e-8].reset_index(drop=True)
        for col in validate.columns:
            self.assertTrue(validate[col].equals(retval[col]))
            
    @mock.patch('Scripts.gws_fetch.subprocess.Popen',side_effect=return_sp_mock)
    @mock.patch('Scripts.gws_fetch.subprocess.call')
    def test_cred_grouping(self,mocked_subprocess, mocked_call):
        #test credible grouping. Test at least two test cases: with empty credible sets, as well as when using proper data.
        #case 1: empty data
        input_ = "fetch_resources/test_grouping.tsv.gz"
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        data=pd.read_csv(input_,compression="gzip",sep="\t")
        data["#variant"]=autils.create_variant_column(data)
        data["locus_id"]=data["#variant"]
        data["pos_rmax"]=data["pos"]
        data["pos_rmin"]=data["pos"]
        data["cs_id"]=np.nan
        data["cs_prob"]=np.nan
        sig_tresh_2 = 0.2
        loc_width=1000
        overlap=False
        ld_treshold=0.2
        ld_panel_path=""
        plink_mem=12000
        prefix=""
        r2_out=None
        _open=mock.mock_open()
        with mock.patch("Scripts.gws_fetch.pd.DataFrame.to_csv"):
            with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=r2_out):
                with mock.patch("Scripts.gws_fetch.open",_open):
                    retval = gws_fetch.credible_set_grouping(data,sig_tresh_2,ld_panel_path,ld_treshold,loc_width,overlap,columns)
        #return value should be empty, have same columns as data 
        self.assertTrue(retval.empty)
        self.assertEqual(data.columns.all(), retval.columns.all())
        #case 2: not empty data, groups should be formed.
        data.loc[1,"cs_id"] = "chr1_2500_A_C_1"
        data.loc[1,"cs_prob"] = 0.999
        data.loc[9,"cs_id"] = "chrX_1500_A_C_1"
        data.loc[9,"cs_prob"] = 0.997
        r2_out=pd.read_csv("fetch_resources/ld_report.csv",sep="\t")
        with mock.patch("Scripts.gws_fetch.pd.DataFrame.to_csv"):
            with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=r2_out):
                with mock.patch("Scripts.gws_fetch.open",_open):
                    retval = gws_fetch.credible_set_grouping(data,sig_tresh_2,ld_panel_path,ld_treshold,loc_width,overlap,columns)
        #validate
        validate=pd.read_csv("fetch_resources/validate_cred.csv",sep="\t").fillna(-1)
        retval=retval.astype(dtype={"pos":np.int64,"pos_rmax":np.int64,"pos_rmin":np.int64}).fillna(-1)
        for col in validate.columns:
            self.assertAlmostEqual(validate[col].all(),retval[col].all())

    def test_cred_merge(self):
        # test merging of data
        # mock: tabix, load_tb_df
        columns={"chrom":"chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}

        df={"chrom":["1","1","1","1"],"pos":[1,2,3,4],"ref":["A","G","C","T"],"alt":["G","C","T","A"],"pval":[0.01,0.001,0.0001,0.00001]}
        df=pd.DataFrame(df)
        cs_df={"chrom":["1","1"],"pos":[5,6],"ref":["G","C"],"alt":["T","A"],"cs_id":["CS_ID_1","CS_ID_2"],"cs_prob":[0.9,0.8]}
        cs_df=pd.DataFrame(cs_df)
        load_retval ={"chrom":["1","1"],"pos":[5,6],"ref":["G","C"],"alt":["T","A"],"pval":[0.1,0.2]}
        load_retval=pd.DataFrame(load_retval)
        validate={"chrom":["1","1","1","1","1","1"],"pos":[1,2,3,4,5,6],"ref":["A","G","C","T","G","C"],
        "alt":["G","C","T","A","T","A"],"pval":[0.01,0.001,0.0001,0.00001,0.1,0.2],"cs_id":[np.nan,np.nan,np.nan,np.nan,"CS_ID_1","CS_ID_2"],"cs_prob":[np.nan,np.nan,np.nan,np.nan,0.9,0.8]}
        validate=pd.DataFrame(validate)
        fname="test_fname"
        with mock.patch("Scripts.gws_fetch.tabix.open") as mock_tabix:
            with mock.patch("Scripts.gws_fetch.load_tb_df",return_value=load_retval) as mock_load_df:
                outdf = gws_fetch.merge_credset(df,cs_df,fname,columns)
        self.assertTrue(mock_tabix.called_once_with(fname))
        self.assertTrue(outdf.equals(validate))
        pass
        
if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()