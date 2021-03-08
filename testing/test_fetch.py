import unittest
import unittest.mock as mock
import sys,os
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
import pandas as pd,numpy as np
from Scripts import gws_fetch
from Scripts import autoreporting_utils as autils
from Scripts.data_access.db import LDData, Variant
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

    def test_simple_filtering(self):
        #test simple filtering
        ##NOTE: We can test only filtering by using get_gws_variants instead of fetch_gws
        input_="testing/fetch_resources/filter_test.tsv.gz"
        args=Arg()
        args.gws_fpath=input_
        args.sig_treshold=0.10
        args.columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval","beta":"beta","af":"maf"}
        #output = gws_fetch.fetch_gws(args)
        output = gws_fetch.get_gws_variants(fname=args.gws_fpath, sign_treshold=args.sig_treshold, columns=args.columns)
        validation=pd.read_csv("testing/fetch_resources/filter_test.tsv.gz",compression="gzip",sep="\t")
        validation=validation.loc[validation["pval"]<=args.sig_treshold,:]
        validation=validation.reset_index(drop=True)
        output=output.reset_index(drop=True)
        for col in validation.columns:
            self.assertEqual( list(output[col].astype(str)) , list(validation[col].astype(str)) )

    def test_simple_grouping(self):
        #test simple grouping function
        #set up data and parameters
        input_="testing/fetch_resources/test_grouping.tsv.gz"
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval","beta":"beta","af":"maf","af_cases":"maf_cases","af_controls":"maf_controls"}
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
        val_path="testing/fetch_resources/group_validate_width_1k.tsv"
        validate=pd.read_csv(val_path,sep="\t")
        for col in output.columns:
            self.assertEqual(list(output[col]),list(validate[col]))
        
    def test_ld_grouping(self):
        """
        Test the ld_grouping algorithm
        Test cases: 1: empty dataframe,
                    2: valid dataframe
        Only test PLINK api, as it does not really matter
        (the API should return similar things regadless of the implementation,
        and those should be tested in their own tests)
        """
        input_="testing/fetch_resources/test_grouping.tsv.gz"
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval","beta":"beta","af":"maf","af_cases":"maf_cases","af_controls":"maf_controls"}
        data=pd.read_csv(input_,compression="gzip",sep="\t")
        data["#variant"]=autils.create_variant_column(data)
        data["locus_id"]=data["#variant"]
        data["pos_rmax"]=data["pos"]
        data["pos_rmin"]=data["pos"]
        sig_treshold=0.05
        sig_treshold_2=0.6
        loc_width=1000
        overlap=False
        df_p1=data[data["pval"] <=sig_treshold].copy()
        df_p2=data[data["pval"] <=sig_treshold_2].copy()
        ld_treshold=0.2
        prefix=""
        r2_out = pd.read_csv("testing/fetch_resources/ld_grouping_report.csv",sep="\t")
        r2_out = r2_out.to_dict('records')
        r2_out = [LDData(Variant(
                                 a['chrom_1'],
                                 a['pos_1'],
                                 a['variant_1'].split("_")[2],
                                 a['variant_1'].split("_")[3]),
                         Variant(
                                 a['chrom_2'],
                                 a['pos_2'],
                                 a['variant_2'].split("_")[2],
                                 a['variant_2'].split("_")[3]),
                         a['r2']) for a in r2_out]
        #1
        class AugmentedMock(mock.Mock):
            def get_ranges(*args):
                return []
        emptydf=pd.DataFrame(columns=df_p1.columns)
        emptydf_2=pd.DataFrame(columns=df_p2.columns)
        with mock.patch("Scripts.gws_fetch.PlinkLD",new_callable=AugmentedMock) as ld_api:
            retval=gws_fetch.ld_grouping(emptydf,emptydf_2,
            sig_treshold_2, loc_width, ld_treshold,
            overlap,prefix,ld_api,columns)
        self.assertTrue(retval.empty)
        #2
        class AugmentedMock(mock.Mock):
            def get_ranges(*args):
                return r2_out
        
        with mock.patch("Scripts.gws_fetch.PlinkLD",new_callable=AugmentedMock) as ld_api:
            retval=gws_fetch.ld_grouping(df_p1,df_p2,
            sig_treshold_2, loc_width, ld_treshold,
            overlap,prefix,ld_api,columns)
        retval=retval.sort_values(by=["#variant"]).astype(object).reset_index(drop=True)
        validate=pd.read_csv("testing/fetch_resources/ld_grouping_validate.csv",sep="\t").astype(object).reset_index(drop=True)
        for col in retval.columns:
            self.assertTrue(retval[col].equals(validate[col]))



    def test_get_gws_vars(self):
        #test the get_gws_variants function
        #case 1: empty data, should return empty dataframe.
        empty_read=[pd.DataFrame(columns=["#chrom", "pos", "ref", "alt", "pval"])]
        columns=autils.columns_from_arguments(["#chrom", "pos", "ref", "alt", "pval"])
        extra_columns = ["beta","maf","maf_cases","maf_controls"]
        with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=empty_read):
            retval=gws_fetch.get_gws_variants("",columns=columns)
        self.assertTrue(retval.empty)
        #case 2: valid data and hits
        valid_read = [pd.read_csv("testing/fetch_resources/test_grouping.tsv.gz",sep="\t",compression="gzip")]
        with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=valid_read):
            retval=gws_fetch.get_gws_variants("",sign_treshold=0.4,columns=columns,extra_cols=extra_columns)
        validate=valid_read[0]
        validate=validate[validate["pval"]<0.4].reset_index(drop=True)
        for col in validate.columns:
            self.assertTrue(validate[col].equals(retval[col]))
        #case 3: valid data, but no hits 
        with mock.patch("Scripts.gws_fetch.pd.read_csv",return_value=valid_read):
            retval=gws_fetch.get_gws_variants("",sign_treshold=5e-8,columns=columns,extra_cols=extra_columns)
        validate=valid_read[0]
        validate=validate[validate["pval"]<5e-8].reset_index(drop=True)
        for col in validate.columns:
            self.assertTrue(validate[col].equals(retval[col]))
            
    def test_cred_grouping(self):
        #test credible grouping. Test at least two test cases: with empty credible sets, as well as when using proper data.
        input_ = "testing/fetch_resources/test_grouping.tsv.gz"
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        data=pd.read_csv(input_,compression="gzip",sep="\t")
        data["#variant"]=autils.create_variant_column(data)
        data["locus_id"]=data["#variant"]
        data["pos_rmax"]=data["pos"]
        data["pos_rmin"]=data["pos"]
        data["cs_id"]=np.nan
        data["cs_prob"]=np.nan
        data["r2_to_lead"]=np.nan
        sig_tresh_2 = 0.2
        loc_width=1000
        overlap=False
        ld_treshold=0.2
        prefix=""
        #case 1: empty data
        r2_out=None
        class AugmentedMock(mock.Mock):
            def get_ranges(*args):
                return r2_out
        with mock.patch("Scripts.gws_fetch.PlinkLD",new_callable=AugmentedMock) as ld_api:
                retval = gws_fetch.credible_set_grouping(data ,ld_treshold ,loc_width,overlap,ld_api,columns)
        #return value should be empty, have same columns as data 
        self.assertTrue(retval.empty)
        self.assertEqual(data.columns.all(), retval.columns.all())
        #case 2: not empty data, groups should be formed.
        data.loc[1,"cs_id"] = "chr1_2500_A_C_1"
        data.loc[1,"cs_prob"] = 0.999
        data.loc[9,"cs_id"] = "chrX_1500_A_C_1"
        data.loc[9,"cs_prob"] = 0.997
        r2_out=pd.read_csv("testing/fetch_resources/ld_report.csv",sep="\t")
        r2_out = r2_out.to_dict('records')
        r2_out = [LDData(Variant(
                                 a['chrom_1'],
                                 a['pos_1'],
                                 a['variant_1'].split("_")[2],
                                 a['variant_1'].split("_")[3]),
                         Variant(
                                 a['chrom_2'],
                                 a['pos_2'],
                                 a['variant_2'].split("_")[2],
                                 a['variant_2'].split("_")[3]),
                         a['r2']) for a in r2_out]
        class AugmentedMock(mock.Mock):
            def get_ranges(*args):
                return r2_out
        with mock.patch("Scripts.gws_fetch.PlinkLD",new_callable=AugmentedMock) as ld_api:
            retval = gws_fetch.credible_set_grouping(data,ld_treshold,loc_width,overlap,ld_api,columns)
        #validate
        validate=pd.read_csv("testing/fetch_resources/validate_cred.csv",sep="\t").fillna(-1).sort_values(by=["locus_id","#variant"]).reset_index(drop=True).astype(object)
        retval=retval.astype(dtype={"pos":np.int64,"pos_rmax":np.int64,"pos_rmin":np.int64}).fillna(-1).sort_values(by=["locus_id","#variant"]).reset_index(drop=True).astype(object)
        for col in validate.columns:
            self.assertAlmostEqual(validate[col].all(),retval[col].all())

    def test_cred_merge(self):
        # test merging of data
        # mock: load_pysam_df
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
        with mock.patch("Scripts.gws_fetch.load_pysam_df",return_value=load_retval) as mock_load_df:
            outdf = gws_fetch.merge_credset(df,cs_df,fname,columns)
        self.assertTrue(outdf.equals(validate))
        pass
        
if __name__=="__main__":
    unittest.main()