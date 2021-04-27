import unittest
import unittest.mock as mock
import sys,os
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
import pandas as pd,numpy as np
from Scripts import gws_fetch
from Scripts import autoreporting_utils as autils
from Scripts.data_access.db import LDData, Variant, LDAccess
from io import StringIO
from typing import List, Dict 

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
        

def create_loci(variants:List[str])->pd.DataFrame:
    cols = [
        "chrom1",
        "pos1",
        "ref1",
        "alt1",
        "pval1",
        "#variant",
        "cs_prob",
        "cs_min_r2",
        "cs_log10bf",
        "good_cs",
        "cs_region",
        "cs_size",
        "cs_id",
        "r2_to_lead"
    ]
    data = [
        [
            a.split(":")[0],
            int(a.split(":")[1]),
            "A",
            "T",
            float(a.split(":")[2]),
            f'chr{a.split(":")[0]}_{ a.split(":")[1] }_A_T',
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan
        ]
        for a in variants
    ]
    return (data,cols)

def create_cs_loci(leads:List[str],variants:Dict[str,str])->pd.DataFrame:
    cols = [
        "chrom1",
        "pos1",
        "ref1",
        "alt1",
        "pval1",
        "#variant",
        "cs_prob",
        "cs_min_r2",
        "cs_log10bf",
        "good_cs",
        "cs_region",
        "cs_size",
        "cs_id",
        "r2_to_lead"
    ]
    lead_data = [
        [
            a.split(":")[0],
            int(a.split(":")[1]),
            "A",
            "T",
            float(a.split(":")[2]),
            f'chr{a.split(":")[0]}_{ a.split(":")[1] }_A_T',
            0.95,
            0.5,
            5.31,
            True,
            f'{a.split(":")[0]}:{int(a.split(":")[1])-1500000}-{int(a.split(":")[1])+1500000}',
            len(variants[a])+1,
            f'chr{a.split(":")[0]}_{a.split(":")[1]}_A_T_1',
            1.0
        ]
        for a in leads
    ]

    var_data = [
        [
            c.split(":")[0],
            int(c.split(":")[1]),
            "A",
            "T",
            5e-6,
            f'chr{c.split(":")[0]}_{int(c.split(":")[1])}_A_T',
            0.3,
            0.5,
            5.31,
            True,
            f'{a.split(":")[0]}:{int(a.split(":")[1])-1500000}-{int(a.split(":")[1])+1500000}',
            len(variants[a])+1,
            f'chr{a.split(":")[0]}_{int(a.split(":")[1])}_A_T_1',
            0.85
        ]
        for (a,b) in variants.items()
        for c in b
    ]
    data = lead_data+var_data
    return (data,cols)#pd.DataFrame(data,columns=cols)

class PosLD(LDAccess):
    def __init__(self, data):
        self.data=data

    def get_range(self, variant, bp_range,ld_threshold):
        #filter by bp range and chrom
        vs = [a for a in self.data if (abs(a.pos - variant.pos)<bp_range)]
        vs = [ a for a in vs if (a.chrom == variant.chrom)]
        #create LD going from 1 to 0 on a space of 500kb
        ld = [
            LDData(
                variant,
                a,
                float(max(1- abs(variant.pos-a.pos)/500000,0))
            )
            for a in vs
        ]
        ld  = [a for a in ld if a.r2>ld_threshold ]
        return ld
        

class TestCredGrouping(unittest.TestCase):
    def test_credible_set_grouping(self):
        """Test credible set grouping
            Say, some groups. Let's say that for LD partners, our simulation data's r2 goes linearly to zero with distance.
            LD goes linearly fron 1 with 0 distance to 0 with 500_000 distance.
            LD calculations are delegated to a mock object that holds the LD data in a dict and returns it.
        """
        cols = {
            "chrom":"chrom1",
            "pos":"pos1",
            "ref":"ref1",
            "alt":"alt1",
            "pval":"pval1"
        }
        loci = [
            "1:2000000:1e-8",
            "1:5000000:5e-9",
            "2:34000000:2e-10"
        ]
        cs_vars = {
            "1:2000000:1e-8":[
                "1:1900000",
                "1:1950000",
                "1:2100000",
                "1:2120000"
            ],
            "1:5000000:5e-9":[
                "1:4750000",
                "1:4500000",
                "1:5500000"
            ],
            "2:34000000:2e-10":[
                "2:33000000",
                "2:35000000"
            ]
        }

        cred_data,c=create_cs_loci(loci,cs_vars)

        ld_vars = [
            "1:1850000:1e-5",
            "1:1960000:1e-5",
            "1:2050000:1e-5",
            "1:4900000:1e-5",
            "1:5250000:1e-5",
            "1:5150000:1e-5",
            "2:33950000:1e-5",
            "2:34050000:1e-5",
            "2:34100000:1e-5"
        ]
        ld_data,c = create_loci(ld_vars)
        list_data = cred_data+ld_data
        data = pd.DataFrame(list_data,columns=c)
        data["locus_id"]=np.nan

        ld_variants = [
            Variant(
                a.split(":")[0],
                int(a.split(":")[1]),
                "A",
                "T"
            )
            for a in ld_vars
        ]
        #run incremental credible grouping
        ld_api = PosLD(ld_variants)
        output = gws_fetch.credible_grouping(data, dynamic_r2=False, ld_threshold=0.2, locus_range=1500000, overlap=False, ld_api=ld_api, columns=cols) 
        validation = [
            ["chr1_1850000_A_T", "chr1_2000000_A_T", np.nan],
            ["chr1_1900000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_1950000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_1960000_A_T", "chr1_2000000_A_T", np.nan],
            ["chr1_2000000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_2050000_A_T", "chr1_2000000_A_T", np.nan],
            ["chr1_2100000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_2120000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_4500000_A_T", "chr1_5000000_A_T", "chr1_5000000_A_T_1"],
            ["chr1_4750000_A_T", "chr1_5000000_A_T", "chr1_5000000_A_T_1"],
            ["chr1_4900000_A_T", "chr1_5000000_A_T", np.nan],
            ["chr1_5000000_A_T", "chr1_5000000_A_T", "chr1_5000000_A_T_1"],
            ["chr1_5150000_A_T", "chr1_5000000_A_T", np.nan],
            ["chr1_5250000_A_T", "chr1_5000000_A_T", np.nan],
            ["chr1_5500000_A_T", "chr1_5000000_A_T", "chr1_5000000_A_T_1"],
            ["chr2_33000000_A_T", "chr2_34000000_A_T", "chr2_34000000_A_T_1"],
            ["chr2_33950000_A_T", "chr2_34000000_A_T", np.nan],
            ["chr2_34000000_A_T", "chr2_34000000_A_T", "chr2_34000000_A_T_1"],
            ["chr2_34050000_A_T", "chr2_34000000_A_T", np.nan],
            ["chr2_34100000_A_T", "chr2_34000000_A_T", np.nan],
            ["chr2_35000000_A_T", "chr2_34000000_A_T", "chr2_34000000_A_T_1"]
        ]
        #only test for grouping now, leaving r2 values etc out of this.
        validation_df = pd.DataFrame(validation,columns=["#variant","locus_id","cs_id"])
        test_out = output[["#variant","locus_id","cs_id"]].sort_values("#variant").reset_index(drop=True)
        if not validation_df.equals(test_out):
            msg = f"Output does not equal validation data!\n validation:\n{validation_df}\noutput:\n{test_out}"
            raise Exception(msg)
        
    def test_overlapping_cs(self):
        """Check that multiple overlapping cs do not "steal" each other's CS variants
        """
        cols = {
            "chrom":"chrom1",
            "pos":"pos1",
            "ref":"ref1",
            "alt":"alt1",
            "pval":"pval1"
        }
        loci = [
            "1:2000000:1e-8",
            "1:2500000:5e-9"
        ]
        cs_vars = {
            "1:2000000:1e-8":[
                "1:1900000",
                "1:1950000",
                "1:2100000",
                "1:2120000"
            ],
            "1:2500000:5e-9":[
                "1:2200000",
                "1:2300000",
                "1:2400000"
            ]
        }
        cred_data,c=create_cs_loci(loci,cs_vars)
        cred_variants = [
            Variant(
                a[0],
                a[1],
                "A",
                "T"
            )
            for a in cred_data
        ]
        data = pd.DataFrame(cred_data,columns=c)
        data["locus_id"]=np.nan
        ld_api = PosLD(cred_variants)
        
        output = gws_fetch.credible_grouping(data, dynamic_r2=False, ld_threshold=0.00, locus_range=1500000, overlap=False, ld_api=ld_api, columns=cols) 
        out = output[["#variant","locus_id","cs_id"]].sort_values(["#variant","locus_id"]).reset_index(drop=True)
        validation= [
            ["chr1_1900000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_1950000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_2000000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_2100000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_2100000_A_T", "chr1_2500000_A_T", np.nan],
            ["chr1_2120000_A_T", "chr1_2000000_A_T", "chr1_2000000_A_T_1"],
            ["chr1_2120000_A_T", "chr1_2500000_A_T", np.nan],
            ["chr1_2200000_A_T", "chr1_2500000_A_T", "chr1_2500000_A_T_1"],
            ["chr1_2300000_A_T", "chr1_2500000_A_T", "chr1_2500000_A_T_1"],
            ["chr1_2400000_A_T", "chr1_2500000_A_T", "chr1_2500000_A_T_1"],
            ["chr1_2500000_A_T", "chr1_2500000_A_T", "chr1_2500000_A_T_1"],
        ]
        validation_df = pd.DataFrame(validation,columns=["#variant","locus_id","cs_id"])
        if not out.equals(validation_df):
            msg = f"Output does not equal validation data!\n validation:\n{validation_df}\noutput:\n{out}"
            raise Exception(msg)


class TestLDGrouping(unittest.TestCase):
    def test_ld_grouping(self):
        """Test non-overlapping LD set grouping
            3 groups with 2 groups being in same chromosome, but in completely separate positions.
            LD goes linearly fron 1 with 0 distance to 0 with 500_000 distance.
            LD from "PosLD" that gives an r2 value based on pos difference.
        """
        cols = {
            "chrom":"chrom1",
            "pos":"pos1",
            "ref":"ref1",
            "alt":"alt1",
            "pval":"pval1"
        }
        loci = [
            "1:2000000:1e-8",
            "1:5000000:5e-9",
            "2:34000000:2e-10"
        ]
        cs_vars = {
            "1:2000000:1e-8":[
                "1:1900000",
                "1:1950000",
                "1:2100000",
                "1:2120000"
            ],
            "1:5000000:5e-9":[
                "1:4750000",
                "1:4500000",
                "1:5500000"
            ],
            "2:34000000:2e-10":[
                "2:33600000",
                "2:34500000"
            ]
        }

        cred_data,c=create_cs_loci(loci,cs_vars)

        ld_vars = [
            "1:1850000:1e-5",
            "1:1960000:1e-5",
            "1:2050000:1e-5",
            "1:4900000:1e-5",
            "1:5250000:1e-5",
            "1:5150000:1e-5",
            "2:33950000:1e-5",
            "2:34050000:1e-5",
            "2:34100000:1e-5"
        ]
        ld_data,c = create_loci(ld_vars)
        list_data = cred_data+ld_data
        data = pd.DataFrame(list_data,columns=c)
        data["locus_id"]=np.nan
        df_p1 = data[data["pval1"]<1e-6].copy()
        df_p2 = data.copy()
        ld_variants = [
            Variant(
                a[0],
                a[1],
                "A",
                "T"
            )
            for a in list_data
        ]
        ld_api = PosLD(ld_variants)
        #run incremental credible grouping
        output = gws_fetch.ld_grouping(df_p1,df_p2,sig_threshold_2 = 1e-6, dynamic_r2=False, ld_threshold=0.15, locus_width=1500000, overlap=False, ld_api=ld_api, columns=cols)
        out = output[["#variant","locus_id"]].sort_values(["#variant","locus_id"]).reset_index(drop=True) 
        validation= [
            ["chr1_1850000_A_T",   "chr1_2000000_A_T"],
            ["chr1_1900000_A_T",   "chr1_2000000_A_T"],
            ["chr1_1950000_A_T",   "chr1_2000000_A_T"],
            ["chr1_1960000_A_T",   "chr1_2000000_A_T"],
            ["chr1_2000000_A_T",   "chr1_2000000_A_T"],
            ["chr1_2050000_A_T",   "chr1_2000000_A_T"],
            ["chr1_2100000_A_T",   "chr1_2000000_A_T"],
            ["chr1_2120000_A_T",   "chr1_2000000_A_T"],
            ["chr1_4750000_A_T",   "chr1_5000000_A_T"],
            ["chr1_4900000_A_T",   "chr1_5000000_A_T"],
            ["chr1_5000000_A_T",   "chr1_5000000_A_T"],
            ["chr1_5150000_A_T",   "chr1_5000000_A_T"],
            ["chr1_5250000_A_T",   "chr1_5000000_A_T"],
            ["chr2_33600000_A_T",  "chr2_34000000_A_T"],
            ["chr2_33950000_A_T",  "chr2_34000000_A_T"],
            ["chr2_34000000_A_T",  "chr2_34000000_A_T"],
            ["chr2_34050000_A_T",  "chr2_34000000_A_T"],
            ["chr2_34100000_A_T",  "chr2_34000000_A_T"]
        ]
        validation_df = pd.DataFrame(validation,columns=["#variant","locus_id"])
        if not out.equals(validation_df):
            msg = f"Output does not equal validation data!\n validation:\n{validation_df}\noutput:\n{out}"
            raise Exception(msg)


    def test_overlapping_groups(self):
        """Test whether the grouping functions properly for overlapping signals.
        Say that there are two "signals" that are overlapping:
        peak 1 with pval 1e-10 on 10Mb
        peak 2 with pval 1e-9 on 10.5Mb
        r2 threshold 0.2 -> variants inside 400 kb are taken to group
        all ld variants in between, so the first peak should steal quite a bit of them from the peak 2
        """
        cols = {
            "chrom":"chrom1",
            "pos":"pos1",
            "ref":"ref1",
            "alt":"alt1",
            "pval":"pval1"
        }
        variants = [
            "1:10000000:1e-10",#peak 1
            "1:10500000:1e-9",#peak 2
            "1:10100000:1e-5",
            "1:10200000:1e-5",
            "1:10300000:1e-5",
            "1:10400000:1e-5",
            "1:10150000:1e-5",
            "1:10250000:1e-5",
            "1:10350000:1e-5",
            "1:10450000:1e-5"
        ]
        data,c = create_loci(variants)
        ld_variants = [
            Variant(
                a[0],
                a[1],
                "A",
                "T"
            )
            for a in data
        ]
        df = pd.DataFrame(data,columns=c)
        df["locus_id"]=np.nan
        ld_api = PosLD(ld_variants)
        df_p1 = df[df["pval1"]<1e-6].copy()
        df_p2 = df.copy()
        output = gws_fetch.ld_grouping(df_p1,df_p2,sig_threshold_2 = 1e-2, dynamic_r2=False, ld_threshold=0.2, locus_width=1000000, overlap=False, ld_api=ld_api, columns=cols)

        out = output[["#variant","locus_id"]].sort_values(["#variant"])

        validation = [
            ["chr1_10000000_A_T",  "chr1_10000000_A_T"],
            ["chr1_10100000_A_T",  "chr1_10000000_A_T"],
            ["chr1_10200000_A_T",  "chr1_10000000_A_T"],
            ["chr1_10300000_A_T",  "chr1_10000000_A_T"],
            ["chr1_10150000_A_T",  "chr1_10000000_A_T"],
            ["chr1_10250000_A_T",  "chr1_10000000_A_T"],
            ["chr1_10350000_A_T",  "chr1_10000000_A_T"],
            ["chr1_10500000_A_T",  "chr1_10500000_A_T"],
            ["chr1_10400000_A_T",  "chr1_10500000_A_T"],
            ["chr1_10450000_A_T",  "chr1_10500000_A_T"]
        ]
        validation_df =  pd.DataFrame(validation,columns=["#variant","locus_id"]).sort_values(["#variant"])

        if not out.equals(validation_df):
            msg = f"Output does not equal validation data!\n validation:\n{validation_df}\noutput:\n{out}"
            raise Exception(msg)

    def test_static_r2(self):
        """Test that a static r2 threshold works as expected
        """
        cols = {
            "chrom":"chrom1",
            "pos":"pos1",
            "ref":"ref1",
            "alt":"alt1",
            "pval":"pval1"
        }
        variants = [
            "1:10000000:1e-10",#peak 1
            "1:10100000:1e-5",
            "1:10200000:1e-5",
            "1:10300000:1e-5",
            "1:10400000:1e-5",
            "1:10500000:1e-5",
            "1:10150000:1e-5",
            "1:10250000:1e-5",
            "1:10350000:1e-5",
            "1:10450000:1e-5"
        ]
        data,c = create_loci(variants)
        ld_variants = [
            Variant(
                a[0],
                a[1],
                "A",
                "T"
            )
            for a in data
        ]
        r2_threshold = 0.201
        df = pd.DataFrame(data,columns=c)
        df["locus_id"]=np.nan
        ld_api = PosLD(ld_variants)
        df_p1 = df[df["pval1"]<1e-6].copy()
        df_p2 = df.copy()
        # r2 threshold is 0.201, which means that everything farther away than 399 kb or so is rejected.
        # Therefore there should be 7 variants in the group incl. lead
        output = gws_fetch.ld_grouping(
            df_p1,
            df_p2,
            sig_threshold_2 = 1e-2,
            dynamic_r2=False,
            ld_threshold=r2_threshold,
            locus_width=1000000,
            overlap=False,
            ld_api=ld_api,
            columns=cols
        )
        out=output[["#variant","locus_id","r2_to_lead"]].sort_values(["#variant"]).reset_index(drop=True)
        validation = [
            ["chr1_10000000_A_T", "chr1_10000000_A_T",1.0],
            ["chr1_10100000_A_T", "chr1_10000000_A_T",0.8],
            ["chr1_10150000_A_T", "chr1_10000000_A_T",0.7],
            ["chr1_10200000_A_T", "chr1_10000000_A_T",0.6],
            ["chr1_10250000_A_T", "chr1_10000000_A_T",0.5],
            ["chr1_10300000_A_T", "chr1_10000000_A_T",0.4],
            ["chr1_10350000_A_T", "chr1_10000000_A_T",0.3],
        ]
        validation_df = pd.DataFrame(validation, columns=["#variant","locus_id","r2_to_lead"])
        passes=[]
        for col in out.columns:
            if out[col].dtype in (float,int):
                try:
                    passes.append(np.allclose(out[col],validation_df[col]))
                except:
                    passes.append(False)
            else:
                    passes.append(out[col].equals(validation_df[col]))
        if not all(passes):
            msg = f"Output does not equal validation data!\n validation:\n{validation_df}\noutput:\n{out}"
            raise Exception(msg)


    def test_r2_chisq(self):
        """Test that dynamic r2
        """
        cols = {
            "chrom":"chrom1",
            "pos":"pos1",
            "ref":"ref1",
            "alt":"alt1",
            "pval":"pval1"
        }
        variants = [
            "1:10000000:1e-10",#peak 1
            "1:10100000:1e-5",
            "1:10200000:1e-5",
            "1:10300000:1e-5",
            "1:10400000:1e-5",
            "1:10500000:1e-5",
            "1:10150000:1e-5",
            "1:10250000:1e-5",
            "1:10350000:1e-5",
            "1:10450000:1e-5"
        ]
        data,c = create_loci(variants)
        ld_variants = [
            Variant(
                a[0],
                a[1],
                "A",
                "T"
            )
            for a in data
        ]
        r2_threshold = 5.0 #with pval 1e-10 and threshold 5.0, r2 threshold will be 0.119 aka 440kb or so
        df = pd.DataFrame(data,columns=c)
        df["locus_id"]=np.nan
        ld_api = PosLD(ld_variants)
        df_p1 = df[df["pval1"]<1e-6].copy()
        df_p2 = df.copy()
        # r2 threshold is 0.201, which means that everything farther away than 399 kb or so is rejected.
        # Therefore there should be 7 variants in the group incl. lead
        output = gws_fetch.ld_grouping(
            df_p1,
            df_p2,
            sig_threshold_2 = 1e-2,
            dynamic_r2=True,
            ld_threshold=r2_threshold,
            locus_width=1000000,
            overlap=False,
            ld_api=ld_api,
            columns=cols
        )
        out=output[["#variant","locus_id","r2_to_lead"]].sort_values(["#variant"]).reset_index(drop=True)
        validation = [
            ["chr1_10000000_A_T", "chr1_10000000_A_T", 1.0],
            ["chr1_10100000_A_T", "chr1_10000000_A_T", 0.8],
            ["chr1_10150000_A_T", "chr1_10000000_A_T", 0.7],
            ["chr1_10200000_A_T", "chr1_10000000_A_T", 0.6],
            ["chr1_10250000_A_T", "chr1_10000000_A_T", 0.5],
            ["chr1_10300000_A_T", "chr1_10000000_A_T", 0.4],
            ["chr1_10350000_A_T", "chr1_10000000_A_T", 0.3],
            ["chr1_10400000_A_T", "chr1_10000000_A_T", 0.2]
        ]
        validation_df = pd.DataFrame(validation, columns=["#variant","locus_id","r2_to_lead"])
        passes=[]
        for col in out.columns:
            if out[col].dtype in (float,int):
                try:
                    passes.append(np.allclose(out[col],validation_df[col]))
                except:
                    passes.append(False)
            else:
                    passes.append(out[col].equals(validation_df[col]))
        if not all(passes):
            msg = f"Output does not equal validation data!\n validation:\n{validation_df}\noutput:\n{out}"
            raise Exception(msg)

if __name__=="__main__":
    unittest.main()