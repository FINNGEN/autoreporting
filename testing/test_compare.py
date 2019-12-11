import unittest
import unittest.mock as mock
import sys,os,json
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
import pandas as pd,numpy as np
from Scripts import compare
from Scripts import autoreporting_utils as autils
from Scripts import gwcatalog_api

def top_merge(df,summary_df,columns):
    summary_df=compare.map_column(summary_df,"map_variant",columns)
    df=compare.map_column(df,"map_variant",columns)
    necessary_columns=["pval","#variant","map_variant","trait","trait_name"]
    report_out_df=pd.merge(df,summary_df.loc[:,necessary_columns],how="left",on="map_variant")
    report_out_df=report_out_df.drop(columns=["map_variant"])
    report_out_df=report_out_df.rename(columns={"#variant_x":"#variant","#variant_y":"#variant_hit","pval_x":columns["pval"],"pval_y":"pval_trait"})
    report_out_df=report_out_df.sort_values(by=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],"#variant"])
    return report_out_df

class TestGws(unittest.TestCase):

    def test_solve_indels(self):
        #test: test only chrom, pos, ref alt. Other columns need not be tested.
        #Necessary data: df which contains our versions of indels, indels
        columns={"chrom":"chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        indels_={"chrom":["1","1","X"], "pos":[2,3,4],"ref":["-","TATA","-"],"alt":["CGTACGTA","-","-"],"pval":[0.4,0.2,0.1],"code":[10,10,10],"trait":["1","2","3"]}
        df_={"chrom":["1","1"],"pos":[1,2],"ref":["A","TTATA"],"alt":["ACGTACGTA","T"]}
        indels=pd.DataFrame(indels_)
        df=pd.DataFrame(df_)
        retval=compare.solve_indels(indels,df,columns)
        #check that the retval is correct.
        #Match alleles to gwascatalog ref=ref and gwas alt=alt
        out_={"chrom":["1","1"],"pos":[1,2],"ref":["A","TTATA"],"ref":["A","TTATA"],"alt":["ACGTACGTA","T"],"pval":[0.4,0.2],"code":[10,10],"trait":["1","2"]}
        out=pd.DataFrame(out_)
        test_cols=["chrom","pos","ref","alt"]
        for col in test_cols:
            self.assertTrue(retval[col].astype(str).equals(out[col].astype(str)) )

    def test_map_column(self):
        chrom=["1","1","2","2"]
        pos=[100,200,300,400]
        ref=["A","T","C","T"]
        alt=["C","G","G","A"]
        df=pd.DataFrame({"#chrom":chrom,"pos":pos,"ref":ref,"alt":alt})
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt"}
        res_ref=["A","A","C","A"]
        res_alt=["C","C","G","T"]
        df2=pd.DataFrame({"#chrom":chrom,"pos":pos,"ref":res_ref,"alt":res_alt})
        df2["map_variant"]=autils.create_variant_column(df2)
        res=compare.map_column(df,"map_variant",columns)
        self.assertEqual(list(df2["map_variant"]),list(res["map_variant"]))

    def test_top_report(self):
        #Needed: dataframe, summary variant dataframe, end result dataframe
        #also need to test using empty summary dataframe as well as empty dataframe
        # test one: empty dataframe, empty summary variant dataframe, should yield an empty dataframe
        cols=["#chrom","pos","ref","alt","pval","#variant","locus_id"]
        summary_cols=["#chrom","pos","ref","alt","pval","#variant","trait","trait_name"]
        end_result_cols=["locus_id", "chr", "start", "end", "enrichment", "most_severe_gene", "most_severe_consequence", "lead_pval", "found_associations_strict", "found_associations_relaxed", "credible_set_variants", "functional_variants_strict", "functional_variants_relaxed"]
        traits=[]
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        df=pd.DataFrame(columns=cols)
        summary_df=pd.DataFrame(columns=summary_cols)
        raport_df = top_merge(df,summary_df,columns)
        res=compare.create_top_level_report(raport_df,traits,columns,"ld",0.5)
        validate=pd.DataFrame(columns=end_result_cols)
        self.assertTrue(res.equals(validate) )
        # test two: populated dataframe, empty summary variant dataframe
        df=pd.read_csv("compare_resources/top_df.csv",sep="\t")
        df["cs_id"]=np.nan
        df["cs_prob"]=np.nan
        raport_df = top_merge(df,summary_df,columns)
        res=compare.create_top_level_report(raport_df,traits,columns,"ld",0.5)
        validate=pd.read_csv("compare_resources/top_result_empty_summary.csv",sep="\t")
        validate=validate.fillna("")
        for col in validate.columns:
            self.assertTrue(res[col].astype(object).equals(validate[col].astype(object)) )
        # test 3: populated dataframe, populated summary variant df, no traits
        summary_df=pd.read_csv("compare_resources/summary_df.csv",sep="\t")
        raport_df = top_merge(df,summary_df,columns)
        res=compare.create_top_level_report(raport_df,traits,columns,"ld",0.5)
        validate=pd.read_csv("compare_resources/top_result_traits_1.csv",sep="\t")
        validate=validate.fillna("")
        for col in validate.columns:
            self.assertTrue(res[col].astype(object).equals(validate[col].astype(object)) )
        # test 4: populated dataframes, common traits
        traits=["TRAIT_1","TRAIT_4"]
        res=compare.create_top_level_report(raport_df,traits,columns,"ld",0.5)
        validate=pd.read_csv("compare_resources/top_result_traits_2.csv",sep="\t")
        validate=validate.fillna("")
        for col in validate.columns:
            self.assertTrue(res[col].astype(object).equals(validate[col].astype(object)) )
    
    def test_load_summaries(self):
        """Test cases:
        Normal
        Empty files
        either file does not exist, should throw a normal FileNotFoundError
        File specified in summary_fpath does not exist (should throw an understandable error message, preferably with information about which file was not found, and in which file was that file declared)
        """
        #normal case: summary file with 2 entries, entries have some variants 
        summ_fpath="compare_resources/summary_fpath_1"
        endpoint_fpath="compare_resources/endpoint_fpath_1"
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        df=compare.load_summary_files(summ_fpath,endpoint_fpath,columns)
        validate=pd.read_csv("compare_resources/loaded_summary_files_1.csv",sep="\t")
        for col in df.columns:
            self.assertEqual(list(df[col]),list(validate[col]))
        #empty files
        summ_fpath="compare_resources/summary_fpath_2"
        endpoint_fpath="compare_resources/endpoint_fpath_2"
        df=compare.load_summary_files(summ_fpath,endpoint_fpath,columns)
        validate=pd.DataFrame(columns=["#chrom", "pos", "ref", "alt", "pval", "#variant", "trait", "trait_name"])
        self.assertTrue(df.equals(validate))
        #either file does not exist
        summ_fpath="compare_resources/summary_fpath_2"
        endpoint_fpath="DOES_NOT_EXIST"
        with self.assertRaises(FileNotFoundError) as notfound:
            df=compare.load_summary_files(summ_fpath,endpoint_fpath,columns)
        #file does not exist
        summ_fpath="compare_resources/summary_fpath_no_file"
        endpoint_fpath="compare_resources/endpoint_fpath_1"
        with self.assertRaises(FileNotFoundError) as notfound:
            df=compare.load_summary_files(summ_fpath,endpoint_fpath,columns)
        #summary and endpoint amounts do not agree
        summ_fpath="compare_resources/summary_fpath_2"
        endpoint_fpath="compare_resources/endpoint_fpath_1"
        with self.assertRaises(RuntimeError) as notfound:
            df=compare.load_summary_files(summ_fpath,endpoint_fpath,columns)

    def test_api_summaries_local(self):
        def post_response(url,headers,data):
            retval=mock.Mock()
            with open("compare_resources/ensembl_response","r") as f:
                val=f.readline().strip("\n")
                retval.json=lambda:json.loads(val)
            retval.status_code=200
            return retval
        # test loading summary statistics using local gwascatalog
        # test with no file found, should error out
        # test with some variants, of which some are found and some are not
        # include a small copy of the database, at most 4 variants.
        localdb_path="compare_resources/5_line_gwascatalog.csv"
        chrom=["1","6"]
        pos=[100000,32441740]
        ref=["A","T"]
        alt=["T","G"]
        pval=[2e-8,1e-9]
        df=pd.DataFrame({"#chrom":chrom,"pos":pos,"ref":ref,"alt":alt,"pval":pval})
        gwascatalog_pad=500
        gwascatalog_pval=5e-8
        database_choice="local"
        gwascatalog_threads=1
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        #case one: wrong path to db
        wrong_path_to_db="wrong/filepath/to_db.tsv"
        with self.assertRaises(FileNotFoundError) as err:
            gwapi=gwcatalog_api.LocalDB(wrong_path_to_db)
            tmp=compare.load_api_summaries(df,gwascatalog_pad,gwascatalog_pval,gwapi,gwascatalog_threads,columns)
        #case two: correct usage of function, should return correct amount of results
        with mock.patch('Scripts.gwcatalog_api.requests.post',side_effect=post_response) as mock_post:
            gwapi=gwcatalog_api.LocalDB(localdb_path)
            gwas_df=compare.load_api_summaries(df,gwascatalog_pad,gwascatalog_pval,gwapi,gwascatalog_threads,columns)
        validate_path="compare_resources/local_gwascatalog_validate.csv"
        validate=pd.read_csv(validate_path,sep="\t",dtype={"chrom":str,"pos":int,"ref":str,"alt":str,"pval":float})
        for col in gwas_df.columns:
            for _,row in gwas_df.iterrows(): 
                if type(gwas_df.loc[_,col])==str:
                    self.assertEqual(gwas_df.loc[_,col],validate.loc[_,col])
                else:
                    self.assertAlmostEqual(gwas_df.loc[_,col],validate.loc[_,col])

if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()