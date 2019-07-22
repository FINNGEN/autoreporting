import unittest
import sys,os
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
import pandas as pd,numpy as np
from Scripts import compare
from Scripts import autoreporting_utils as autils

class TestGws(unittest.TestCase):

    def test_solve_indels(self):
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
        for col in out.columns:
            self.assertTrue(retval[col].equals(out[col]) )

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
        end_result_cols=["locus_id","chr","start","end","matching_pheno_gwas_catalog_hits","other_gwas_hits"]
        traits=[]
        columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"}
        df=pd.DataFrame(columns=cols)
        summary_df=pd.DataFrame(columns=summary_cols)
        res=compare.create_top_level_report(df,summary_df,traits,columns)
        validate=pd.DataFrame(columns=end_result_cols)
        self.assertTrue(res.equals(validate) )
        # test two: populated dataframe, empty summary variant dataframe
        df=pd.read_csv("compare_resources/top_df.csv",sep="\t")
        res=compare.create_top_level_report(df,summary_df,traits,columns)
        validate=pd.read_csv("compare_resources/top_result_empty_summary.csv",sep="\t")
        validate=validate.fillna("")
        for col in res.columns:
            self.assertTrue(res[col].equals(validate[col].astype(object)) )
        # test 3: populated dataframe, populated summary variant df, no traits
        summary_df=pd.read_csv("compare_resources/summary_df.csv",sep="\t")
        res=compare.create_top_level_report(df,summary_df,traits,columns)
        validate=pd.read_csv("compare_resources/top_result_traits_1.csv",sep="\t")
        validate=validate.fillna("")
        for col in res.columns:
            self.assertTrue(res[col].equals(validate[col].astype(object)) )
        # test 4: populated dataframes, common traits
        traits=["TRAIT_1","TRAIT_4"]
        res=compare.create_top_level_report(df,summary_df,traits,columns)
        validate=pd.read_csv("compare_resources/top_result_traits_2.csv",sep="\t")
        validate=validate.fillna("")
        for col in res.columns:
            self.assertTrue(res[col].equals(validate[col].astype(object)) )
    
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
            tmp=compare.load_api_summaries(df,gwascatalog_pad,gwascatalog_pval,database_choice,gwascatalog_threads,wrong_path_to_db,columns)
        #case two: correct usage of function, should return correct amount of results
        gwas_df=compare.load_api_summaries(df,gwascatalog_pad,gwascatalog_pval,database_choice,gwascatalog_threads,localdb_path,columns)
        validate_path="compare_resources/local_gwascatalog_validate.csv"
        validate=pd.read_csv(validate_path,sep="\t")
        for col in gwas_df.columns:
            self.assertEqual(list(gwas_df[col].astype(str)),list(validate[col].astype(str)))

if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()