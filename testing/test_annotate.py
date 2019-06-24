import unittest
import sys,os,io
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
import pandas as pd,numpy as np
from Scripts import annotate

class Arg():
    def __init__(self):
        pass

class TestAnnotate(unittest.TestCase):
    def test_rename_dict(self):
        lst=["A","B","C"]
        result={"A":"Prefix_A","B":"Prefix_B","C":"Prefix_C"}
        d=annotate.create_rename_dict(lst,"Prefix_")
        self.assertEqual(d,result)
        lst=[]
        result={}
        d=annotate.create_rename_dict(lst,"Prefix_")
        self.assertEqual(d,result)

    def test_calculate_enrichment(self):
        #load "gnomad" df
        g_df=pd.read_csv("annotate_resources/enrichment_df",sep="\t")
        g_df=g_df.apply(pd.to_numeric,errors="ignore")
        fi_af_col="AF_fin"
        count_nfe_lst=["AC_1","AC_2","AC_3"]
        number_nfe_lst=["AN_1","AN_2","AN_3"]
        result=pd.Series([2.0,1e6,0.0])
        test=annotate.calculate_enrichment(g_df,fi_af_col,count_nfe_lst,number_nfe_lst)
        for idx,_val in result.iteritems():
            self.assertAlmostEqual(result[idx],test[idx]) 

    def test_annotate(self):
        #define args
        #needed files: gnomad_genome_tabixed file, exome gnomad, finngen annotations
        #also, gws file
        #test so that we either have 0 or more lines matching to the annotation files,
        #to test if they work with zero or more lines
        #best to make a df with lines 1-2 corresponding to genomes df, 3-4 to exomes, 
        #and 4-6 to finngen annotations
        #then, test it, choosing only some lines for the df file
        #this might be best to do with stringio, loading the relevant lines and passing that as the path argument 
        #for annotate
        args=Arg()
        args.gnomad_genome_path="annotate_resources/gnomad_genomes.tsv.gz"
        args.gnomad_exome_path="annotate_resources/gnomad_exomes.tsv.gz"
        args.finngen_path="annotate_resources/finngen_anno.tsv.gz"
        args.annotate_out="annotate_resources/test_out.csv"
        args.batch_freq=False
        args.column_labels=["#chrom", "pos", "ref", "alt", "pval"]
        try:
            with open("annotate_resources/annotate_df.tsv","r") as f:
                #test case lines
                lines=f.readlines()
                test_cases=[[0,1,2,3,4],[0,3,4,5,6],[0,1,2,5,6],[0,1,2,3,4,5,6]]
                for tcase in test_cases:
                    #load the lines that appear in test case
                    test_lines=[lines[x] for x in tcase]
                    tmp="".join(test_lines)
                    with io.StringIO(tmp) as args.annotate_fpath:
                        annotate.annotate(args)
                        with open(args.annotate_out,"r") as f:
                            df=pd.read_csv(f,sep="\t")
                            #TODO: check whether the file is correct
        except:
            self.assertTrue(False)



        

if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()