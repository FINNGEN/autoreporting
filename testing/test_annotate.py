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

    def test_rename_dict(self):
        prefix="prefix"
        values=[1,2,3,4,5]
        values=[str(a) for a in values]
        rename=annotate.create_rename_dict(values,prefix)
        validate={"1":"prefix1","2":"prefix2","3":"prefix3","4":"prefix4","5":"prefix5"}
        self.assertEqual(rename,validate)

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
        args.functional_path=""
        args.fg_ann_version="r3"
        args.prefix=""
        correct_value_path="annotate_resources/ann_validate"
        args.batch_freq=False
        args.column_labels=["#chrom", "pos", "ref", "alt", "pval"]
        try:
            with open("annotate_resources/annotate_df.tsv","r") as f:
                #test case lines
                lines=f.readlines()
                test_cases=[[0,1,2,3,4],[0,3,4,5,6],[0,1,2,5,6],[0,1,2,3,4,5,6]]
                for i, tcase in enumerate(test_cases):
                    #load the lines that appear in test case
                    test_lines=[lines[x] for x in tcase]
                    tmp="".join(test_lines)
                    with io.StringIO(tmp) as args.annotate_fpath:
                        in_df = pd.read_csv(args.annotate_fpath,sep="\t")
                        out = annotate.annotate(df=in_df,gnomad_genome_path=args.gnomad_genome_path, gnomad_exome_path=args.gnomad_exome_path, batch_freq=args.batch_freq, finngen_path=args.finngen_path,fg_ann_version=args.fg_ann_version,
                            functional_path=args.functional_path, prefix=args.prefix, column_labels=args.column_labels).astype(object)
                        df2=pd.read_csv("{}_{}.csv".format(correct_value_path,i+1),sep="\t").astype(object)
                        pd.testing.assert_frame_equal(out,df2)
        except:
            self.assertTrue(False)



        

if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()