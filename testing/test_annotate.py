import unittest
import unittest.mock as mock
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
        g_df=pd.read_csv("testing/annotate_resources/enrichment_df",sep="\t")
        g_df=g_df.apply(pd.to_numeric,errors="ignore")
        fi_af_col="AF_fin"
        count_nfe_lst=["AC_1","AC_2","AC_3"]
        number_nfe_lst=["AN_1","AN_2","AN_3"]
        result=pd.Series([2.0,1e6,0.0])
        test=annotate.calculate_enrichment(g_df,fi_af_col,count_nfe_lst,number_nfe_lst)
        for idx,_val in result.iteritems():
            self.assertAlmostEqual(result[idx],test[idx]) 
    """
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
        args.gnomad_genome_path="testing/annotate_resources/gnomad_genomes.tsv.gz"
        args.gnomad_exome_path="testing/annotate_resources/gnomad_exomes.tsv.gz"
        args.finngen_path="testing/annotate_resources/finngen_anno.tsv.gz"
        args.annotate_out="testing/annotate_resources/test_out.csv"
        args.functional_path=None
        args.previous_release_path=None
        args.prefix=""
        correct_value_path="testing/annotate_resources/ann_validate"
        args.batch_freq=False
        columns={"chrom":"#chrom", "pos":"pos", "ref":"ref", "alt":"alt", "pval":"pval", "beta":"beta", "af":"af"}
        try:
            with open("testing/annotate_resources/annotate_df.tsv","r") as f:
                #test case lines
                lines=f.readlines()
                test_cases=[[0,1,2,3,4],[0,3,4,5,6],[0,1,2,5,6],[0,1,2,3,4,5,6]]
                for i, tcase in enumerate(test_cases):
                    #load the lines that appear in test case
                    test_lines=[lines[x] for x in tcase]
                    tmp="".join(test_lines)
                    with io.StringIO(tmp) as args.annotate_fpath:
                        in_df = pd.read_csv(args.annotate_fpath,sep="\t")
                        out = annotate.annotate(df=in_df,gnomad_genome_path=args.gnomad_genome_path, gnomad_exome_path=args.gnomad_exome_path, batch_freq=args.batch_freq, finngen_path=args.finngen_path,
                            functional_path=args.functional_path, previous_release_path=args.previous_release_path, prefix=args.prefix, columns=columns).astype(object)
                        df2=pd.read_csv("{}_{}.csv".format(correct_value_path,i+1),sep="\t").astype(object)
                        for col in df2.columns:
                            pd.testing.assert_series_equal(out[col],df2[col])
        except:
            self.assertTrue(False)
    """
    def test_func_anno(self):
        """Test functional annotation
        Cases:
            Empty data
            Missing gzip file
            Missing tabix file
            Data available
        """
        columns={"chrom":"#chrom", "pos":"pos", "ref":"ref", "alt":"alt", "pval":"pval", "beta":"beta", "af":"af"}
        #empty data
        input_data = pd.DataFrame()
        retval = annotate.functional_annotate(input_data,None,columns)
        self.assertTrue(retval.empty)
        self.assertEqual(retval.columns,["#variant"])
        #missing gzip file
        input_data = "testing/annotate_resources/annotate_df.tsv"
        annotation_file = "nofile" 
        ann_df = pd.read_csv(input_data,sep="\t")
        with self.assertRaises(FileNotFoundError):
            retval = annotate.functional_annotate(ann_df,annotation_file,columns)
        #missing tabix file
        annotation_file = "testing/annotate_resources/annotate_df.tsv"
        with self.assertRaises(FileNotFoundError):
            retval = annotate.functional_annotate(ann_df,annotation_file,columns)
        #with proper data
        #prepare data
        rename_dict = {
            "#chrom":"chrom"
        }
        ann_df = pd.read_csv(input_data,sep="\t")
        ann_df["consequence"] = "missense_variant"
        func_mock_df=ann_df.rename(columns=rename_dict)[[
            "chrom",
            "pos",
            "ref",
            "alt",
            "consequence"
        ]].astype(dtype={"chrom":str})
        validation = ann_df[["#variant","consequence"]].rename(columns={"consequence":"functional_category"})
        #mock the file loading so it's not actually loaded, so we can use any other file
        with mock.patch("Scripts.annotate.load_tb_df",return_value=func_mock_df):
            out = annotate.functional_annotate(ann_df,"testing/annotate_resources/finngen_anno.tsv.gz",columns)
        self.assertTrue(out.equals(validation))

    def test_fg_anno(self):
        """Test annotate.finngen_annotate
        Cases:
            Empty data
            Missing gzip file
            Missing tabix file
            Data available
        """
        columns={"chrom":"#chrom", "pos":"pos", "ref":"ref", "alt":"alt", "pval":"pval", "beta":"beta", "af":"af"}
        #empty data
        input_data = pd.DataFrame()
        retval = annotate.finngen_annotate(input_data,None,False,columns)
        self.assertTrue(retval.empty)
        self.assertEqual(retval.columns,["#variant"])
        #missing gzip file
        input_data = "testing/annotate_resources/annotate_df.tsv"
        annotation_file = "nofile" 
        ann_df = pd.read_csv(input_data,sep="\t")
        with self.assertRaises(FileNotFoundError):
            retval = annotate.finngen_annotate(ann_df,annotation_file,False,columns)
        #missing tabix file
        annotation_file = "testing/annotate_resources/annotate_df.tsv"
        with self.assertRaises(FileNotFoundError):
            retval = annotate.finngen_annotate(ann_df,annotation_file,False,columns)
        #proper data
        fg_file = "testing/annotate_resources/finngen_anno.tsv.gz"
        out = annotate.finngen_annotate(ann_df,fg_file,False,columns)
        validation = {"#variant":["chr23_23_G_C","chr23_230_C_G"],
            "most_severe_gene":["GENE3","GENE3"],
            "most_severe_consequence":["non_coding_transcript_exon_variant","missense_variant"]}
        validation = pd.DataFrame(validation)
        self.assertTrue(validation.equals( out))


    def test_gnomad_gen(self):
        """Test annotate.gnomad_gen_annotate
        Cases:
            Empty data
            Missing gzip file
            Missing tabix file
            Data available
        """
        #Empty data
        columns={"chrom":"#chrom", "pos":"pos", "ref":"ref", "alt":"alt", "pval":"pval", "beta":"beta", "af":"af"}
        #empty data
        input_data = pd.DataFrame()
        retval = annotate.gnomad_gen_annotate(input_data,None,columns)
        self.assertTrue(retval.empty)
        self.assertEqual(retval.columns,["#variant"])
        #missing gzip file
        input_data = "testing/annotate_resources/annotate_df.tsv"
        annotation_file = "nofile" 
        ann_df = pd.read_csv(input_data,sep="\t")
        with self.assertRaises(FileNotFoundError):
            retval = annotate.gnomad_gen_annotate(ann_df,annotation_file,columns)
        #missing tabix file
        annotation_file = "testing/annotate_resources/annotate_df.tsv"
        with self.assertRaises(FileNotFoundError):
            retval = annotate.gnomad_gen_annotate(ann_df,annotation_file,columns)
        #proper data
        gnomad_file = "testing/annotate_resources/gnomad_genomes.tsv.gz"
        out = annotate.gnomad_gen_annotate(ann_df,gnomad_file,columns)
        validation = pd.read_csv("testing/annotate_resources/gnomad_gen_validate.tsv",sep="\t")
        self.assertTrue(validation.equals( out))
        

    def test_gnomad_exo(self):
        """Test annotate.gnomad_exo_annotate
        Cases:
            Empty data
            Missing gzip file
            Missing tabix file
            Data available
        """
        #Empty data
        columns={"chrom":"#chrom", "pos":"pos", "ref":"ref", "alt":"alt", "pval":"pval", "beta":"beta", "af":"af"}
        #empty data
        input_data = pd.DataFrame()
        retval = annotate.gnomad_exo_annotate(input_data,None,columns)
        self.assertTrue(retval.empty)
        self.assertEqual(retval.columns,["#variant"])
        #missing gzip file
        input_data = "testing/annotate_resources/annotate_df.tsv"
        annotation_file = "nofile" 
        ann_df = pd.read_csv(input_data,sep="\t")
        with self.assertRaises(FileNotFoundError):
            retval = annotate.gnomad_exo_annotate(ann_df,annotation_file,columns)
        #missing tabix file
        annotation_file = "testing/annotate_resources/annotate_df.tsv"
        with self.assertRaises(FileNotFoundError):
            retval = annotate.gnomad_exo_annotate(ann_df,annotation_file,columns)
        #proper data
        gnomad_file = "testing/annotate_resources/gnomad_exomes.tsv.gz"
        out = annotate.gnomad_exo_annotate(ann_df,gnomad_file,columns)
        validation = pd.read_csv("testing/annotate_resources/gnomad_exo_validate.tsv",sep="\t")
        self.assertTrue(validation.equals( out))

if __name__=="__main__":
    unittest.main()