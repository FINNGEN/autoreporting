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

    

if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()