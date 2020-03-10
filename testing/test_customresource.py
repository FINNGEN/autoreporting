import unittest
import unittest.mock as mock
import sys,os,json, requests
import pandas as pd, numpy as np
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
from Scripts.data_access import custom_catalog
from io import StringIO

def create_data() -> pd.DataFrame:
    cols = ["chrom","pos","ref","alt","pval","beta","se","trait","study_doi","trait_name"]
    chroms = ["1"]*100
    pos = list(range(1,101))
    ref = ["A","C","G","T"]*25
    alt = ["C","A","T","G"]*25
    pval = [a/1000000 for a in pos]
    beta = [0.1,0.2,-0.1,-0.2]*25
    se=["NA"]*100
    trait = ["test_trait"]*100
    study_doi=["NA"]*100
    trait_name = ["test_trait"]*100
    data=pd.DataFrame({ cols[0]:chroms, cols[1]:pos, cols[2]:ref, cols[3]:alt, cols[4]:pval, cols[5]:beta, cols[6]:se, cols[7]:trait, cols[8]:study_doi, cols[9]:trait_name })
    return data

def create_catalog() -> custom_catalog.CustomCatalog:
    data=create_data()
    #create dummy tsv
    file_obj = StringIO()
    data.to_csv(file_obj,sep="\t",index=False)
    file_obj.seek(0)
    return custom_catalog.CustomCatalog(file_obj, 1.0, 0)

class TestCustomCat(unittest.TestCase):
    def test_init(self):
        try:
            catalog = create_catalog()
        except:
            self.assertTrue(False)
            raise

    def test_get(self):
        #test the get function 
        catalog = create_catalog()
        data=create_data()
        range_start=10
        range_end=20
        pval=1.0
        chromosome="1"
        validation_data = data.loc[data["chrom"]== chromosome,:].copy()
        validation_data=validation_data.loc[(validation_data["pos"] >=range_start) & (validation_data["pos"] <=range_end) ,:]
        validation_data=validation_data.loc[validation_data["pval"]<=pval,:]
        validation_data = validation_data.replace("NA",np.nan).reset_index(drop=True)
        validation_data=validation_data.rename(columns={"study_doi":"study_link"})
        out = pd.DataFrame(catalog.get_associations(chromosome,range_start,range_end)).reset_index(drop=True)
        for col in out.columns:
            pd.testing.assert_series_equal(out[col], validation_data[col])

    def test_get_trait(self):
        trait="TRAIT"
        catalog=create_catalog()
        self.assertEqual(catalog.get_trait(trait),trait)

if __name__=="__main__":
    unittest.main()