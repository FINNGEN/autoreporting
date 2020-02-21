import abc
from typing import List, Text, Dict,Any
import pandas as pd, numpy as np
from gwcatalog_api import ExtDB


class CustomCatalog(ExtDB):
    def __init__(self, fname):
        self.fname=fname
        self.data = pd.read_csv(fname,sep="\t",na_values="NA")
        #load data into dataframe
        #fix the following column types: pval, beta as float, pos as int, chrom as str
        self.data["pval"] = pd.to_numeric(self.data["pval"],errors="coerce")
        self.data["beta"] = pd.to_numeric(self.data["beta"],errors="coerce")
        self.data["pos"] = pd.to_numeric(self.data["pos"],errors="coerce") 
        self.data = self.data.astype({"#chrom":"str"})
        self.data=self.data.dropna(axis="index",subset=["#chrom","pos","ref","alt","pval"])

    def get_associations(self,chromosome: str,start: int,end: int,pval: float,size: int=0)-> List[Dict[str,Any]]:
        #filter by chromosome, pos, pval
        tmpdata=self.data.loc[self.data["#chrom"] == chromosome,:]
        tmpdata=tmpdata.loc[(tmpdata["pos"]>= start) &(tmpdata["pos"]<= end),: ]
        tmpdata=tmpdata.loc[tmpdata["pval"]<=pval,:]
        return tmpdata.to_dict("records")
        #return results

    def get_trait(self, trait_code: str) -> str:
        """Return trait name
        Return trait code as name as the custom data resource only has a trait defined, no trait code
        """
        return trait_code