import abc
from typing import List, Text, Dict,Any
import pandas as pd, numpy as np
from data_access.db import ExtDB

class CustomCatalog(ExtDB):
    def __init__(self, fname: str, pval_threshold: float, padding: int):
        self.fname = fname
        self.pval_threshold = float(pval_threshold)
        self.pad = int(padding)
        self.data = pd.read_csv(fname,sep="\t",na_values="NA")
        #load data into dataframe
        #fix the following column types: pval, beta as float, pos as int, chrom as str
        self.data["pval"] = pd.to_numeric(self.data["pval"],errors="coerce")
        self.data["beta"] = pd.to_numeric(self.data["beta"],errors="coerce")
        self.data["pos"] = pd.to_numeric(self.data["pos"],errors="coerce") 
        self.data = self.data.astype({"chrom":"str"})
        self.data=self.data.dropna(axis="index",subset=["chrom","pos","ref","alt","pval"])

    def get_associations(self,chromosome: str,start: int,end: int)-> List[Dict[str,Any]]:
        start=max(0,int(start)-self.pad)
        end=int(end)+self.pad
        #filter by chromosome, pos, pval
        tmpdata=self.data.loc[self.data["chrom"] == chromosome,:]
        tmpdata=tmpdata.loc[(tmpdata["pos"]>= start) &(tmpdata["pos"]<= end),: ]
        tmpdata=tmpdata.loc[tmpdata["pval"]<=self.pval_threshold,:]
        tmpdata["trait_name"] = tmpdata["trait"].apply(lambda x: self.get_trait(x))
        rename_d = {"study_doi":"study_link"}
        return tmpdata.rename(columns=rename_d).to_dict("records")
        #return results

    def get_trait(self, trait_code: str) -> str:
        """Return trait name
        Return trait code as name as the custom data resource only has a trait defined, no trait code
        """
        return trait_code
    
    def associations_for_regions(self, regions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Return associations for a list of regions of type {"chrom": str, "min": int, "max": int }
        Args:
            regions (List[Dict[str, Any]]): The list of regions for which associations are queried
        """
        out= []
        for region in regions:
            assocs=self.get_associations(region["chrom"],region["min"],region["max"])
            out.extend(assocs)
        return out