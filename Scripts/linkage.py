import abc
from typing import List, Text, Dict,Any
import pandas as pd, numpy as np
from gwcatalog_api import try_request

class LDAccess(object):
    """
    Abstract object for getting LD
    """
    @abc.abstractmethod
    def get_range(self, chrom: str,pos: int,ref: str,alt: str, window: int)-> pd.DataFrame:
        """
        Return LD between a variant and a chromosomal range
        In: chromosome, variant ID, range in bases
        Out: Dataframe with LD information
        """
        return

    @abc.abstractmethod
    def get_ranges(self, variants: pd.DataFrame, window: int) -> List[ Dict[str, Any ]]:
        """
        Return LD for multiple variant ranges
        In: variant data, i.e. a dataframe with columns [chr, pos, ref, alt, #variant], a window
        Out: Dataframe with LD information
        """
        return


def parse_variant(variant):
    parsed_variant="chr{}".format( "_".join(variant.split(":")) )
    return parsed_variant


class OnlineLD(LDAccess):
    def __init__(self):
        self.url="http://api.finngen.fi/api/ld"
    
    def get_range(self, chrom,pos,ref,alt,window):
        variant="{}:{}:{}:{}".format(chrom, pos, ref, alt)
        params={"variant":variant,"panel":"sisu3","variant":variant,"window":window}
        data=try_request(self.url,params=params)
        if type(data) == type(None):
            print("LD data not found.") #TODO: handle this better. Perhaps return more information from try_request
            return None
        #parse data
        ld_data=data.json()["ld"]
        ld_data=pd.DataFrame(ld_data)
        ld_data=ld_data[ ["variation1", "variation2", "r2"] ]
        return ld_data

    def get_ranges(self, variants, window):
        #The api works with requesting chr:pos:ref:alt keys
        data=variants[ [ "chr", "pos", "ref", "alt" ] ]
        ld_data=pd.DataFrame()
        for idx, row in variants.iterrows():
            variant_ld = self.get_range(row["chr"],row["pos"],row["ref"],row["alt"],window)
            ld_data=pd.concat([ld_data,variant_ld],ignore_index=True,sort=False)
        #parse the LD data. parsing single columns at a time because result_type="expand" adds horrible (8x running time) overhead.
        ld_data["variant_1"] = ld_data["variation1"].apply(parse_variant)
        ld_data["variant_2"] = ld_data["variation2"].apply(parse_variant)
        ld_data["chrom_1"] = ld_data["variation1"].apply(lambda x: x.split(":")[0])
        ld_data["pos_1"] = ld_data["variation1"].apply(lambda x: x.split(":")[1])
        ld_data["chrom_2"] = ld_data["variation2"].apply(lambda x: x.split(":")[0])
        ld_data["pos_2"] = ld_data["variation2"].apply(lambda x: x.split(":")[1])
        ld_data=ld_data.astype(dtype={"pos_1":int,"pos_2":int})
        returncolumns=["chrom_1","pos_1","variant_1","chrom_2","pos_2","variant_2","r2"]
        ld_data=ld_data[returncolumns]
        return ld_data


class PlinkLD(LDAccess):
    def __init__(self,memory):
        self.memory=memory

    def get_range(self,chrom,pos,ref,alt,window):
        raise NotImplementedError("Not implemented yet. Use get_ranges.")

    def get_ranges(self, variants, window):
        #assume columns are: [chrom, pos, ref, alt,#variant], and window is int
        #calculate per-chromosome to save up on memory. Performance should be similar, if not faster.
        
        raise NotImplementedError("Not implemented yet")

# class LDStoreLD(LDAccess):
#     """
#     Might be a bit tricky, as the ld panel has to be separated into chromosomes.
#     """
#     def get_range(self,chrom,pos,ref,alt,window):
#         raise NotImplementedError("Not implemented yet. Use get_ranges.")

#     def get_ranges(self, variants, window):
#         raise NotImplementedError("Not implemented yet")


c=OnlineLD()
data=pd.DataFrame({ "chr":["X"], "pos":[130355577], "ref": ["A"], "alt":["C"]  })
print(c.get_ranges(data,1500000))