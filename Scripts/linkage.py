import abc
from typing import List, Text, Dict,Any
import pandas as pd, numpy as np
from gwcatalog_api import try_request

class LDAccess(object):
    """
    Abstract object for getting LD
    """
    @abc.abstractmethod
    def get_range(self, chrom: str,pos: int,ref: str,alt: str, window: int)-> List[Dict[str,Any]]:
        """
        Return LD between a variant and a chromosomal range
        In: chromosome, variant ID, range in bases
        Out: list of dicts where each dict has the LD information
        """
        return


class OnlineLD(LDAccess):
    def __init__(self):
        self.url="http://api.finngen.fi/api/ld"
    
    def get_range(self, chrom,pos,ref,alt,window):
        variant="chr{}:{}:{}:{}".format(chrom, pos, ref, alt)
        params={"variant":variant,"panel":"sisu3","variant":variant,"window":window}
        data=try_request(self.url,params=params)
        if type(data) == type(None):
            raise Exception("Fatal error: LD api not responding.")
        return data


class PlinkLD(LDAccess):
    def get_range(self,chrom,pos,ref,alt,window):
        pass

c=OnlineLD()
print(c.get_range("X",130355577,"A","C",100000).json())