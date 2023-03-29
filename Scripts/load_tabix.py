import gzip
from typing import Dict, List,Generator,NamedTuple
from data_access.db import Variant
import pysam
import os
import sys
import contextlib

class TabixOptions(NamedTuple):
    fname: str
    c:str
    p:str
    r:str
    a:str

class TabixResource:
    """Load variants from tabix-indexed resource
    Assumptions:
    - The resource has c,p,r,a in it.

    """
    def __init__(self,opts:TabixOptions):
        #check that this file and its tabix file exist
        self.fname = opts.fname
        if not os.path.exists(opts.fname):
            raise Exception(f"Resource {opts.fname} does not to exist")
        self.fileobject = pysam.TabixFile(opts.fname,encoding="utf-8")
        self.cpra = [
            opts.c,
            opts.p,
            opts.r,
            opts.a
        ]
        #if header is stored in tabix-compatible format, use that
        if len(self.fileobject.header) != 0:
            self.header = self.fileobject.header[-1].split("\t")
        #else load from head of file
        else:
            with gzip.open(opts.fname) as f:
                header = f.readline().decode()
                self.header = header.strip("\n").split("\t")

    def load_region(self, sequence:str, start:int, end:int, data_columns:List[str])-> Dict[Variant,List[str]]:
        """Load data columns for a region
        NOTE: needs better return format, this is just quite bad to work with
        """
        
        #ensure all columns are in header
        if any([a for a in data_columns if a not in self.header]):
            raise Exception(f"columns {[a for a in data_columns if a not in self.header]} not in header! Header: {self.header}")
        hdi = {a:i for i,a in enumerate(self.header)}
        col_idx = [hdi[a] for a in data_columns ]
        out = {}
        
        iter = self.fileobject.fetch(sequence, max(start-1,0),end)
        for l in iter:
            try:
                cols = l.split("\t")
                vid = Variant(
                    cols[hdi[self.cpra[0]]],
                    int(cols[hdi[self.cpra[1]]]),
                    cols[hdi[self.cpra[2]]],
                    cols[hdi[self.cpra[3]]],
                )
                datacols = [cols[i] for i in col_idx]
                out[vid] = datacols
            except:
                print("Exception handling for tabix not yet implemented",file=sys.stderr)
                raise
         
        return out

    def close(self):
        self.fileobject.close()

 
@contextlib.contextmanager
def tb_resource_manager(fname, c,p,r,a) -> Generator:
    """Context manager helper for TabixResource
    """
    opts = TabixOptions(fname,c,p,r,a)
    obj = TabixResource(opts)
    try:
        yield obj
    finally:
        obj.close()
