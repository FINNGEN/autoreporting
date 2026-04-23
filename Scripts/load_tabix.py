import gzip
import shlex,subprocess, glob, time
from typing import Dict, List,Generator,NamedTuple
from data_access.db import Variant
import pysam
import os
import sys
import contextlib

MAX_RETRIES=7

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
        if opts.fname.startswith("gs://"):
            print(f"file {opts.fname} is from google cloud")
            self.local = False
            self.refresh_token()
            self.token_refresh_time = time.time()
        else: #local file
            self.local=True
            if not os.path.exists(opts.fname):
                raise Exception(f"Resource {opts.fname} does not to exist")
        self.fileobject = pysam.TabixFile(opts.fname,encoding="utf-8")
        self.cpra = [
            opts.c,
            opts.p,
            opts.r,
            opts.a
        ]
        self.sequences = set(self.fileobject.contigs)
        #if header is stored in tabix-compatible format, use that
        if len(self.fileobject.header) != 0:
            self.header = self.fileobject.header[-1].split("\t")
        #else load from head of file
        else:
            with gzip.open(opts.fname) as f:
                header = f.readline().decode()
                self.header = header.strip("\n").split("\t")
                self.hdi = {a:i for i,a in enumerate(self.header)}

    def load_region(self, sequence:str, start:int, end:int, data_columns:List[str])-> Dict[Variant,List[str]]:
        """Load data columns for a region
        NOTE: needs better return format, this is just quite bad to work with
        """
        ## if not local, update refresh token
        if not self.local:
            # refresh tokens every 20 minutes
            current_time = time.time()
            if (current_time - self.token_refresh_time) > 1200.0:
                print("Refreshing GCS_AUTH_TOKEN",file=sys.stderr)
                self.refresh_token()
                self.token_refresh_time = current_time
        #ensure all columns are in header
        if any([a for a in data_columns if a not in self.header]):
            raise Exception(f"columns {[a for a in data_columns if a not in self.header]} not in header! Header: {self.header}")
        if sequence not in self.sequences:
            raise ValueError((f"Error in tabix file loading: The sequence {sequence} is not in tabix-indexed file {self.fname}. List of available sequences: {[a for a in self.sequences]}.\n"
                              f"Error occurred with region {sequence}:{start}-{end}"))
        
        col_idx = [self.hdi[a] for a in data_columns]
        out = {}

        tries = 0
        while True:
            try:
                iter = self.fileobject.fetch(sequence, max(start-1,0),end)
                for l in iter:
                    cols = l.split("\t")
                    vid = Variant(
                        cols[self.hdi[self.cpra[0]]],
                        int(cols[self.hdi[self.cpra[1]]]),
                        cols[self.hdi[self.cpra[2]]],
                        cols[self.hdi[self.cpra[3]]],
                    )
                    datacols = [cols[i] for i in col_idx]
                    out[vid] = datacols
                break
            except:
                print(f"Error loading data from region {sequence}:{start}-{end}",file=sys.stderr)
                #assume error is in accessing over gcp, so we retry
                if tries > MAX_RETRIES:
                    raise Exception(f"Accessing region {sequence}:{start}-{end} from file {self.fname} failed after {tries} tries!")
                
                print(f"Error when accessing data from tabix file {self.fname}, assuming error is with access over GCP. Waiting {2**tries} seconds, recycling fileobject.",file=sys.stderr)
                time.sleep(2**tries)
                self.restart_fileobject()
                tries +=1
        return out

    def close(self):
        self.fileobject.close()

    def refresh_token(self):
        tmp = subprocess.run(shlex.split("gcloud auth print-access-token"),capture_output=True,encoding="utf-8")
        os.environ["GCS_OAUTH_TOKEN"] = tmp.stdout.strip()

    def restart_fileobject(self):
        #close current fobj
        self.fileobject.close()
        self.refresh_token()
        #try to reopen the fileobject, with exponential backoff.
        tries = 0
        while True:
            try:
                self.fileobject = pysam.TabixFile(self.fname,encoding="utf-8")
                break
            except:
                if tries > MAX_RETRIES:
                    raise Exception(f"Could not reopen tabix fileobject {self.fname} with multiple tries.")
                print(f"Error refreshing fileobject for tabixfile {self.fname}, waiting {2**tries} seconds.",file=sys.stderr)
                time.sleep(2**tries)
                tries += 1

 
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
