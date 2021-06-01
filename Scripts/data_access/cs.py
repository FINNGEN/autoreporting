import gzip
import os.path as path
from typing import List
from db import Variant, CSVariant, CS, CSAccess

class CSGzReader(CSAccess):
    """
    Read credible set data from bgzipped susie outputs
    """

    def __init__(self, snp_path: str, cred_path: str):
        """
        Initialize CSGzReager, i.e. make sure the files exist and that they contain the correct columns.
        """
        if (not path.isfile(snp_path)) or (not path.isfile(cred_path)):
            raise Exception(f"One of the credible set filepaths does not exist: '{snp_path}' or '{cred_path}'") 
        self.snp_path = snp_path
        self.cred_path = cred_path
        #make sure both have correct headers
        required_snp_headers = [
            "chromosome",
            "position",
            "allele1",
            "allele2",
            "prob",
            "region",
            "cs",
            "lead_r2"
        ]
        required_cred_headers = [
            "cs_log10bf",
            "cs_min_r2",
            "cs_size",
            "cs_number",
            "cs_region",
            "good_cs"
        ]
        with gzip.open(self.snp_path,"rt") as snp:
            headers = snp.readline().strip("\n").split("\t")
            if not all([a in headers for a in required_snp_headers]):
                missing_cols = [a for a in required_snp_headers if a not in headers]
                raise Exception(f"Required columns {','.join(missing_cols)} not in credible set snp file: {snp_path}")
        with gzip.open(self.cred_path,"rt") as cred:
            headers = cred.readline().strip("\n").split("\t")
            if not all([a in headers for a in required_cred_headers]):
                missing_cols = [a for a in required_cred_headers if a not in headers]
                raise Exception(f"Required columns {','.join(missing_cols)} not in credible set cred file: {cred_path}")
        
    
    def get_cs(self) -> List[CS]:
        raise NotImplementedError

class CSSummaryReader(CSAccess):
    """
    Read credible set data from susie result summaries
    """
    def __init__(self, snp_path: str, cred_path: str):
        """
        Initialize CSSummaryReader, i.e. make sure the files exist and that they contain the correct columns.
        Args:
            snp_path(str):
            cred_path(str)
        """
        if (not path.isfile(snp_path)) or (not path.isfile(cred_path)):
            raise Exception(f"One of the credible set filepaths does not exist: '{snp_path}' or '{cred_path}'") 
        required_cred_headers=[
            "region",
            "cs",
            "v",
            "cs_specific_prob",
            "cs_log10bf",
            "cs_min_r2",
            "cs_size",
            "good_cs"
        ]
        required_snp_headers=[
            "region",
            "cs",
            "v",
            "cs_specific_prob",
            "lead_r2"
        ]
        self.snp_path = snp_path
        self.cred_path = cred_path
        with open(self.snp_path,"r") as  snp:
            headers = snp.readline().strip("\n").split("\t")
            if not all([a in headers for a in required_snp_headers]):
                missing_cols = [a for a in required_snp_headers if a not in headers]
                raise Exception(f"Required columns {','.join(missing_cols)} not in credible set snp file: {snp_path}")
        with open(self.cred_path,"r") as cred:
            headers = cred.readline().strip("\n").split("\t")
            if not all([a in headers for a in required_cred_headers]):
                missing_cols = [a for a in required_cred_headers if a not in headers]
                raise Exception(f"Required columns {','.join(missing_cols)} not in credible set cred file: {cred_path}")

    def get_cs(self) -> List[CS]:

        raise NotImplementedError