import gzip
import os.path as path
from typing import List, Dict
from collections import defaultdict
import pandas as pd
from data_access.db import Variant, CSVariant, CS, CSAccess

def _v_parser(s:str)->Variant:
    """
    parse variant format in cs summaries
    """
    data=s.split(":")
    return Variant(data[0],int(data[1]),data[2],data[3])

class CSFullReader(CSAccess):
    """
    Read credible set data from bgzipped susie outputs
    """

    def __init__(self, snp_path: str, cred_path: str):
        """
        Initialize CSFullReader, i.e. make sure the files exist and that they contain the correct columns.
        """
        if (not path.isfile(snp_path)) or (not path.isfile(cred_path)):
            raise Exception(f"One of the credible set filepaths does not exist: '{snp_path}' or '{cred_path}'") 
        self.snp_path = snp_path
        self.cred_path = cred_path
        #make sure both have correct headers
        required_snp_headers = [
            "v",
            "cs_specific_prob",
            "region",
            "cs",
            "lead_r2"
        ]
        required_cred_headers = [
            "cs_log10bf",
            "cs_min_r2",
            "cs",
            "cs_size",
            "region",
            "low_purity"
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
        """
        Read credible set data from full susie results
        """
        #snps, group them using their regions and cs numbers
        #also find out the lead variant for each of those groups
        #then, read in the credible sets, and assign those variants and lead vars to those.
        variants = defaultdict(list)
        with gzip.open(self.snp_path,"rt") as snp_file:
            headers = snp_file.readline().strip().split("\t")
            h_order = {a:i for (i,a) in enumerate(headers)}
            for dataline in snp_file:
                data=dataline.split("\t")
                #if the snp is not part of a credible set, skip it
                snp_cs = int(data[h_order["cs"]])
                if snp_cs != -1:
                    snp = _v_parser(data[h_order["v"]])
                    snp_prob = float(data[h_order["cs_specific_prob"]])
                    snp_r2 = float(data[h_order["lead_r2"]])
                    snp_region = data[h_order["region"]]
                    variants[(snp_region,snp_cs)].append(CSVariant(
                        snp,
                        snp_prob,
                        snp_r2
                    ))

        cs = []
        with gzip.open(self.cred_path,"rt") as cred_file:
            headers = cred_file.readline().strip().split("\t")
            h_order = {a:i for (i,a) in enumerate(headers)}
            for dataline in cred_file:
                data=dataline.strip().split("\t")
                reg = data[h_order["region"]]
                n = int(data[h_order["cs"]])
                bayes = float(data[h_order["cs_log10bf"]]) 
                min_r2 = float(data[h_order["cs_min_r2"]])
                size = int(data[h_order["cs_size"]])
                low_purity = (data[h_order["low_purity"]].lower() == "true")
                good_cs = not low_purity
                try:
                    cs_vars = variants[(reg,n)]
                    cs_lead = max(cs_vars, key=lambda x: x.prob)
                    cs.append(
                        CS(
                            cs_vars,
                            cs_lead.variant,
                            reg,
                            n,
                            bayes,
                            min_r2,
                            size, 
                            good_cs
                        )
                    )
                except:
                    raise Exception(f"CS with region {reg} and number {n} not found in credible set data!")
        return cs



class CSSummaryReader(CSAccess):
    """
    Read credible set data from susie result summaries
    members:
        snp_path: path to .snp.filter.tsv file
        cred_path: path to .cred.summary.tsv file
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
            "cs_specific_prob"#,
            #"lead_r2" #lead_r2 is not present in the summaries for some reason
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
        """
        Read credible set information from cred and snp files
        """
        #read credible sets
        credsets = []
        with open(self.cred_path,"r") as cred_file:
            headers = cred_file.readline().strip().split("\t")
            h_order = {a:i for (i,a) in enumerate(headers)}
            for dataline in cred_file:
                data=dataline.strip().split("\t")
                credsets.append(
                    CS(
                        [],
                        _v_parser(data[h_order["v"]]),
                        data[h_order["region"]],
                        int(data[h_order["cs"]]),
                        float(data[h_order["cs_log10bf"]]),
                        float(data[h_order["cs_min_r2"]]),
                        int(data[h_order["cs_size"]]),
                        bool(data[h_order["good_cs"]])
                    )
                )
        # create credible set group dictionary
        cred_set_assign = {(a.region,a.number):i for i,a in enumerate(credsets)} 
        # read snp data
        with open(self.snp_path, "r") as snp_file:
            headers = snp_file.readline().strip().split("\t")
            h_order = {a:i for (i,a) in enumerate(headers)}
            for dataline in snp_file:
                data = dataline.strip().split("\t")

                snp_v = _v_parser(data[h_order["v"]])
                snp_reg = data[h_order["region"]]
                snp_num = int(data[h_order["cs"]])
                snp_prob = float(data[h_order["cs_specific_prob"]])
                snp_r2 = float("nan")
                if snp_v == credsets[cred_set_assign[(snp_reg,snp_num)]].lead:
                    snp_r2 = 1.0
                snp = CSVariant(snp_v,snp_prob,snp_r2)
                #assign to correct credible set
                try:
                    credsets[cred_set_assign[(snp_reg,snp_num)]].variants.append(snp)
                except:
                    raise Exception(f"snp {snp} with region {snp_reg} and number {snp_num} could not be placed into a credible set: Credible set was not found!")
        return credsets


def cs_to_df(cs:List[CS],columns: Dict[str,str])->pd.DataFrame():
    """Kludge to convert from list of CS to a dataframe containing the same information.
    """
    output_columns = [
        columns["chrom"],
        columns["pos"],
        columns["ref"],
        columns["alt"],
        "cs_id",
        "cs_region",
        "cs_number",
        "cs_prob",
        "cs_log10bf",
        "cs_min_r2",
        "cs_size",
        "good_cs",
        "r2_to_lead"
    ]
    out = pd.DataFrame(columns=output_columns)
    for c in cs:
        var_data = [ [
            a.variant.chrom,
            a.variant.pos,
            a.variant.ref,
            a.variant.alt,
            f"chr{c.lead.chrom}_{c.lead.pos}_{c.lead.ref}_{c.lead.alt}_{c.number}",
            c.region,
            c.number,
            a.prob,
            c.bayes,
            c.min_r2,
            c.size,
            c.good_cs,
            a.lead_r2
        ] for a in c.variants]
        data = pd.DataFrame(var_data,columns=output_columns)
        out = pd.concat([out,data],join="inner",ignore_index=True)
    return out

