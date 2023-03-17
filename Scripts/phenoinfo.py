
from grouping_model import GroupingOptions, PhenoInfo
from enum import Enum, unique
from typing import NamedTuple, Optional
import sys 

class PhenoInfoOptions(NamedTuple):
    infofile: str
    phenotype:str


def get_phenotype_data(options:Optional[PhenoInfoOptions])->PhenoInfo:
    if options == None:
        return PhenoInfo("NA","NA","NA","NA","NA")
    with open(options.infofile) as f:
        hdr = f.readline().strip("\r\n").split("\t")
        hdi = {a:i for i,a in enumerate(hdr)}
        data = [a.strip("\r\n").split("\t") for a in f.readlines()]
        try:
            row = [a for a in data if a[hdi["phenocode"]] == options.phenotype][0]
            return PhenoInfo(row[hdi["phenocode"]],row[hdi["name"]],int(row[hdi["num_cases"]]),int(row[hdi["num_controls"]]),row[hdi["category"]])
        except:
            print(f"endpoint {options.phenotype} info was not loaded successfully from phenotype info file {options.infofile}",file=sys.stderr)
            return PhenoInfo("NA","NA","NA","NA","NA")