#!/usr/bin/env python3
from data_access.cs import CSFullReader, CSSummaryReader
from data_access.db import CSAccess


def csfactory(cred_set_file:str)->CSAccess:
    """
    Return an access object to the credible set data
    """
    snp_file = cred_set_file
    if ".bgz" in cred_set_file:
        cred_file = snp_file.replace("snp","cred")
        return CSFullReader(snp_file,cred_file)
    elif ".snp.filter.tsv" in cred_set_file:
        cred_file = snp_file.replace(".snp.filter",".cred.summary")
        return CSSummaryReader(snp_file, cred_file)
    else:
        raise ValueError(f"credsetfile format not recognized for file {cred_set_file}.")