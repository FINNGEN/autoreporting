#!/usr/bin/env python3

#create Factory that contains all databases, so they can be queried without caring about their implementations

#How about just an object that has e.g. a dict for those
from typing import List, Dict, Any
from data_access.db import ExtDB
from data_access import *

#TODO: change into ExtDB because it kinda already is
class CompoundDB(ExtDB):
    def __init__(self, conf: List[ExtDB]):
        self.db=[]
        for db in conf:
            self.db.append(db)
        self.get_trait=None
        self.get_associations=None

    def associations_for_regions(self, regions: List[Dict[str, Any]] ) -> List[Dict[str, Any]]:
        #for every database, query the regions and append the results
        results = []
        for db in self.db:
            results.extend( db.associations_for_regions( regions ) )
        return results

def db_factory(use_gwascatalog, custom_dataresource, database_choice, localdb_path, gwas_width, gwas_pval ,gwas_threads) -> ExtDB:
    gwapi=None
    customresource=None
    if use_gwascatalog:
        if database_choice=="local":
            gwapi=gwcatalog_api.LocalDB(localdb_path,gwas_pval,gwas_width)
            gwas_threads=1
        elif database_choice=="summary_stats":
            gwapi=gwcatalog_api.SummaryApi(gwas_pval, gwas_width,gwas_threads)
        else:
            gwapi=gwcatalog_api.GwasApi(gwas_pval, gwas_width,gwas_threads)
    if custom_dataresource != "":
        customresource = custom_catalog.CustomCatalog(custom_dataresource,gwas_pval,gwas_width)
    dataresources = [a for a in [gwapi, customresource] if a != None]
    datadb = CompoundDB(dataresources)
    return datadb