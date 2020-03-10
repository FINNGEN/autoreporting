#!/usr/bin/env python3

#create Factory that contains all databases, so they can be queried without caring about their implementations

#How about just an object that has e.g. a dict for those
from typing import List, Dict, Any
from db import ExtDB

class DataFactory(object):
    def __init__(self, conf: List[ExtDB]):
        self.db=[]
        for db in conf:
            self.db.append(db)


    def query_variants(self, regions: List[Dict[str, Any]] ) -> List[Dict[str, Any]]  :
        #for every database, query the regions and append the results
        results = []
        for db in self.db:
            results.extend( db.associations_by_regions( regions ) )
        return results

def construct_datafactory(use_gwascatalog, customresource_path, database_choice, localdb_path) -> DataFactory:
    gwapi=None
    customresource=None
    if use_gwascatalog:
        if database_choice=="local":
            gwapi=gwcatalog_api.LocalDB(args.localdb_path)
            args.gwascatalog_threads=1
        elif args.database_choice=="summary_stats":
            gwapi=gwcatalog_api.SummaryApi()
        else:
            gwapi=gwcatalog_api.GwasApi()
    if args.custom_dataresource != "":
        customdataresource = custom_catalog.CustomCatalog(args.custom_dataresource)
    dataresources = [a for a in [gwapi, customdataresource] if a != None]
    factory = DataFactory(dataresources)
    return factory