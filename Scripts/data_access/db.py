#!/usr/bin/env python3

## File that contains abstract classes for different DAOs
import json, requests
import time
import abc
from typing import List, Text, Dict,Any
from io import StringIO
import pandas as pd, numpy as np

class ExtDB(object):

    @abc.abstractmethod
    def get_associations(self,chromosome: str,start: int,end: int,pval: float,size: int)-> List[Dict[str,Any]]:
        """ Return associations of range chr:start-end that have pval smaller than pval. Get results in at most size sized chunks.
            Args: chromosome start end pval size
            Returns: List of Dictionaries with elements "chrom":chromosome "pos":position "ref":ref_allele "alt":alt_allele "pval":p-value "trait":phenotype_code
        """
        return

    @abc.abstractmethod
    def get_trait(self, trait_code : str)-> str:
        """ Return trait given trait code
            Args: trait_code
            Returns: Trait name
        """
        return

class LDAccess(object):
    """
    Abstract object for getting LD
    """

    @abc.abstractmethod
    def get_ranges(self, variants: pd.DataFrame, window: int, ld_threshold: float) -> pd.DataFrame:#List[ Dict[str, Any ]]:
        """
        Return LD for multiple variant ranges
        In: variant data, i.e. a dataframe with columns [chr, pos, ref, alt, #variant], a window,ld threshold
        Out: Dataframe with LD information
        """
        return