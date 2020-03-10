#!/usr/bin/env python3

## File that contains abstract classes for different DAOs
import abc
from typing import List, Text, Dict,Any
from io import StringIO
import pandas as pd

class ExtDB(object):
    """Abstract base class for association searches
    """
    @abc.abstractmethod
    def get_associations(self, chromosome: str, start: int, end: int)-> List[Dict[str, Any]]:
        """ Return associations of range chr:start-end that have pval smaller than pval. Get results in at most size sized chunks.
            Args: chromosome start end pval size
            Returns: List of Dictionaries with elements "chrom":chromosome "pos":position "ref":ref_allele "alt":alt_allele "pval":p-value "trait":phenotype_code
        """
        return

    @abc.abstractmethod
    def get_trait(self, trait_code: str)-> str:
        """ Return trait given trait code
            Args: trait_code
            Returns: Trait name
        """
        return

    @abc.abstractmethod
    def associations_for_regions(self, regions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Return associations for a list of regions of type {"chrom": str, "min": int, "max": int }
        Args:
            regions (List[Dict[str, Any]]): The list of regions for which associations are queried
        """

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