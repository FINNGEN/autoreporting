#!/usr/bin/env python3

## File that contains abstract classes for different DAOs
import abc
from typing import List, Text, Dict,Any, Optional, NamedTuple
from io import StringIO
import pandas as pd #type: ignore

class ExtDB(object):
    """Abstract base class for association searches
    """

    @abc.abstractmethod
    def associations_for_regions(self, regions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Return associations for a list of regions of type {"chrom": str, "min": int, "max": int }
        Args:
            regions (List[Dict[str, Any]]): The list of regions for which associations are queried
        """

class LDInput(NamedTuple):
    """Variant class for LDAccess api
        variant: str
        chrom: str
        pos: int
        ref: str
        alt: str
    """
    variant: str
    chrom: str
    pos: int
    ref: str
    alt: str

    def __eq__(self, other): 
        if not isinstance(other, LDInput):
            return NotImplemented
        return self.variant == other.variant

class LDData(NamedTuple):
    """LD information class for LDAccess api
        variant1: str
        variant2: str
        chrom1: str
        chrom2: str
        pos1: int
        pos2: int
        r2: float
    """
    variant1: str
    variant2: str
    chrom1: str
    chrom2: str
    pos1: int
    pos2: int
    r2: float

class LDAccess(object):
    """
    Abstract object for getting LD
    """

    @abc.abstractmethod
    def get_ranges(self, variants: List[LDInput], window: int, ld_threshold: Optional[float]) -> List[LDData]:
        """Return LD for multiple variant ranges
        Args: variant data, i.e. a dataframe with columns [chr, pos, ref, alt, #variant], a window,ld threshold
            variants (List[LDInput]): List of input variants ()
            window (int):
            ld_threshold (Optional[float]): Optional LD R^2 threshold. Only variants that are in greater correlation than the threshold are reported. 
        Returns: 
            (List[LDData]):List of variant associations
        """
        return