#!/usr/bin/env python3

## File that contains abstract classes for different DAOs
import abc
from typing import List, Text, Dict,Any, Optional, NamedTuple
from io import StringIO
import pandas as pd #type: ignore
from autoreporting_utils import Region

class ExtDB(object):
    """Abstract base class for association searches
    """

    @abc.abstractmethod
    def associations_for_regions(self, regions: List[Region]) -> List[Dict[str, Any]]:
        """Return associations for a list of regions
        Args:
            regions (List[Region]): The list of regions for which associations are queried
        """

class Variant(NamedTuple):
    """Variant class
        chrom: str
        pos: int
        ref: str
        alt: str
    """
    chrom: str
    pos: int
    ref: str
    alt: str


class LDData(NamedTuple):
    """LD information class for LDAccess api
        variant1: Variant
        variant2: Variant
        r2: float
    """
    variant1: Variant
    variant2: Variant
    r2: float
    def to_flat(self) -> Dict[str,Any]:
        """Helper method to flatten LDData
        Returns:
            Dict[str,Any]: Dictionary with keys (variant1,chrom1,pos1,ref1,alt1,variant2,chrom2,pos2,ref2,alt2,r2)
        """
        return {
            "chrom1":self.variant1.chrom,
            "pos1":self.variant1.pos,
            "ref1":self.variant1.ref,
            "alt1":self.variant1.alt,
            "chrom2":self.variant2.chrom,
            "pos2":self.variant2.pos,
            "ref2":self.variant2.ref,
            "alt2":self.variant2.alt,
            "r2":self.r2
        }

class LDAccess(object):
    """
    Abstract object for getting LD
    """

    @abc.abstractmethod
    def get_ranges(self, variants: List[Variant], window: int, ld_threshold: Optional[float]) -> List[LDData]:
        """Return LD for multiple variant ranges
        Args: variant data, i.e. a dataframe with columns [chr, pos, ref, alt, #variant], a window,ld threshold
            variants (List[LDInput]): List of input variants ()
            window (int):
            ld_threshold (Optional[float]): Optional LD R^2 threshold. Only variants that are in greater correlation than the threshold are reported. 
        Returns: 
            (List[LDData]):List of variant associations
        """
        return

class Location(NamedTuple):
    """Chromosomal position
    """
    chromosome: str
    position: int

class VariantData(NamedTuple):
    """Potentially multiallelic variant
    """
    variant: Variant
    other_alts: List[str] #posibly empty
    rsid: int

    def biallelic(self) -> bool:
        return len(self.other_alts) == 0

class AlleleDB(object):
    """
    Abstract object for getting alleles for c:p
    """
    @abc.abstractmethod
    def get_alleles(self, positions: List[Location]) -> List[VariantData]:
        """Get alleles for chromosomal positions
        Args:
            positions (List[Location]): List of genetic locations
        Returns:
            (List[VariantData]): As many of those locations 
        """
        return
