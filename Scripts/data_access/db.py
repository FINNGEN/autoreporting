#!/usr/bin/env python3

## File that contains abstract classes for different DAOs
import abc
from typing import List, Dict,Any, Optional, NamedTuple
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
    def get_range(self, variant: Variant, bp_range: int, ld_threshold: Optional[float]) -> List[LDData]:
        """Return LD for single variant range
        Args:
            variant (Variant): Variant for which to get the LD neighbourhood
            bp_range (int): LD calculation range. Get LD results for variants closer than this many basepairs away from the lead var.
            ld_threshold (Optional[float]): Optional LD R^2 threshold. Only variants that are in greater correlation than the threshold are reported.
        Returns:
            (List[LDData]):List of variant associations
        """
        pass

    def get_ranges(self, variants: List[Variant], bp_range: int,
                   thresholds: Optional[List[float]] = None,
                   bp_ranges: Optional[List[int]] = None,
                   workers: int = 1) -> Dict[Variant, List[LDData]]:
        """Return LD for a batch of variants, keyed by variant.

        Default serial implementation; subclasses (e.g. TabixLD) may parallelize.
        bp_ranges, when given, overrides bp_range per-variant; thresholds likewise.
        """
        out: Dict[Variant, List[LDData]] = {}
        for i, v in enumerate(variants):
            r = bp_ranges[i] if bp_ranges is not None else bp_range
            t = thresholds[i] if thresholds is not None else None
            out[v] = self.get_range(v, r, t)
        return out

    def make_fetcher(self, workers: int = 1) -> "LDFetcher":
        """Return a reusable fetcher for lazy, batched LD prefetch (see LDFetcher). Subclasses
        with expensive per-process setup (e.g. TabixLD's worker pool) override this to keep
        that setup alive across batches; the base fetcher is serial."""
        return LDFetcher(self)


class LDFetcher:
    """Fetch LD for a batch of leads, keyed by variant, reusing setup across calls.

    The grouping loop fetches LD lazily in small batches instead of prefetching every lead up
    front, which bounds peak memory to one batch and skips leads that get consumed as partners
    before they would be fetched. This base implementation is serial; close() is a no-op.
    """
    def __init__(self, ld: "LDAccess"):
        self.ld = ld

    def fetch(self, leads: List[Variant], bp_range: int) -> Dict[Variant, List[LDData]]:
        return {v: self.ld.get_range(v, bp_range, None) for v in leads}

    def close(self):
        pass

class Location(NamedTuple):
    """Chromosomal position
    """
    chromosome: str
    position: int

class Rsid(NamedTuple):
    """Rsid number
    """
    location: Location
    rsid: str

class RsidVar(NamedTuple):
    variant: Variant
    rsid: str

class VariantData(NamedTuple):
    """Potentially multiallelic variant
    """
    variant: Variant
    other_alts: List[str] #possibly empty
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
        pass

class CSVariant(NamedTuple):
    variant:Variant
    prob: float
    lead_r2:float

class CS(NamedTuple):
    variants:List[CSVariant]
    lead:Variant
    region: str
    number:int
    bayes:float
    min_r2:float
    size:int
    good_cs: bool

class CSAccess(object):
    """
    Abstract object for getting CS information
    """
    @abc.abstractmethod
    def get_cs(self)->List[CS]:
        """Returns all credible sets from a datasource
        Returns:
            (List[CS]): List of credible sets. Each credible set contains variants.
        """
        pass