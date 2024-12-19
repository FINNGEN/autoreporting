import argparse,shlex,subprocess, os
from subprocess import Popen, PIPE
from typing import List, Dict, NamedTuple
"""
Utility functions that are used in the scripts, put here for keeping the code clearer
"""

def filebasename(s):
    if s != "":
        return os.path.basename(s).split(".")[0]
    return ""


class Region(NamedTuple):
    chrom: str
    start: int
    end: int

    def overlaps(self, other: 'Region')->bool:
        """Check if two Regions overlap
        """
        #check that both regions are valid
        if self.end < self.start:
            raise Exception(f"Region {self} is invalid: start is larger than end!")
        if other.end < other.start:
            raise Exception(f"Region {other} is invalid: start is larger than end!")
        if self.chrom == other.chrom:
            if (self.start <= other.end) and (other.start <= self.end):
                return True
        return False

def prune_regions(regions:List[Region])->List[Region]:
    """merge overlapping regions, so that there are less overlapping regions
    Args:
        regions (List[Region]): List of regions to merge
    Returns:
        (List[Region]): List of non-overlapping regions
    """
    out=[]
    #create chromosome divided intervals
    chromosomes = set([r.chrom for r in regions])
    cdict:dict[str,List[Region]] = {a:[] for a in sorted(chromosomes)}
    for r in regions:
        cdict[r.chrom].append(r)
    for _chrom, c_regions in cdict.items():
        sorted_regions = sorted(c_regions, key=lambda x:x.start)
        out.append(sorted_regions[0])
        for s_r in sorted_regions:
            #check for overlap
            if out[-1].overlaps(s_r):
                #get max region
                out[-1] = Region(out[-1].chrom, min(out[-1].start,s_r.start),max(out[-1].end,s_r.end))
            else:
                out.append(s_r)
    return out