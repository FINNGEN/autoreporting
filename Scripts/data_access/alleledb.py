"""Concrete implementation of an AlleleDB
"""
import os
import gzip
import pysam
from typing import List, Dict
from functools import partial
from data_access.db import AlleleDB, Location, VariantData

def _partial_filter(chrom: str, pos: int, cpra: List[str]):
    return (cpra[0] == chrom) and (int(cpra[1]) == pos)

def remdups(lst):
    seen = set()
    seen_add = seen.add
    return [x for x in lst if not (x in seen or seen_add(x))]

class FGAlleleDB(AlleleDB):
    def __init__(self, fg_annotation_file):
        if not os.path.exists(fg_annotation_file):
            raise FileNotFoundError("File {} not found".format(fg_annotation_file))
        try:
            self.handle = pysam.TabixFile(fg_annotation_file)
        except:
            print("Error with pysam handle creation:")
            raise
        with gzip.open(fg_annotation_file, "rt") as f:
            self.header = f.readline().strip("\n").split("\t")
            self.variant_idx = self.header.index("#variant")

    def __delete__(self):
        self.handle.close()

    def get_alleles(self, positions: List[Location])-> List[VariantData]:
        output = []
        for pos in positions:
            rows = self.handle.fetch(pos.chromosome, pos.position-1,pos.position)
            rows = [a.strip('\n').split('\t') for a in rows]
            # different possibilites: 0, 1 or more rows returned
            # 0 means no match
            # 1 means biallelic and match
            # more means not biallelic 
            if len(rows) == 0:
                continue
            elif len(rows) == 1:
                variant = rows[0][self.variant_idx]
                cpra = variant.split(":")
                c,p,r,a = (
                    cpra[0],
                    int(cpra[1]),
                    cpra[2],
                    cpra[3]
                )
                if c == pos.chromosome and p == pos.position:
                    output.append(
                        VariantData(
                            c,
                            p,
                            r,
                            [a],
                            True
                        )
                    )
                else:
                    continue
            else:
                variants = [r[self.variant_idx] for r in rows]
                variants = [v.split(":") for v in variants]
                c,p = (pos.chromosome, pos.position) 
                filtered_variants = list(filter(partial(_partial_filter, c, p), variants))
                ref_alleles = remdups([a[2] for a in filtered_variants])
                alt_alleles = remdups([a[3] for a in filtered_variants])
                #if there are multiple reference alleles for a single variant then things are bad
                if len(ref_alleles) > 1:
                    continue
                output.append(
                    VariantData(
                        filtered_variants[0][0],
                        int(filtered_variants[0][1]),
                        ref_alleles[0],
                        alt_alleles,
                        len(alt_alleles) == 1
                    )
                )
        return output

