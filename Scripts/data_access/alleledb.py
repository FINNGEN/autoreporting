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
"""
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

    def __del__(self):
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
                #if there are multiple reference alleles for a single variant then things are bad, since we can't identify the variants
                if len(ref_alleles) > 1:
                    continue
                output.append(
                    VariantData(
                        filtered_variants[0][0],
                        int(filtered_variants[0][1]),
                        ref_alleles[0],
                        alt_alleles,
                        len(alt_alleles) == 1,
                        ""
                    )
                )
        return output
"""

class VCFAlleleDB(AlleleDB):
    """Allele db based on vcf file from dbsnp. NOTE: File can be on ftp server, IF the index file exists. 
    E.g. https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
    """

    def __init__(self, vcf_file):
        if not os.path.exists(vcf_file):
            raise FileNotFoundError("File {} not found".format(vcf_file))
        self.handle = None
        try:
            self.handle = pysam.TabixFile(vcf_file)
        except:
            print("Error with pysam handle creation:")
            raise
        self.file = vcf_file
        self.header = self.handle.header[-1].strip('\n').split('\t')


    def get_alleles(self, positions: List[Location])-> List[VariantData]:
        output = []
        for pos in positions:
            #vcf is numeric
            chrom = pos.chromosome.replace("X","23").replace("Y","24").replace("MT","25").replace("M","25")
            rows = self.handle.fetch(chrom, pos.position-1,pos.position)
            rows = [a.strip('\n').split('\t') for a in rows]
            rows = [a for a in rows if int(a[1]) == pos.position]#indels behave weirdly when using tabix

            #get alleles from rows
            for r in rows:
                output.append(
                    VariantData(
                        r[0],
                        int(r[1]),
                        r[3],
                        r[4].split(','),
                        len(r[4].split(',')) == 1,
                        int(r[2].replace('rs',''))
                    )
                )
        return output

    def __del__(self):
        self.handle.close()