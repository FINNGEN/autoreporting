"""Concrete implementation of an AlleleDB
"""
import os
import gzip
import pysam
from typing import List, Dict
from functools import partial
from data_access.db import AlleleDB, Location, VariantData

#NOTE: Is this used anywhere?
def _partial_filter(chrom: str, pos: int, cpra: List[str]):
    return (cpra[0] == chrom) and (int(cpra[1]) == pos)
#NOTE: Is this used anywhere?
def remdups(lst):
    seen = set()
    seen_add = seen.add
    return [x for x in lst if not (x in seen or seen_add(x))]

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