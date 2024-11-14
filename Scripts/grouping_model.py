from annotation_model import Annotation, AnnotationSource
from typing import Any, List, Dict, NamedTuple, Optional, Tuple
from data_access.db import Variant, CS
from enum import Enum, unique
import abc
from functools import cmp_to_key
from time_decorator import timefunc

@unique
class Grouping(Enum):
    NONE = ""
    RANGE = "Range"
    LD = "LD"
    CS = "Credible sets"

@unique
class LDMode(Enum):
    CONSTANT = "Constant"
    DYNAMIC = "Dynamic"




class SummstatColumns(NamedTuple):
    c:str
    p:str
    r:str
    a:str
    pval:str
    beta:str#this is basically included always so let's include it

class GroupingOptions(NamedTuple):
    summary_statistic_path:str
    grouping_mode: Grouping
    column_names: SummstatColumns
    range: int
    ld_mode: LDMode
    r2_threshold: float
    p1_threshold: float
    p2_threshold: float
    overlap: bool

class PhenoInfo(NamedTuple):
    name: Optional[str]
    longname: Optional[str]
    cases: Optional[int]
    controls: Optional[int]
    category: Optional[str]

class Var(NamedTuple):
    id:Variant
    pval: float
    beta: float
    r2_to_lead: Optional[float]

class LocusVariants(NamedTuple):
    lead:Var
    cs:Optional[List[Var]]
    ld_partners:Optional[List[Var]]

class Locus:

    @abc.abstractmethod
    def variants(self) -> List[Variant]:
        pass

    @abc.abstractmethod
    def type(self)->Grouping:
        pass

    @abc.abstractmethod
    def get_vars(self) -> LocusVariants:
        pass

class CSInfo(NamedTuple):
    region:str
    lead:Variant
    number:int
    log10bf:float
    min_r2:float
    size:int
    good_cs:bool


#CSLocus
# attributes: lead_variant, a Variant
# attributes: cs, which is the variants in the credible set. Contains also the lead variant data
# attributes: info, which contains the credible set information that is not included in the cs variants,
# such as region information, log10bf, things like that.
# attributes: ld_partners: Ld partner variants.
class CSLocus(Locus):
    def __init__(self,lead_variant: Var, cs:List[Var],csinfo:CSInfo, ld_partners: Optional[List[Var]]=None):
        self.lead = lead_variant
        self.cs = cs
        self.info=csinfo
        self.ld_partners = ld_partners

    def add_ld_partners(self, ld_partners:List[Var]):
        self.ld_partners = ld_partners
    
    def variants(self) -> List[Variant]:
        output = set()
        output.add(self.lead.id)
        for v in self.cs:
            output.add(v.id)
        if self.ld_partners:
            for v in self.ld_partners:
                output.add(v.id)
        return list(output)

    def type(self)->Grouping:
        return Grouping.CS

    def get_vars(self):
        return LocusVariants(self.lead,self.cs,self.ld_partners)
        



class PeakLocus(Locus):
    def __init__(self, peak: Var, ld_partners: Optional[List[Var]]=None,locus_type=Grouping):
        self.lead  = peak
        self.ld_partners = ld_partners
        self.locus_type = locus_type
    
    def variants(self)->List[Variant]:
        output = set()
        output.add(self.lead.id)
        if self.ld_partners:
            for v in self.ld_partners:
                output.add(v.id)
        return list(output)

    def get_vars(self):
        return LocusVariants(self.lead,None,self.ld_partners)

    def type(self)->Grouping:
        return self.locus_type

def _varcmp(a,b):
    if (a.chrom == b.chrom) and (a.pos==b.pos):
        return 0
    elif (a.chrom < b.chrom) or ((a.chrom == b.chrom) and (a.pos < b.pos)):
        return -1
    return 1


class PhenoData:
    def __init__(self,phenoinfo:PhenoInfo,loci:List[Locus]):
        self.phenoinfo = phenoinfo
        self.loci = loci
        self.annotations: Dict[str,Dict[Variant,Annotation]] = {}
    
    @timefunc
    def get_variants(self)->List[Variant]:
        
        out:set[Variant] = set()
        for l in self.loci:
            out = out.union(l.variants())
        return sorted(list(out),key=cmp_to_key(_varcmp))

    def add_annotation(self,annotation_name,annotation:Dict[Variant,Annotation]):
        self.annotations[annotation_name]=annotation