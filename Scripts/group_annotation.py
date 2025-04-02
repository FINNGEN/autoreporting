import re,sys
from autoreporting_utils import Region
from data_access.datafactory import db_factory
from typing import NamedTuple, Optional, List, Dict

from grouping_model import CSInfo, PhenoData, Var, Locus
from data_access.db import CSAccess, Variant, CS
from annotation_model import Annotation, AnnotationSource, Value, VariantAnnotation
from load_tabix import tb_resource_manager, TabixResource, TabixOptions
# implement all concrete annotation sources here
# care needs to be taken to make sure that loading an annotation for a variant is a fallible operation, i.e. 
# the annotation might not actually exist for that variant
from time_decorator import timefunc

FUNCTIONAL_CATEGORIES=[
        "transcript_ablation",
        "splice_donor_variant",
        "stop_gained",
        "splice_acceptor_variant",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant"
]


def generate_chrom_ranges(variants:List[Variant])->Dict[str, Dict[str,int]]:
    chr_d = {}
    for v in variants:
        if v.chrom not in chr_d:
            chr_d[v.chrom] = {"start":v.pos-1,"end":v.pos}
        else:
            chr_d[v.chrom]["start"] = max(min(chr_d[v.chrom]["start"],v.pos-1),0)
            chr_d[v.chrom]["end"] = max(chr_d[v.chrom]["end"],v.pos)
    return chr_d


class TabixAnnotation(AnnotationSource):
    def __init__(self,opts: TabixOptions):
        self.cpra = [opts.c,opts.p,opts.r,opts.a]
        self.fname = opts.fname
        self.variant_switch = 100_000
        with tb_resource_manager(self.fname,opts.c,opts.p,opts.r,opts.a) as tb_resource:
            self.sequences= tb_resource.sequences
            if any([a for a in (self.cpra) if a not in tb_resource.header ]):
                raise Exception(f"TabixAnnotation input file did not contain all required columns! Missing columns: {[a for a in self.cpra if a not in tb_resource.header ]}. Supplied columns:{self.cpra}")

    def _chrom_to_source(self, chrom:str)->str:
        """If chromosome needs modification to match tabix source, e.g. adding chr, override this
        """
        return chrom

    def _chrom_from_source(self, chrom:str)->str:
        """If chromoosme needs modification from tabix source, e.g. removing chr, override this
        """
        return chrom

    def _create_annotation(self, cols:List[str], hdi:Dict[str,int])->Annotation:
        """Override this to create annotation
        """
        raise Exception("Override _create_annotation for the annotation")
    
    @timefunc
    def annotate_variants(self, variants: set[Variant]) -> Dict[Variant, Annotation]:
        """Annotate a list of variants with the previous release's p-value and beta
        """
        output = {}
        if len(variants) > self.variant_switch:
            #figure out chroms,starts,ends for each one region in chromosome
            chr_d = generate_chrom_ranges(variants)
            chr_v_sets:Dict[str,set[Variant]] = {key:set() for key in chr_d.keys()}
            for v in variants:
                chr_v_sets[v.chrom].add(v)
            missing_sequences = [self._chrom_to_source(a) for a in chr_d.keys() if self._chrom_to_source(a) not in self.sequences]
            if missing_sequences:
                msg = (f"Warning in annotation tabix file loading: The sequence(s) {missing_sequences} are missing from annotation {self.fname}. List of available sequences in annotation: {[a for a in self.sequences]}.\n"
                        "Skipping the missing sequences.")
                print(msg,file=sys.stderr)
            with tb_resource_manager(self.fname,self.cpra[0],self.cpra[1],self.cpra[2],self.cpra[3]) as tabix_resource:
                hdr = tabix_resource.header
                hdi = {a:i for i,a in enumerate(hdr)}
                for chrom_name, region in chr_d.items():
                    if chrom_name not in self.sequences:
                        print(f"Skipping missing sequence {chrom_name} for annotation {self.fname}",file=sys.stderr)
                        continue
                    chrom_v_set = chr_v_sets[chrom_name]
                    for l in tabix_resource.fileobject.fetch(self._chrom_to_source(chrom_name),region["start"],region["end"]):
                        cols = l.split("\t")
                        v = Variant(self._chrom_from_source(cols[hdi[self.cpra[0]]]),int(cols[hdi[self.cpra[1]]]),cols[hdi[self.cpra[2]]],cols[hdi[self.cpra[3]]] )
                        if v in chrom_v_set:
                            output[v] = self._create_annotation(cols,hdi)
        else:
            with tb_resource_manager(self.fname,self.cpra[0],self.cpra[1],self.cpra[2],self.cpra[3]) as tabix_resource:
                hdr = tabix_resource.header
                hdi = {a:i for i,a in enumerate(hdr)}
                sequences_in_variants = set([a.chrom for a in variants])
                missing_sequences = [self._chrom_to_source(a) for a in sequences_in_variants if self._chrom_to_source(a) not in self.sequences]
                if missing_sequences:
                    msg = (f"Warning in annotation tabix file loading: The sequence(s) {missing_sequences} are missing from annotation {self.fname}. List of available sequences in annotation: {[a for a in self.sequences]}.\n"
                            "Skipping the missing sequences.")
                    print(msg,file=sys.stderr)
                for original_v in variants:
                    if original_v.chrom not in self.sequences:
                        continue
                    for l in tabix_resource.fileobject.fetch(self._chrom_to_source(original_v.chrom),max(original_v.pos-1,0),original_v.pos):
                        cols = l.split("\t")
                        v = Variant(self._chrom_from_source(cols[hdi[self.cpra[0]]]),int(cols[hdi[self.cpra[1]]]),cols[hdi[self.cpra[2]]],cols[hdi[self.cpra[3]]] )
                        if v in variants:
                            output[v] = self._create_annotation(cols,hdi)
        return output
    @staticmethod
    def get_name() -> str:
        raise Exception("annotation name not overridden for generic class TabixAnnotation")

class PreviousReleaseOptions(NamedTuple):
    fname: str
    c: str
    p: str
    r: str
    a: str
    pval_col: str
    beta_col: str


class PreviousReleaseAnnotation(TabixAnnotation):
    def __init__(self, opts: PreviousReleaseOptions) -> None:
        super().__init__(TabixOptions(opts.fname,opts.c,opts.p,opts.r,opts.a))
        self.pval_col = opts.pval_col
        self.beta_col = opts.beta_col
        self.variant_switch = 50_000
        #validate source file
        with tb_resource_manager(self.fname,opts.c,opts.p,opts.r,opts.a) as tb_resource:
            if any([a for a in [self.pval_col,self.beta_col] if a not in tb_resource.header ]):
                raise Exception(f"Previous release annotation file did not contain all required columns! Missing columns: {[a for a in [self.pval_col,self.beta_col] if a not in tb_resource.header ]}. Supplied columns:{[self.pval_col,self.beta_col]}")
        

    def _create_annotation(self,cols:List[str],hdi:Dict[str,int])->Annotation:
        #get pval and beta
        pval = tryfloat(cols[hdi[self.pval_col]])
        beta = tryfloat(cols[hdi[self.beta_col]])
        return [{"pval_previous_release":pval,"beta_previous_release":beta}]

    @staticmethod
    def get_name():
        return "previous_release"

    @staticmethod
    def get_output_columns():
        return ["pval_previous_release","beta_previous_release"]

class ExtraColAnnotation(TabixAnnotation):
    def __init__(self, opts: TabixOptions,extra_columns:List[str]) -> None:
        super().__init__(opts)
        self.extra_columns = extra_columns
        self.variant_switch = 50_000
        with tb_resource_manager(self.fname,opts.c,opts.p,opts.r,opts.a) as tb_resource:
            if any([a for a in extra_columns if a not in tb_resource.header ]):
                raise Exception(f"Summary statistic file did not contain all extra columns! Missing columns: {[a for a in self.extra_columns if a not in tb_resource.header ]}. Supplied columns:{self.extra_columns}")
        

    def _create_annotation(self,cols:List[str],hdi:Dict[str,int])->Annotation:
        annotation:VariantAnnotation = {a:cols[hdi[a]] for a in self.extra_columns}
        return [annotation]

    @staticmethod
    def get_name():
        return "extra_columns"


class CSAnnotation(AnnotationSource):
    def __init__(self, cs_access: CSAccess):
        self.cs_access = cs_access
        csdata = self.cs_access.get_cs()
        #CS can be indexed with region-number
        self.cs_dict:Dict[str,CS] = {f"{c.region}-{c.number}":c for c in csdata}
        self.variant_dict:Dict[Variant,str] = {}
        for csid,cs in self.cs_dict.items():
            self.variant_dict[cs.lead] = csid
            for v in cs.variants:
                self.variant_dict[v.variant] = csid
        
    @timefunc
    def annotate_variants(self, variants: set[Variant]) -> Dict[Variant, Annotation]:
        """Annotation has the following pieces of data:
        region:
        cs number:
        log10bf:
        min r2:
        size:
        good cs:
        lead var:
        prob:
        lead r2:
        """
        
        output = {}
        for v in variants:
            var_annot = []
            cs:Optional[CS] = self.cs_dict.get(self.variant_dict.get(v,None),None)
            if cs is not None:
                annotation:Dict[str,Value] = {}
                annotation["region"] = cs.region
                annotation["lead_variant"] = cs.lead
                annotation["number"] = cs.number
                annotation["bayes"] = cs.bayes
                annotation["min_r2"] = cs.min_r2
                annotation["size"] = cs.size
                annotation["good_cs"] = cs.good_cs
                specific_var = [a for a in cs.variants if a.variant == v][0]
                annotation["prob"] = specific_var.prob
                annotation["lead_r2"] = specific_var.lead_r2
                var_annot.append(annotation)
            if var_annot:
                output[v] = var_annot

        return output

    def get_cs_info(self, lead_variant:Variant)->CSInfo:
        correct_cs = [a for a in self.csdata if a.lead == lead_variant][0]
        return CSInfo(
            correct_cs.region,
            correct_cs.lead,
            correct_cs.number,
            correct_cs.bayes,
            correct_cs.min_r2,
            correct_cs.size,
            correct_cs.good_cs
        )
        
    @staticmethod
    def get_name():
        return "cs"

    @staticmethod
    def get_output_columns() -> List[str]:
        return [
            "good_cs",
            "region",
            "lead_variant",
            "number",
            "bayes",
            "min_r2",
            "size",
            "prob",
            "lead_r2"
        ]

def tryfloat(value: str) -> float:
    return float(value) if value != "NA" else float("nan")

def tryint(value: str) -> Optional[int]:
    return int(value) if value !="NA" else None

class FunctionalAnnotation(TabixAnnotation):
    def __init__(self,fname:str):
        opts = TabixOptions(fname,"chrom","pos","ref","alt")
        super().__init__(opts)
        self.columntypes = {
            "enrichment_nfsee":tryfloat,
            "fin.AF":tryfloat,
            "fin.AN":tryint,
            "fin.AC":tryint,
            "fin.homozygote_count":tryint,
            "fet_nfsee.odds_ratio":tryfloat,
            "fet_nfsee.p_value":tryfloat,
            "nfsee.AC":tryint,
            "nfsee.AN":tryint,
            "nfsee.AF":tryfloat,
            "nfsee.homozygote_count":tryint
        }
        self.variant_switch = 1_000_000
        #make sure all columns are in 
        with tb_resource_manager(self.fname,opts.c,opts.p,opts.r,opts.a) as tb_resource:
            if any([a for a in [b for b in self.columntypes.keys()] if a not in tb_resource.header ]):
                raise Exception(f"Functional annotation file did not contain all required columns! Missing columns: {[a for a in [b for b in self.columntypes.keys()] if a not in tb_resource.header ]}. Supplied columns:{[b for b in self.columntypes.keys()]}")

    def _chrom_to_source(self,chrom:str)->str:
        """If chromosome needs modification to match tabix source, e.g. adding chr, override this
        """
        return "chr"+chrom.replace("23","X").replace("23","Y")

    def _chrom_from_source(self,chrom:str)->str:
        """If chromosome needs modification from tabix source, e.g. removing chr, override this
        """
        return chrom.replace("chr","").replace("X","23").replace("Y","24")


    def _create_annotation(self,cols:List[str], hdi: Dict[str,int])->Annotation:
        try:
            annotation:Dict[str,Value] = {key:valuetype(cols[hdi[key]]) for key, valuetype in self.columntypes.items()}
        except:
            print(cols)
            raise
        return [annotation]

    @staticmethod
    def get_name():
        return "functional_annotation"

    @staticmethod
    def get_output_columns() -> List[str]:
        return [
            "enrichment_nfsee",
            "fin.AF",
            "fin.AN",
            "fin.AC",
            "fin.homozygote_count",
            "fet_nfsee.odds_ratio",
            "fet_nfsee.p_value",
            "nfsee.AC",
            "nfsee.AN",
            "nfsee.AF",
            "nfsee.homozygote_count"
        ]

class FGAnnotation(TabixAnnotation):
    def __init__(self,fname:str):
        super().__init__(TabixOptions(fname,"chr","pos","ref","alt"))
        self.out_columns = [
            "most_severe_gene",
            "most_severe_consequence",
            "FG_INFO",
            "n_INFO_gt_0_6",
            "functional_category",
            "rsids"
        ]
        self.colnames = {
            "most_severe_gene":"gene_most_severe",
            "most_severe_consequence":"most_severe",
            "FG_INFO":"INFO",
            "rsids":"rsid"
        }
        self.variant_switch = 1_000_000

    def _create_annotation(self, cols: List[str], hdi: Dict[str, int])->Annotation:
        most_severe_gene = cols[hdi[self.colnames["most_severe_gene"]]]
        most_severe_consequence = cols[hdi[self.colnames["most_severe_consequence"]]]
        functional_category = cols[hdi[self.colnames["most_severe_consequence"]]] if cols[hdi[self.colnames["most_severe_consequence"]]] in FUNCTIONAL_CATEGORIES else "NA"
        fg_info = cols[hdi[self.colnames["FG_INFO"]]]
        rsid = cols[hdi[self.colnames["rsids"]]]
        infocol_names = [a for a in hdi.keys() if a.startswith("INFO_")]
        infocol_values = [tryfloat(cols[hdi[a]]) for a in infocol_names]
        n_info_gt_0_6 = len([a for a in infocol_values if a > 0.6])/len(infocol_values)
        return [{
            self.out_columns[0]:most_severe_gene,
            self.out_columns[1]:most_severe_consequence,
            self.out_columns[2]:fg_info,
            self.out_columns[3]:n_info_gt_0_6,
            self.out_columns[4]:functional_category,
            self.out_columns[5]:rsid
        }]

    @staticmethod
    def get_name() -> str:
        return "finngen"

    @staticmethod
    def get_output_columns() -> List[str]:
        return ["most_severe_gene",
            "most_severe_consequence",
            "FG_INFO",
            "n_INFO_gt_0_6",
            "functional_category",
            "rsids"
        ]

def calculate_enrichment(nfe_AC:List[int],nfe_AN:List[int],fi_af:float):
    """Calculate enrichment, clipped to range [0.0, 1e6]
    In case fi_af and AC are both 0, return nan.
    If only AC is 0, return 1_000_000.
    Else return calculated enrichment, clipped to range [0.0, 1e6].
    """
    AC = sum(nfe_AC)
    if AC == 0 and fi_af == 0.0:
        return float("nan")
    elif AC == 0:
        return 1_000_000.0
    return max(min(sum(nfe_AN)*fi_af/AC,1_000_000.0),0.0)

class GnomadGenomeAnnotation(TabixAnnotation):
    def __init__(self,fname: str):
        super().__init__(TabixOptions(fname,"#CHROM","POS","REF","ALT"))
        self.variant_switch = 1_000_000
        self.columns=["AF_fin",
        "AF_nfe",
        "AF_nfe_est",
        "AF_nfe_nwe",
        "AF_nfe_onf",
        "AF_nfe_seu"]
        self.nfe_ac = ["AC_nfe_est","AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
        self.nfe_an = ["AN_nfe_est","AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
        self.nfe_est_ac = ["AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
        self.nfe_est_an = ["AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
        self.af_fin = "AF_fin"
        self.colrename = lambda x: f"GENOME_{x}"
        with tb_resource_manager(self.fname,self.cpra[0],self.cpra[1],self.cpra[2],self.cpra[3]) as res:
            hdr = res.header
            colnames = list(set(self.nfe_ac + self.nfe_an + self.nfe_est_ac + self.nfe_est_an + [self.af_fin]))
            if any([a for a in colnames if a not in hdr]):
                raise Exception(f"columns {[a for a in colnames if a not in hdr]} not found in gnomad genome annotation header! Is this the correct file for gnomad genome annotation?" )

    def _create_annotation(self, cols: List[str], hdi: Dict[str, int])->Annotation:
        nfe_ac_values = [int(cols[hdi[a]]) for a in self.nfe_ac]
        nfe_an_values = [int(cols[hdi[a]]) for a in self.nfe_an]
        nfe_est_ac_values = [int(cols[hdi[a]]) for a in self.nfe_est_ac]
        nfe_est_an_values = [int(cols[hdi[a]]) for a in self.nfe_est_an]
        af_fin = float(cols[hdi[self.af_fin]])
        enrichment_nfe = calculate_enrichment(nfe_ac_values,nfe_an_values,af_fin)
        enrichment_nfe_est = calculate_enrichment(nfe_est_ac_values,nfe_est_an_values,af_fin)
        #fill in values to annotation
        annotation:Dict[str,Value] = { self.colrename(a):tryfloat(cols[hdi[a]]) for a in self.columns }
        annotation["GENOME_FI_enrichment_nfe"] = enrichment_nfe
        annotation["GENOME_FI_enrichment_nfe_est"] = enrichment_nfe_est
        return [annotation]

    @staticmethod
    def get_name():
        return "gnomad_genome"

    @staticmethod
    def get_output_columns() -> List[str]:
        return ["GENOME_AF_fin",
        "GENOME_AF_nfe",
        "GENOME_AF_nfe_est",
        "GENOME_AF_nfe_nwe",
        "GENOME_AF_nfe_onf",
        "GENOME_AF_nfe_seu",
        "GENOME_FI_enrichment_nfe",
        "GENOME_FI_enrichment_nfe_est"]

class GnomadExomeAnnotation(TabixAnnotation):
    def __init__(self,fname: str):
        super().__init__(TabixOptions(fname,"#CHROM","POS","REF","ALT"))
        self.variant_switch = 50_000
        self.columns=["AF_nfe_bgr",
        "AF_fin",
        "AF_nfe",
        "AF_nfe_est",
        "AF_nfe_swe",
        "AF_nfe_nwe",
        "AF_nfe_onf",
        "AF_nfe_seu"]
        self.nfe_ac=["AC_nfe_bgr","AC_nfe_est","AC_nfe_onf","AC_nfe_seu","AC_nfe_swe"]
        self.nfe_an=["AN_nfe_bgr","AN_nfe_est","AN_nfe_onf","AN_nfe_seu","AN_nfe_swe"]
        self.nfe_est_ac=["AC_nfe_bgr","AC_nfe_onf","AC_nfe_seu","AC_nfe_swe"]
        self.nfe_est_an=["AN_nfe_bgr","AN_nfe_onf","AN_nfe_seu","AN_nfe_swe"]
        self.nfe_swe_ac=["AC_nfe_bgr","AC_nfe_est","AC_nfe_onf","AC_nfe_seu"]
        self.nfe_swe_an=["AN_nfe_bgr","AN_nfe_est","AN_nfe_onf","AN_nfe_seu"]
        self.nfe_est_swe_ac=["AC_nfe_bgr","AC_nfe_onf","AC_nfe_seu"]
        self.nfe_est_swe_an=["AN_nfe_bgr","AN_nfe_onf","AN_nfe_seu"]
        self.af_fin = "AF_fin"
        self.colrename = lambda x: f"EXOME_{x}"

    def _create_annotation(self, cols: List[str], hdi: Dict[str, int])->Annotation:
        nfe_ac_values = [int(cols[hdi[a]]) for a in self.nfe_ac]
        nfe_an_values = [int(cols[hdi[a]]) for a in self.nfe_an]
        nfe_est_ac_values = [int(cols[hdi[a]]) for a in self.nfe_est_ac]
        nfe_est_an_values = [int(cols[hdi[a]]) for a in self.nfe_est_an]
        nfe_swe_ac_values = [int(cols[hdi[a]]) for a in self.nfe_swe_ac]
        nfe_swe_an_values = [int(cols[hdi[a]]) for a in self.nfe_swe_an]
        nfe_est_swe_ac_values = [int(cols[hdi[a]]) for a in self.nfe_est_swe_ac]
        nfe_est_swe_an_values = [int(cols[hdi[a]]) for a in self.nfe_est_swe_an]
        af_fin = tryfloat(cols[hdi[self.af_fin]])
        enrichment_nfe = calculate_enrichment(nfe_ac_values,nfe_an_values,af_fin)
        enrichment_nfe_est = calculate_enrichment(nfe_est_ac_values,nfe_est_an_values,af_fin)
        enrichment_nfe_swe = calculate_enrichment(nfe_swe_ac_values,nfe_swe_an_values,af_fin)
        enrichment_nfe_est_swe = calculate_enrichment(nfe_est_swe_ac_values,nfe_est_swe_an_values,af_fin)
        #fill in values to annotation
        annotation:Dict[str,Value] = {self.colrename(a):tryfloat(cols[hdi[a]]) for a in self.columns }
        annotation["EXOME_FI_enrichment_nfe"] = enrichment_nfe
        annotation["EXOME_FI_enrichment_nfe_est"] = enrichment_nfe_est
        annotation["EXOME_FI_enrichment_nfe_swe"] = enrichment_nfe_swe
        annotation["EXOME_FI_enrichment_nfe_est_swe"] = enrichment_nfe_est_swe
        return [annotation]

    @staticmethod
    def get_name():
        return "gnomad_exome"
    
    @staticmethod
    def get_output_columns() -> List[str]:
        return [
            "EXOME_AF_nfe_bgr",
            "EXOME_AF_fin",
            "EXOME_AF_nfe",
            "EXOME_AF_nfe_est",
            "EXOME_AF_nfe_swe",
            "EXOME_AF_nfe_nwe",
            "EXOME_AF_nfe_onf",
            "EXOME_AF_nfe_seu",
            "EXOME_FI_enrichment_nfe",
            "EXOME_FI_enrichment_nfe_est",
            "EXOME_FI_enrichment_nfe_swe",
            "EXOME_FI_enrichment_nfe_est_swe",
        ]

class Gnomad4Annotation(TabixAnnotation):
    def __init__(self,fname: str):
        super().__init__(TabixOptions(fname,"#chr","pos","ref","alt"))
        self.variant_switch = 500_000


    def _create_annotation(self, cols: List[str], hdi: Dict[str, int]) -> Annotation:
        annotation:Dict[str,Value] = {}
        annotation["GNOMAD_AF_fin"] = tryfloat(cols[hdi["AF_fin"]])
        annotation["GNOMAD_AF_nfe"] = tryfloat(cols[hdi["AF_nfe"]])
        annotation["GNOMAD_FI_enrichment_nfe"] = tryfloat(cols[hdi["enrichment_nfe"]])
        return [annotation]

    def _chrom_to_source(self,chrom:str)->str:
        """If chromosome needs modification to match tabix source, e.g. adding chr, override this
        """
        return chrom.replace("23","X").replace("24","Y")

    def _chrom_from_source(self,chrom:str)->str:
        """If chromosome needs modification from tabix source, e.g. removing chr, override this
        """
        return chrom.replace("X","23").replace("Y","24")

    @staticmethod
    def get_name():
        return "gnomad_4"

    @staticmethod
    def get_output_columns() -> List[str]:
        return ["GNOMAD_AF_fin",
        "GNOMAD_AF_nfe",
        "GNOMAD_FI_enrichment_nfe"]

def map_alleles(a1,a2):
    """
    Flips alleles to the A strand if necessary and orders them lexicogaphically
    Author: Pietro, with small changes by Arto
    """
    allele_dict={"T":"A","C":"G","G":"C"}
     # check if the A variant is present
    if 'A' not in a1 + a2 and 'a' not in a1+a2:
        # for both/ref and alt map them to the A strand and order each one lexicographically
        a1 = ''.join([allele_dict[elem.upper()] for elem in a1])
        a2 = ''.join([allele_dict[elem.upper()] for elem in a2])
    # further sorting
    return sorted([a1,a2])

def construct_mapping_variant(v:Variant)->Variant:
    try:
        mapped_alleles = map_alleles(v.ref,v.alt)
        return Variant(v.chrom,v.pos,mapped_alleles[0],mapped_alleles[1])
    except:
        return v

class CatalogAnnotation(AnnotationSource):
    def __init__(self,use_gwascatalog, custom_dataresource, database_choice, localdb_path, gwas_width, gwas_pval ,gwas_threads, allele_filepath):
        self.inner = db_factory(use_gwascatalog, custom_dataresource, database_choice, localdb_path, gwas_width, gwas_pval ,gwas_threads, allele_filepath)
        self.cols = [
            "#variant_hit","pval","trait","trait_name","study_link"
        ]
        self.col_rename = {a:a for a in self.cols}
        self.col_rename["pval"] = "pval_trait"

    @timefunc
    def annotate_variants(self,variants:List[Variant]) -> Dict[Variant,Annotation]:
        #generate chromosomal regions, that might be the best for now.
        # The gwas catalog/custom catalog data is not so big, so it will greatly limit 
        # the amount of results coming out.
        # In the long term, this abomination should be replaced with a sane setup, without pandas.
        # TODO: variant matching etc, and the following rigamarole
        #what columns are output?
        #trait and traitname, pval_trait, study_link
        chr_d = generate_chrom_ranges(variants)
        regs = [Region(a,b["start"],b["end"]) for a,b in chr_d.items()]

        data=self.inner.associations_for_regions(regs)
        #format data according to variant
        #variants are identified by "chrom","pos","ref","alt"
        output:Dict[Variant,Annotation] = {}

        mapping_variant_set = {construct_mapping_variant(v):v for v in variants}
        mset='^[acgtACGT-]+$'
        reg = re.compile(mset)
        for item in data:
            #filter out invalid alleles
            if not (reg.match(item["ref"]) and reg.match(item["alt"])):
                continue
            v = Variant(item["chrom"],item["pos"],item["ref"],item["alt"])
            mapping_v = construct_mapping_variant( v )

            if mapping_v in mapping_variant_set:
                item["#variant_hit"] = f"chr{item['chrom']}_{item['pos']}_{item['ref']}_{item['alt']}"
                entry = {self.col_rename[a]:item[a] for a in self.cols}
                if mapping_variant_set[mapping_v] in output:
                    output[mapping_variant_set[mapping_v]].append(entry)
                else:
                    output[mapping_variant_set[mapping_v]] = [entry]
        return output

    @staticmethod
    def get_name():
        return "compound_catalog"

    @staticmethod
    def get_output_columns() -> List[str]:
        return  [
            "#variant_hit","pval_trait","trait","trait_name","study_link"
        ]

def annotate(pheno:PhenoData,annotations:List[AnnotationSource]) -> PhenoData:
    vars = pheno.get_variants()
    for asource in annotations:
        print(f"annotating variants with {asource.get_name()}")
        annotation = asource.annotate_variants(vars)
        annotation_name = asource.get_name()
        pheno.add_annotation(annotation_name,annotation)
    return pheno
