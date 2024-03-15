from annotation_model import Value
from group_annotation import CSAnnotation, CatalogAnnotation, ExtraColAnnotation, FGAnnotation, FunctionalAnnotation, GnomadExomeAnnotation, GnomadGenomeAnnotation, PreviousReleaseAnnotation, Gnomad4Annotation
from grouping_model import CSInfo, CSLocus, Grouping, LDMode, PeakLocus, PhenoData, PhenoInfo, SummstatColumns, Var
from grouping import ld_threshold
from data_access.db import Variant
from typing import Any, Dict, Optional, TextIO, NamedTuple, List, Tuple, TypeVar, cast
from collections import defaultdict
from copy import deepcopy
import scipy.stats as stats # type: ignore
import math
from time_decorator import timefunc

class VariantReportOptions(NamedTuple):
    summstat_options: SummstatColumns
    extra_columns: List[str]
    squash_multiple_annotations:bool

class TopReportOptions(NamedTuple):
    strict_group:float
    significance_threshold:float
    grouping_method: Grouping
    efo_traits: List[str]
    extra_columns: List[str]
    ld_mode: LDMode
    r2_threshold: float
    


def format_value(value:Value)->str:
    if value is None:
        return "NA"
    elif type(value) == int:
        return str(value)
    elif type(value) == float:
        if math.isnan(value):
            return "NA"
        return f"{value:.6g}"
    elif type(value) == str:
        return value
    elif type(value) == bool:
        return str(value)
    elif type(value) == Variant:
        return f"chr{value.chrom}_{value.pos}_{value.ref}_{value.alt}"
    else:
        raise Exception(f"unsupported type of value: {type(value)} for value {value}")
@timefunc
def generate_variant_report(data:PhenoData,output:TextIO, options: VariantReportOptions):
    """Generate and write the variant report from available data.
    Inputs:
        data: PhenoData is the grouped and annotated data
        output: TextIO is where we write outputs
        options: VariantReportOptions contains all relevant meta-information for e.g. column names etc.
    """
    ### List all columns that are created
    annotation_columns = [
        a for colsource in 
            (
                #GnomadGenomeAnnotation.get_output_columns(),
                #GnomadExomeAnnotation.get_output_columns(),
                Gnomad4Annotation.get_output_columns(),
                FunctionalAnnotation.get_output_columns(),
                FGAnnotation.get_output_columns(),
                PreviousReleaseAnnotation.get_output_columns(),
                CatalogAnnotation.get_output_columns()
            ) 
        for a in colsource
    ]
    outputcolumns = [
        options.summstat_options.c,
        options.summstat_options.p,
        options.summstat_options.r,
        options.summstat_options.a,
        options.summstat_options.pval,
        options.summstat_options.beta,
        "r2_to_lead"] + options.extra_columns+[
        "cs_id",
        "cs_region",
        "cs_number",
        "cs_prob",
        "cs_log10bf",
        "cs_min_r2",
        "cs_size",
        "good_cs",
        "#variant",
        "locus_id",
        "pos_rmax",
        "pos_rmin",
        "phenotype",
        "longname",
        "category",
        "n_cases",
        "n_controls"] + annotation_columns
    #print header
    header = "\t".join(outputcolumns)+"\n"
    # whether there is cs data
    cs_available = False
    if CSAnnotation.get_name() in data.annotations:
        cs_available = True
    output.write(header)
    ### For each locus, and each variant in them, calculate values & write output
    for locus in data.loci:
        vars:List[Var] = []
        locus_v_lists = locus.get_vars()
        lead = locus_v_lists.lead
        ld_partners = locus_v_lists.ld_partners if locus_v_lists.ld_partners is not None else []
        cs_vars = locus_v_lists.cs if locus_v_lists.cs is not None else []
        vars.extend(ld_partners)
        vars.extend(cs_vars)
        vars.append(lead)
        # remove duplicates
        vars = sorted(list(set(vars)),key = lambda x: x.pval)
        #precalculate locus-wide values
        pos_rmin = min([a.id.pos for a in vars])
        pos_rmax = max([a.id.pos for a in vars])
        locus_id = f"chr{lead.id.chrom}_{lead.id.pos}_{lead.id.ref}_{lead.id.alt}"

        for v in vars:
            list_of_cols = []
            #init cols with by default no value for a column
            cols:Dict[str,Value] = defaultdict(lambda:None)
            #variant cols
            cols[options.summstat_options.c] = v.id.chrom
            cols[options.summstat_options.p] = v.id.pos
            cols[options.summstat_options.r] = v.id.ref
            cols[options.summstat_options.a] = v.id.alt
            cols[options.summstat_options.pval] = v.pval
            cols[options.summstat_options.beta] = v.beta
            cols["r2_to_lead"] = v.r2_to_lead

            #cs
            #two ways of working:
            # in case of a PeakLocus, we just add the cs data, which is simply defined
            # if it's a cslocus, we find the corresponding cs (lead is cs lead), and all the other cs
            # if this variant is part of cs, add all cs data
            # else add only cs_number
            if cs_available:
                if isinstance(locus,CSLocus):
                    #check if variant is in cs
                    annot = data.annotations[CSAnnotation.get_name()]
                    if v in set(cs_vars+[lead]):
                        try:
                            cs_data = [a for a in annot[v.id] if a["lead_variant"]==lead.id][0]
                        except:
                            raise Exception(f"Variant {v.id} is in credible set but does not have cs data in annotation. Check implementation.")
                        cs_lead = cs_data["lead_variant"] if isinstance(cs_data["lead_variant"],Variant) else None
                        cs_number = cs_data["number"] if isinstance(cs_data["number"],int) else None
                        cols["cs_id"] = f"chr{cs_lead.chrom}_{cs_lead.pos}_{cs_lead.ref}_{cs_lead.alt}_{cs_number}" if isinstance(cs_lead,Variant) else None
                        cols["cs_region"] = cs_data["region"]
                        cols["cs_number"] = cs_data["number"]
                        cols["cs_prob"] = cs_data["prob"]
                        cols["cs_log10bf"] = cs_data["bayes"]
                        cols["cs_min_r2"] = cs_data["min_r2"]
                        cols["cs_size"] = cs_data["size"]
                        cols["good_cs"] = cs_data["good_cs"]
                    else:
                        #if the variant has cs annotation, only set cs_number
                        if v.id in annot:
                            cols["cs_number"] = annot[v.id][0]["number"]
                else:
                    #fetch cs annotations for each variant, no checking for which cs it belongs to.
                    annot = data.annotations[CSAnnotation.get_name()]
                    if v.id in annot:
                        #NOTE: this picks the largest prob cs for this variant. In ideal circumstances, a variant
                        # would always have only one cs which it belongs to, but 
                        # in practice this is not always enforced in region selection. 
                        cs_data = sorted(annot[v.id],key = lambda x:float(x["prob"]))[-1] #type:ignore
                        cs_lead = cs_data["lead_variant"] if isinstance(cs_data["lead_variant"],Variant) else None
                        cs_number = cs_data["number"] if isinstance(cs_data["number"],int) else None
                        cols["cs_id"] = f"chr{cs_lead.chrom}_{cs_lead.pos}_{cs_lead.ref}_{cs_lead.alt}_{cs_number}" if isinstance(cs_lead,Variant) else None
                        cols["cs_region"] = cs_data["region"]
                        cols["cs_number"] = cs_data["number"]
                        cols["cs_prob"] = cs_data["prob"]
                        cols["cs_log10bf"] = cs_data["bayes"]
                        cols["cs_min_r2"] = cs_data["min_r2"]
                        cols["cs_size"] = cs_data["size"]
                        cols["good_cs"] = cs_data["good_cs"]
            #locus specs
            cols["#variant"] = f"chr{v.id.chrom}_{v.id.pos}_{v.id.ref}_{v.id.alt}" 
            cols["pos_rmin"] = pos_rmin
            cols["pos_rmax"] = pos_rmax
            cols["locus_id"] = locus_id
            #phenoinfo
            if data.phenoinfo:
                cols["phenotype"] = data.phenoinfo.name
                cols["longname"] = data.phenoinfo.longname
                cols["category"] = data.phenoinfo.category
                cols["n_cases"] = data.phenoinfo.cases
                cols["n_controls"] = data.phenoinfo.controls
            
            ### Annotations

            #single annotations
            single_annotation_dict = {
                #GnomadGenomeAnnotation.get_name():GnomadGenomeAnnotation.get_output_columns(),
                #GnomadExomeAnnotation.get_name():GnomadExomeAnnotation.get_output_columns(),
                Gnomad4Annotation.get_name():Gnomad4Annotation.get_output_columns(),
                FunctionalAnnotation.get_name():FunctionalAnnotation.get_output_columns(),
                FGAnnotation.get_name():FGAnnotation.get_output_columns(),
                PreviousReleaseAnnotation.get_name():PreviousReleaseAnnotation.get_output_columns(),
                ExtraColAnnotation.get_name():options.extra_columns,
            }
            for ann_name, ann_columns in single_annotation_dict.items():
                if ann_name in data.annotations:
                    if v.id in data.annotations[ann_name]:
                        for c in ann_columns:
                            cols[c] = data.annotations[ann_name][v.id][0][c]
            #multiple annotation, everything else is already ready
            #catalog annotation, squash if necessary, else we might have to keep a variable for lines in hand.
            ann_name = CatalogAnnotation.get_name()
            ann_columns = CatalogAnnotation.get_output_columns()
            if ann_name in data.annotations:
                anno = data.annotations[ann_name]
                if v.id in anno:
                    #get length
                    for annotation in anno.get(v.id,[]):
                        nucols = deepcopy(cols)
                        for c in ann_columns:
                            nucols[c] = annotation[c]
                        list_of_cols.append(nucols)
                else:
                    #only emit one row
                    list_of_cols.append(cols)
            else:
                #only emit one row
                list_of_cols.append(cols)
            #output the data
            #squashing not supported yet
            if options.squash_multiple_annotations:
                raise NotImplementedError
            else:
                for columns in list_of_cols:
                    column_values = [format_value(columns[colname]) for colname in outputcolumns]
                    output.write("\t".join(column_values)+"\n")

get_varid = lambda variant: f"chr{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"

@timefunc
def generate_top_report(data:PhenoData,output:TextIO, options: TopReportOptions):
    """Generate and write the top report (group aggregation report) from available data.
    Inputs:
        data: PhenoData is the grouped and annotated data
        output: TextIO is where we write outputs
        options: TopReportOptions contains all relevant meta-information for e.g. column names etc.
    """
    ### List all columns that are created
    lead_cols = [
        "lead_r2_threshold",
        "lead_beta_previous_release",
        "lead_pval_previous_release",
        "lead_most_severe_consequence",
        "lead_most_severe_gene",
        "lead_enrichment",
        "lead_mlogp",
        "lead_beta",
        "lead_sebeta",
        "lead_af_alt",
        "lead_af_alt_cases",
        "lead_af_alt_controls"
    ] 
    if ExtraColAnnotation.get_name() in data.annotations:
        for c in options.extra_columns:
            lead_c_name = f"lead_{c}"
            if lead_c_name not in lead_cols:
                lead_cols.append(lead_c_name)
    columns = [
        "phenotype",
        "phenotype_abbreviation",
        "locus_id",
        "rsids",
        "Cases",
        "Controls",
        "chrom",
        "pos",
        "ref",
        "alt",
        "pval"] + lead_cols + [
        "gnomAD_functional_category",
        "gnomAD_enrichment_nfsee",
        "gnomAD_fin.AF",
        "gnomAD_fin.AN",
        "gnomAD_fin.AC",
        "gnomAD_fin.homozygote_count",
        "gnomAD_fet_nfsee.odds_ratio",
        "gnomAD_fet_nfsee.p_value",
        "gnomAD_nfsee.AC",
        "gnomAD_nfsee.AN",
        "gnomAD_nfsee.AF",
        "gnomAD_nfsee.homozygote_count",
        "cs_id",
        "cs_size",
        "cs_log_bayes_factor",
        "cs_number",
        "cs_region",
        "good_cs",
        "credible_set_min_r2_value",
        "start",
        "end",
        "found_associations_strict",
        "found_associations_relaxed",
        "credible_set_variants",
        "functional_variants_strict",
        "functional_variants_relaxed",
        "specific_efo_trait_associations_strict",
        "specific_efo_trait_associations_relaxed",
        "n_ld_partners_0_8",
        "n_ld_partners_0_6"
    ]
    #phenoinfo
    phenotype = data.phenoinfo.longname
    pheno_abbr = data.phenoinfo.name
    cases = data.phenoinfo.cases
    controls = data.phenoinfo.controls

    output.write("\t".join(columns)+"\n")
    for locus in data.loci:
        #column values for the row
        cols:Dict[str,Value] = defaultdict(lambda:None)

        locus_vars = locus.get_vars()
        lead = locus_vars.lead
        cols["locus_id"] = get_varid(lead.id)
        ld_partners = locus_vars.ld_partners if locus_vars.ld_partners is not None else []
        cs_vars = locus_vars.cs if locus_vars.cs is not None else []

        #phenotype cols
        cols["phenotype"] = phenotype
        cols["phenotype_abbreviation"] = pheno_abbr
        cols["Cases"] = cases
        cols["Controls"] = controls
        cols["chrom"] = lead.id.chrom
        #lead variant cols
        cols["pos"] = lead.id.pos
        cols["ref"] = lead.id.ref
        cols["alt"] = lead.id.alt
        cols["pval"] = lead.pval
        cols["lead_beta"] = lead.beta
        cols["lead_r2_threshold"] = ld_threshold(options.r2_threshold,options.ld_mode,lead.pval)
        #extra columns
        if ExtraColAnnotation.get_name() in data.annotations:
            extra_anno = data.annotations[ExtraColAnnotation.get_name()]
            if lead.id in extra_anno:
                for c in options.extra_columns:
                    cols[f"lead_{c}"] = extra_anno[lead.id][0][c]
        #previous release cols
        if PreviousReleaseAnnotation.get_name() in data.annotations:
            prev_anno = data.annotations[PreviousReleaseAnnotation.get_name()]
            if lead.id in prev_anno:
                for c in PreviousReleaseAnnotation.get_output_columns():
                    cols[f"lead_{c}"] = prev_anno[lead.id][0][c]
        #lead variant most severe gene, consequence
        if FGAnnotation.get_name() in data.annotations:
            fg_ann = data.annotations[FGAnnotation.get_name()]
            if lead.id in fg_ann:
                cols["lead_most_severe_gene"] = fg_ann[lead.id][0]["most_severe_gene"]
                cols["lead_most_severe_consequence"] = fg_ann[lead.id][0]["most_severe_consequence"]
                cols["gnomAD_functional_category"] = fg_ann[lead.id][0]["functional_category"]
                cols["rsids"] = fg_ann[lead.id][0]["rsids"]
        #lead enrichment
        if Gnomad4Annotation.get_name() in data.annotations:
            gnomad_ann = data.annotations[Gnomad4Annotation.get_name()]
            if lead.id in gnomad_ann:
                cols["lead_enrichment"] = gnomad_ann[lead.id][0]["GNOMAD_FI_enrichment_nfe"]
        if FunctionalAnnotation.get_name() in data.annotations:
            func_ann = data.annotations[FunctionalAnnotation.get_name()]
            func_add_cols = [
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
            if lead.id in func_ann:
                for c in func_add_cols:
                    cols[f"gnomAD_{c}"] = func_ann[lead.id][0][c]
        #cs data
        if CSAnnotation.get_name() in data.annotations:
            cs_ann = data.annotations[CSAnnotation.get_name()]
            cs_ann_cols = [
                "size",
                "bayes",
                "number",
                "region",
                "good_cs",
                "min_r2"
            ]
            cs_ann_col_d = {
                "size":"cs_size",
                "bayes":"cs_log_bayes_factor",
                "number":"cs_number",
                "region":"cs_region",
                "good_cs":"good_cs",
                "min_r2":"credible_set_min_r2_value"
            }
            cols["cs_id"] = get_varid(lead.id)+"_"+str(cs_ann[lead.id][0]["number"])

            if lead.id in cs_ann:
                for c in cs_ann_cols:
                    cols[cs_ann_col_d[c]] = cs_ann[lead.id][0][c]
        
        
        #locus start and end
        cols["start"] = min([a.id.pos for a in ld_partners+cs_vars+[lead]])
        cols["end"] = max([a.id.pos for a in ld_partners+cs_vars+[lead]])
        #create strict group
        if options.grouping_method == Grouping.LD:
            relaxed_set:set[Var] = set(cs_vars+ld_partners+[lead])
            strict_set:set[Var] = set([a for a in relaxed_set if 
                (a.r2_to_lead is not None) and (a.r2_to_lead >= options.strict_group) and (a.pval <= options.significance_threshold) ])
        elif options.grouping_method == Grouping.CS:
            strict_set = set(cs_vars+[lead])
            relaxed_set = set(ld_partners+cs_vars+[lead])
        else:
            relaxed_set = set(ld_partners+ cs_vars+[lead])
            strict_set = set([a for a in relaxed_set if a.pval <= options.significance_threshold])        
        
        #functional strict & relaxed
        if FGAnnotation.get_name() in data.annotations:
            fg_ann = data.annotations[FGAnnotation.get_name()]
            var_ids = [a for a in relaxed_set if a.id in fg_ann]
            func_vars_relaxed = set([a for a in var_ids if fg_ann[a.id][0]["functional_category"]!="NA"])
            func_vars_strict = [a for a in func_vars_relaxed if a in strict_set]
            func_data_strict = sorted([(
                get_varid(a.id),
                fg_ann[a.id][0]["functional_category"],
                fg_ann[a.id][0]["most_severe_gene"],
                a.r2_to_lead
            ) for a in func_vars_strict if a.r2_to_lead is not None],key=lambda x:x[3],reverse=True)
            #functional set of variants, relaxed
            func_data_relaxed = sorted([(
                get_varid(a.id),
                fg_ann[a.id][0]["functional_category"],
                fg_ann[a.id][0]["most_severe_gene"],
                a.r2_to_lead
            ) for a in func_vars_relaxed if a.r2_to_lead is not None],key=lambda x:x[3],reverse=True)
            cols["functional_variants_strict"] = ";".join([f"{a[0]}|{a[1]}|{a[2]}|{a[3]}" for a in func_data_strict])
            cols["functional_variants_relaxed"] = ";".join([f"{a[0]}|{a[1]}|{a[2]}|{a[3]}" for a in func_data_relaxed])
        
        #credset vars
        if options.grouping_method == Grouping.CS:
            try:
                cs_ann = data.annotations[CSAnnotation.get_name()]
            except:
                raise Exception("No cs annotation even though using cs grouping, this is a bug")
            try:
                cs_vars_extended = set(cs_vars+[lead])
                cs_data = {a.id:cs_ann[a.id][0]  for a in cs_vars_extended}
                _temp_lst:List[Tuple[Var,float]] = [(a,cast(float,cs_data[a.id]["prob"]) ) for a in cs_vars_extended ]
                sorted_cs_vars = sorted(_temp_lst,key=lambda x:x[1],reverse=True)
                cols["credible_set_variants"] = ";".join([f"{get_varid(a[0].id)}|{a[1]:.3g}|{a[0].r2_to_lead:.3g}" for a in sorted_cs_vars])
            except Exception as e:
                raise Exception("Error: bug",e)
        
        #all traits & strict traits
        #Use CatalogAnnotation
        if CatalogAnnotation.get_name() in data.annotations:
            cat_ann = data.annotations[CatalogAnnotation.get_name()]
            #strict
            trait_vars_strict = [a for a in strict_set if a.id in cat_ann]
            all_traits_strict:Dict[str,float] = {}
            for v in trait_vars_strict:
                for item in cat_ann[v.id]:
                    if cast(float,v.r2_to_lead) > all_traits_strict.get(cast(str,item["trait_name"]),0.0):
                        all_traits_strict[cast(str,item["trait_name"])] = cast(float,v.r2_to_lead)
            trait_data_strict = sorted(
                [
                    (trait_name,r2) for trait_name,r2 in all_traits_strict.items()
                ],
                key=lambda x:x[1],reverse=True
            )
            trait_vars_relaxed = [a for a in relaxed_set if a.id in cat_ann]
            all_traits_relaxed:dict[str,float] = {}
            for v in trait_vars_relaxed:
                for item in cat_ann[v.id]:
                    if cast(float,v.r2_to_lead) > all_traits_relaxed.get(cast(str,item["trait_name"]),0.0):
                        all_traits_relaxed[cast(str,item["trait_name"])] = cast(float,v.r2_to_lead)
            trait_data_relaxed = sorted(
                [
                    (trait_name,r2) for trait_name,r2 in all_traits_relaxed.items()
                ],
                key=lambda x:x[1],reverse=True
            )
            efo_traits = set(options.efo_traits)
            specific_traits_strict = [(a,b) for a,b in trait_data_strict if a in efo_traits]
            specific_traits_relaxed = [(a,b) for a,b in trait_data_relaxed if a in efo_traits]
            other_traits_strict = [(a,b) for a,b in trait_data_strict if a not in efo_traits]
            other_traits_relaxed = [(a,b) for a,b in trait_data_relaxed if a not in efo_traits]
            cols["specific_efo_trait_associations_strict"] = ";".join([f"{a}|{b:.3g}" for a,b in specific_traits_strict])
            cols["specific_efo_trait_associations_relaxed"] = ";".join([f"{a}|{b:.3g}" for a,b in specific_traits_relaxed])
            cols["found_associations_strict"] = ";".join([f"{a}|{b:.3g}" for a,b in other_traits_strict])
            cols["found_associations_relaxed"] = ";".join([f"{a}|{b:.3g}" for a,b in other_traits_relaxed])
        
        #n ld partners LD > 0.8
        r2_thresh_0_8 = 0.8
        r2_thresh_0_6 = 0.6
        #TODO: make sure that the cs vars are the same as they should be
        cs_set = set(cs_vars)
        opt_to_f = lambda x:x if x != None else 0.0
        n_ld_gt_0_8 = len([a for a in ld_partners if a not in cs_set and opt_to_f(a.r2_to_lead)>r2_thresh_0_8])
        n_ld_gt_0_6 = len([a for a in ld_partners if a not in cs_set and opt_to_f(a.r2_to_lead)>r2_thresh_0_6])
        cols["n_ld_partners_0_8"] = n_ld_gt_0_8
        cols["n_ld_partners_0_6"] = n_ld_gt_0_6
        try:
            output.write(
                "\t".join([format_value(cols[c]) for c in columns])+"\n"
            )
        except:
            raise
