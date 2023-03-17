#!/usr/bin/env python3
from grouping_report import generate_top_report, generate_variant_report,TopReportOptions,VariantReportOptions
import argparse,shlex,subprocess
import pandas as pd 
import numpy as np
from data_access import datafactory, csfactory
from data_access.linkage import PlinkLD, OnlineLD
from grouping import form_groups
from grouping_model import Grouping, LDMode, PhenoInfo, PhenoData, SummstatColumns,GroupingOptions,CSInfo,Var,Locus,PeakLocus,CSLocus
from group_annotation import CSAnnotation, ExtraColAnnotation, FunctionalAnnotation, PreviousReleaseAnnotation, PreviousReleaseOptions, FGAnnotation, TabixOptions, annotate, GnomadExomeAnnotation, GnomadGenomeAnnotation, CatalogAnnotation
from phenoinfo import get_phenotype_data, PhenoInfoOptions
from load_tabix import TabixResource, tb_resource_manager
from data_access.csfactory import csfactory
from data_access.db import Variant, CSAccess
import os
from typing import Optional, List

def groupingModeFactory(method:str,group:bool):
    if not group:
        return Grouping.NONE
    if method=="simple":
        return Grouping.RANGE
    elif method=="ld":
        return Grouping.LD
    elif method=="cred":
        return Grouping.CS
    else:
        raise Exception(f"Unknown grouping method {method} passed to --grouping-method")

def ldModeFactory(ld_chisq):
    if ld_chisq != None:
        return LDMode.DYNAMIC
    return LDMode.CONSTANT

def main(args):
    print("input file: {}".format(args.gws_fpath))
    ### Construct options & resources
    # grouping mode
    gr_mode = groupingModeFactory(args.grouping_method,args.grouping)
    #summstat columns
    column_names = SummstatColumns(args.column_labels[0],args.column_labels[1],args.column_labels[2],args.column_labels[3],args.column_labels[4],args.column_labels[5])
    #range
    loc_range = args.loc_width*1000
    # ld mode
    ld_mode = ldModeFactory(args.dynamic_r2_chisq)
    if ld_mode == LDMode.DYNAMIC:
        r2_threshold = args.dynamic_r2_chisq
    else:
        r2_threshold = args.ld_r2
    # ld resource
    ld_api=None
    if args.grouping_method != "simple":
        if args.ld_api_choice == "plink":
            ld_api = PlinkLD(args.ld_panel_path,args.plink_mem)
        elif args.ld_api_choice == "online":
            ld_api = OnlineLD(url="http://api.finngen.fi/api/ld")
        else:
            raise ValueError("Wrong argument for --ld-api:{}".format(args.ld_api_choice))
    
    gr_opts = GroupingOptions(args.gws_fpath,
        gr_mode,
        column_names,
        loc_range,
        ld_mode,
        r2_threshold,
        args.sig_treshold,
        args.sig_treshold_2,
        args.overlap
    )

    ### Annotation options
    # extra columns come from summstat_resource
    #previous release
    prevrel_opts = PreviousReleaseOptions(args.previous_release_path,
        column_names.c,
        column_names.p,
        column_names.r,
        column_names.a,
        column_names.pval,
        column_names.beta)
    prevrel_annotation = PreviousReleaseAnnotation(prevrel_opts)
    extra_cols_annotation = ExtraColAnnotation(TabixOptions(
        args.gws_fpath,
        column_names.c,
        column_names.p,
        column_names.r,
        column_names.a,
    ),args.extra_cols) if args.extra_cols else None
    # functional annotation
    functional_annotation = None
    if args.functional_path:
        functional_annotation = FunctionalAnnotation(args.functional_path)
    # finngen annotation
    fg_annotation = None
    if args.finngen_path:
        fg_annotation = FGAnnotation(args.finngen_path)
    # gnomad exome annotation
    gnomad_genome_annotation = None
    if args.gnomad_genome_path:
        gnomad_genome_annotation = GnomadGenomeAnnotation(args.gnomad_genome_path)
    # gnomad genome annotations
    gnomad_exome_annotation = None
    if args.gnomad_exome_path:
        gnomad_exome_annotation = GnomadExomeAnnotation(args.gnomad_exome_path)
    #gwas catalog annotation
    catalog_annotation = CatalogAnnotation(args.use_gwascatalog,
                                                    args.custom_dataresource,
                                                    args.database_choice,
                                                    args.localdb_path,
                                                    args.gwascatalog_pad,
                                                    args.gwascatalog_pval,
                                                    args.gwascatalog_threads,
                                                    args.allele_db_file)

    ### Report options
    top_report_options = TopReportOptions(args.strict_group_r2,
        args.sig_treshold,
        gr_mode,
        args.efo_traits,
        args.extra_cols,
        ld_mode,
        r2_threshold)
    variant_report_options = VariantReportOptions(column_names,args.extra_cols,False)

    ### create resources for grouping
    # cs resource
    if args.cred_set_file:
        cs_access = csfactory(args.cred_set_file)
    else:
        cs_access=None
    cs_annotation = None if cs_access == None else CSAnnotation(cs_access)

    #join all annotation resources
    annotation_resources = [a for a in (
            prevrel_annotation,
            extra_cols_annotation,
            functional_annotation,
            fg_annotation,
            gnomad_genome_annotation,
            gnomad_exome_annotation,
            catalog_annotation,
            cs_annotation
        ) if a != None
    ]

    #create phenotype info
    pheno_options = None
    if args.pheno_info_file != "":
        pheno_options = PhenoInfoOptions(args.pheno_info_file,args.pheno_name)
    phenotype_info = get_phenotype_data(pheno_options)

    #summary statistic resource
    with tb_resource_manager(args.gws_fpath,
        gr_opts.column_names.c,
        gr_opts.column_names.p,
        gr_opts.column_names.r,
        gr_opts.column_names.a) as summstat_resource:
        
        ### group
        loci = form_groups(summstat_resource,gr_opts,cs_access,ld_api)
        if os.path.exists("./quickstart_df.tsv"):
            fetch_df = pd.read_csv("./quickstart_df.tsv",sep="\t")
        list_of_loci  = pandas_df_to_loci(fetch_df,args.grouping,args.grouping_method,cs_access)
        
        #compare
        #TESTING
        #TODO: delet this
        locd = {a.lead.id:a for a in loci}
        locd2 = {a.lead.id:a for a in list_of_loci}
        keys = [a.id for a in sorted([b.lead for b in loci],key=lambda x:x.pval)]
        for key in keys:
            l1  = locd[key]
            if key not in locd2:
                print(f"Key {key} not in l2: Instead the following were there: {[a for a in locd2.keys()]}")
                continue
            l2 = locd2[key]
            if l1.lead.id != l2.lead.id:
                print(f"leads different: {l1.lead} {l2.lead}")
            l1cs = [a.id for a in l1.cs]
            l2cs = [a.id for a in l2.cs]
            if set(l1cs) != set(l2cs):
                print(f"cs different: {l1cs} {l2cs}")
            l1ld = [a.id for a in l1.ld_partners]
            l2ld = [a.id for a in l2.ld_partners]
            if set(l1ld) != set(l2ld):
                print(f"ld partners different for {key}")
                l1s = set(l1ld)
                l2s = set(l2ld)
                l1_ = len(l1ld)
                l2_ = len(l2ld)
                li = len(l1s.intersection(l2s))
                diff1 = len(l1s.difference(l2s))
                diff2 = len(l2s.difference(l1s))
                print(f"set intersection {li}, diff1 {diff1}, diff2 {diff2}, 1size {l1_} 2size {l2_}")
                for v in l1.ld_partners:
                    print(v)

        raise NotImplementedError
        ### create phenotype
        phenodata = PhenoData(phenotype_info,loci)
        ### annotate
        phenodata = annotate(phenodata,annotation_resources)
        ### create report
        with open(args.report_out,"w") as report_file, open(args.top_report_out,"w") as top_file:
            generate_variant_report(phenodata,report_file,variant_report_options)
            generate_top_report(phenodata,top_file,top_report_options)
    

def pandas_df_to_loci(pandas_df,group:bool,grouping_method:str,csaccess:Optional[CSAccess])->List[Locus]:
    output = []
    #get locus ids
    csd = {}
    if csaccess:
        for c in csaccess.get_cs():
            csd[c.lead] = CSInfo(
                c.region,
                c.lead,
                c.number,
                c.bayes,
                c.min_r2,
                c.size,
                c.good_cs
            )
    if not group:
        #amount of locus ids should be the same as the amount of variants
        #therefore, we can just take the locus ids, make them into variants, take pval and beta and that's it
        for row in pandas_df[["locus_id","pval","beta"]].itertuples():
            identifier = row.locus_id
            s = identifier.replace("chr","").split("_")
            output.append(
                PeakLocus(
                    Var(Variant(s[0],int(s[1]),s[2],s[3]),row.pval,row.beta,None),
                    None,
                    Grouping.NONE
                )
            )
    else:
        locus_ids = list(pandas_df["locus_id"].unique())
        for l in locus_ids:
            if grouping_method == "simple":
                #lead var is lead var, and the range vars is every other var
                lead_var_row = pandas_df[(pandas_df["locus_id"]==l)&(pandas_df["#variant"]==l)].iloc[0]
                s = l.replace("chr","").split("_")
                lead_var = Var(Variant(s[0],int(s[1]),s[2],s[3] ),float(lead_var_row["pval"]),float(lead_var_row["beta"]),None)
                range_var_df = pandas_df[(pandas_df["locus_id"]==l)&(pandas_df["#variant"] != l)][["#variant","locus_id","pval","beta"]]
                range_vals = [(a._1.replace("chr","").split("_"),a.pval,a.beta) for a in range_var_df.itertuples()]
                range_vars = [Var(Variant(a[0][0],int(a[0][1]),a[0][2],a[0][3]),float(a[1]),float(a[2]),None ) for a in range_vals]
                grouping = Grouping.RANGE
                output.append(PeakLocus(lead_var, range_vars,grouping))
            elif grouping_method == "ld":
                lead_var_row = pandas_df[(pandas_df["locus_id"]==l)&(pandas_df["#variant"]==l)].iloc[0]
                s = l.replace("chr","").split("_")
                lead_var = Var(Variant(s[0],int(s[1]),s[2],s[3] ),float(lead_var_row["pval"]),float(lead_var_row["beta"]),1.0)
                ld_var_df = pandas_df[(pandas_df["locus_id"]==l) & (pandas_df["#variant"] != l)][["#variant","locus_id","pval","beta","r2_to_lead"]]
                ld_vals = [(a._1.replace("chr","").split("_"),a.pval,a.beta,a.r2_to_lead) for a in ld_var_df.itertuples()]
                ld_partners = [Var(Variant(a[0][0],int(a[0][1]),a[0][2],a[0][3]),float(a[1]),float(a[2]),float(a[3]) ) for a in ld_vals]
                grouping = Grouping.LD
                output.append(PeakLocus(lead_var, ld_partners,grouping))
            elif grouping_method == "cred":
                lead_var_row = pandas_df[(pandas_df["locus_id"]==l)&(pandas_df["#variant"]==l)].iloc[0]
                s = l.replace("chr","").split("_")
                lead_var = Var(Variant(s[0],int(s[1]),s[2],s[3] ),float(lead_var_row["pval"]),float(lead_var_row["beta"]),1.0)
                cs_var_df = pandas_df[(pandas_df["locus_id"]==l)& (pandas_df["cs_id"] == lead_var_row["cs_id"])][["#variant","locus_id","pval","beta","r2_to_lead"]]
                cs_vals = [(a._1.replace("chr","").split("_"),a.pval,a.beta,a.r2_to_lead) for a in cs_var_df.itertuples()]
                cs = [Var(Variant(a[0][0],int(a[0][1]),a[0][2],a[0][3]),float(a[1]),float(a[2]),float(a[3]) ) for a in cs_vals]
                ld_var_df = pandas_df[(pandas_df["locus_id"]==l)& (pandas_df["cs_id"] != lead_var_row["cs_id"])& (pandas_df["#variant"] != l)][["#variant","locus_id","pval","beta","r2_to_lead"]]
                ld_vals = [(a._1.replace("chr","").split("_"),a.pval,a.beta,a.r2_to_lead) for a in ld_var_df.itertuples()]
                ld_partners = [Var(Variant(a[0][0],int(a[0][1]),a[0][2],a[0][3]),float(a[1]),float(a[2]),float(a[3]) ) for a in ld_vals]
                cs_info = csd[lead_var.id]
                output.append(CSLocus(lead_var,cs,cs_info,ld_partners))
    return output

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="FINNGEN automatic hit reporting tool")
    
    #gws_fetch
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed summary statistic file to be processed")
    parser.add_argument("--sign-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix")
    parser.add_argument("--fetch-out",dest="fetch_out",type=str,default="fetch_out.tsv",help="GWS output filename, default is fetch_out.tsv")
    parser.add_argument("--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("--alt-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    
    #r2 static bound or dynamic per peak
    r2_group = parser.add_mutually_exclusive_group()
    r2_group.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    r2_group.add_argument("--dynamic-r2-chisq",type=float,nargs="?",const=5.0,default=None,help="If flag is passed, r2 threshold is set per peak so that leadvar_chisq*r2=value (default 5).")
    
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    parser.add_argument("--overlap",dest="overlap",action="store_true",help="Are groups allowed to overlap")
    parser.add_argument("--ignore-region",dest="ignore_region",type=str,default="",help="Ignore the given region, e.g. HLA region, from analysis. Give in CHROM:BPSTART-BPEND format.")
    parser.add_argument("--credible-set-file",dest="cred_set_file",type=str,default="",help="bgzipped SuSiE credible set file.")
    parser.add_argument("--ld-api",dest="ld_api_choice",type=str,default="plink",help="LD interface to use. Valid options are 'plink' and 'online'.")
    parser.add_argument("--pheno-name",dest="pheno_name",type=str,default="",help="Phenotype name")
    parser.add_argument("--pheno-info-file",dest="pheno_info_file",type=str,default="",help="Phenotype information file path")
    parser.add_argument("--extra-cols",dest="extra_cols",nargs="*",default=[],help="extra columns in the summary statistic you want to add to the results")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL","BETA"),nargs=6,default=["#chrom","pos","ref","alt","pval","beta"],help="Names for data file columns. Default is '#chrom pos ref alt pval beta'.")
    
    #annotate
    parser.add_argument("--gnomad-genome-path",dest="gnomad_genome_path",type=str,help="Gnomad genome annotation file filepath")
    parser.add_argument("--gnomad-exome-path",dest="gnomad_exome_path",type=str,help="Gnomad exome annotation file filepath")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default="",help="Finngen annotation file filepath")
    parser.add_argument("--functional-path",dest="functional_path",type=str,default="",help="File path to functional annotations file")
    parser.add_argument("--previous-release-path",dest="previous_release_path",type=str,default="",help="File path to previous release summary statistic file")
    #parser.add_argument("--annotate-out",dest="annotate_out",type=str,default="annotate_out.tsv",help="Annotation output filename, default is annotate_out.tsv")
    
    #compare results
    parser.add_argument("--use-gwascatalog",action="store_true",help="Add flag to use GWAS Catalog for comparison.")
    parser.add_argument("--custom-dataresource",type=str,default="",help="Custom dataresource path.")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--report-out",dest="report_out",type=str,default="report_out.tsv",help="Comparison report output path")
    parser.add_argument("--ld-report-out",dest="ld_report_out",type=str,default="ld_report_out.rsv",help="LD check report output path")
    parser.add_argument("--gwascatalog-pval",default=5e-8,type=float,help="P-value cutoff for GWASCatalog searches")
    parser.add_argument("--gwascatalog-width-kb",dest="gwascatalog_pad",type=int,default=0,help="gwascatalog range padding")
    parser.add_argument("--gwascatalog-threads",dest="gwascatalog_threads",type=int,default=4,help="Number of concurrent queries to GWAScatalog API. Default 4. Increase if the gwascatalog api takes too long.")
    parser.add_argument("--ldstore-threads",type=int,default=4,help="Number of threads to use with ldstore. Default 4")
    parser.add_argument("--ld-threshold",type=float,default=0.9,help="ld threshold for including ld associations in ld report")
    parser.add_argument("--cache-gwas",action="store_true",help="save gwascatalog results into gwas_out_mapping.tsv and load them from there if it exists. Use only for testing.")
    parser.add_argument("--local-gwascatalog",dest='localdb_path',type=str,help="Path to local GWAS Catalog file.")
    parser.add_argument("--db",dest="database_choice",type=str,choices=['local','gwas','summary_stats'],default="gwas",help="Database to use for comparison. use 'local','gwas' or 'summary_stats'.")
    parser.add_argument("--gwascatalog-allele-file",dest="allele_db_file",default="",help="GWAS Catalog alleles taken from here. Use dbSNP variation VCF file (current gwcat build hg38p13, b153).")

    #top report creation
    parser.add_argument("--top-report-out",dest="top_report_out",type=str,default="top_report.tsv",help="Top level report filename.")
    parser.add_argument("--strict-group-r2",dest="strict_group_r2",type=float,default=0.5,help="R^2 threshold for including variants in strict groups in top report")
    parser.add_argument("--efo-codes",dest="efo_traits",type=str,nargs="+",default=[],help="Specific EFO codes to look for in the top level report")
    
    args=parser.parse_args()
    if args.prefix!="":
        args.prefix=args.prefix+"."
    main(args)
