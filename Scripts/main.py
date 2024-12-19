#!/usr/bin/env python3
from grouping_report import generate_top_report, generate_variant_report,TopReportOptions,VariantReportOptions
import argparse
from data_access.linkage import PlinkLD, OnlineLD, TabixLD
from grouping import form_groups
from grouping_model import Grouping, LDMode, PhenoData, SummstatColumns,GroupingOptions
from group_annotation import CSAnnotation, ExtraColAnnotation, FunctionalAnnotation, PreviousReleaseAnnotation, PreviousReleaseOptions, FGAnnotation, annotate, CatalogAnnotation, Gnomad4Annotation
from phenoinfo import get_phenotype_data, PhenoInfoOptions
from load_tabix import tb_resource_manager,TabixOptions
from data_access.csfactory import csfactory
from typing import Optional, List
import os


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
        elif args.ld_api_choice == "tabix":
            ld_api = TabixLD(args.ld_panel_path)
        else:
            raise ValueError("Wrong argument for --ld-api:{}".format(args.ld_api_choice))
    args.sig_treshold_2 = max(args.sig_treshold, args.sig_treshold_2)
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

    ### Annotation resources
    prevrel_annotation = None
    if args.previous_release_path:
        prevrel_opts = PreviousReleaseOptions(args.previous_release_path,
            column_names.c,
            column_names.p,
            column_names.r,
            column_names.a,
            column_names.pval,
            column_names.beta,
            args.previous_release_additional_columns)
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
    #gnomad_genome_annotation = None
    #if args.gnomad_genome_path:
    #    gnomad_genome_annotation = GnomadGenomeAnnotation(args.gnomad_genome_path)
    # gnomad genome annotations
    #gnomad_exome_annotation = None
    #if args.gnomad_exome_path:
    #    gnomad_exome_annotation = GnomadExomeAnnotation(args.gnomad_exome_path)
    # gnomad 4 annotation
    gnomad_four_annotation = None
    if args.gnomad_path:
        gnomad_four_annotation = Gnomad4Annotation(args.gnomad_path)
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
            #gnomad_genome_annotation,
            #gnomad_exome_annotation,
            gnomad_four_annotation,
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
        
        ### create phenotype
        phenodata = PhenoData(phenotype_info,loci)
        ### annotate
        phenodata = annotate(phenodata,annotation_resources)
        ### create report

        report_fname = create_fname(args.report_out,args.prefix)
        top_fname = create_fname(args.top_report_out,args.prefix)
        with open(report_fname,"w") as report_file, open(top_fname,"w") as top_file:
            generate_variant_report(phenodata,report_file,variant_report_options,annotation_resources)
            generate_top_report(phenodata,top_file,top_report_options,annotation_resources)
    

def create_fname(path:str,prefix:str)->str:
    head,tail = os.path.split(path)
    return os.path.join(head,f"{prefix}{tail}")

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="FINNGEN automatic hit reporting tool")
    
    #gws_fetch
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed summary statistic file to be processed")
    parser.add_argument("--sign-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix")
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
    parser.add_argument("--ld-api",dest="ld_api_choice",type=str,default="plink",help="LD interface to use. Valid options are 'plink', 'online' and 'tabix'.")
    parser.add_argument("--pheno-name",dest="pheno_name",type=str,default="",help="Phenotype name")
    parser.add_argument("--pheno-info-file",dest="pheno_info_file",type=str,default="",help="Phenotype information file path")
    parser.add_argument("--extra-cols",dest="extra_cols",nargs="*",default=[],help="extra columns in the summary statistic you want to add to the results")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL","BETA"),nargs=6,default=["#chrom","pos","ref","alt","pval","beta"],help="Names for data file columns. Default is '#chrom pos ref alt pval beta'.")
    
    #annotate
    parser.add_argument("--gnomad-path",dest="gnomad_path",type=str,help="Gnomad 4 annotation filepath")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default="",help="Finngen annotation file filepath")
    parser.add_argument("--functional-path",dest="functional_path",type=str,default="",help="File path to functional annotations file")
    parser.add_argument("--previous-release-path",dest="previous_release_path",type=str,default="",help="File path to previous release summary statistic file")
    parser.add_argument("--previous-release-extra-cols",dest="previous_release_additional_columns",nargs = "*",default=[],help="Additional columns apart from pval and beta to take from previous release summary statistic file")
    parser.add_argument("--use-gwascatalog",action="store_true",help="Add flag to use GWAS Catalog for comparison.")
    parser.add_argument("--custom-dataresource",type=str,default="",help="Custom dataresource path.")
    parser.add_argument("--report-out",dest="report_out",type=str,default="report_out.tsv",help="Comparison report output path")
    parser.add_argument("--gwascatalog-pval",default=5e-8,type=float,help="P-value cutoff for GWASCatalog searches")
    parser.add_argument("--gwascatalog-width-kb",dest="gwascatalog_pad",type=int,default=0,help="gwascatalog range padding")
    parser.add_argument("--gwascatalog-threads",dest="gwascatalog_threads",type=int,default=4,help="Number of concurrent queries to GWAScatalog API. Default 4. Increase if the gwascatalog api takes too long.")
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
