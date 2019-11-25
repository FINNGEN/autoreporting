#!/usr/bin/env python3

import argparse,shlex,subprocess
import pandas as pd 
import numpy as np
import gws_fetch, compare, annotate,autoreporting_utils
from linkage import PlinkLD, OnlineLD

def main(args):
    print("input file: {}".format(args.gws_fpath))
    args.fetch_out = "{}{}".format(args.prefix,args.fetch_out)
    args.annotate_out = "{}{}".format(args.prefix,args.annotate_out)
    args.report_out = "{}{}".format(args.prefix,args.report_out)
    args.top_report_out = "{}{}".format(args.prefix,args.top_report_out)
    args.ld_report_out = "{}{}".format(args.prefix,args.ld_report_out)
    args.sig_treshold_2=max(args.sig_treshold_2,args.sig_treshold)
    ld_api=None
    if args.ld_api_choice == "plink":
        ld_api = PlinkLD(args.ld_panel_path,args.plink_mem)
    elif args.ld_api_choice == "online":
        ld_api = OnlineLD(url="http://api.finngen.fi/api/ld")
    else:
        raise ValueError("Wrong argument for --ld-api:{}".format(args.ld_api_choice))
    ###########################
    ###Filter and Group SNPs###
    ###########################
    print("filter & group SNPs")
    args.annotate_fpath=args.fetch_out
    args.compare_fname=args.annotate_out
    fetch_df = gws_fetch.fetch_gws(gws_fpath=args.gws_fpath, sig_tresh_1=args.sig_treshold, prefix=args.prefix, group=args.grouping, grouping_method=args.grouping_method, locus_width=args.loc_width,
        sig_tresh_2=args.sig_treshold_2, ld_r2=args.ld_r2, overlap=args.overlap, column_labels=args.column_labels,
        ignore_region=args.ignore_region, cred_set_file=args.cred_set_file,ld_api=ld_api)
    
    #write fetch_df as a file, so that other parts of the script work
    fetch_df.to_csv(path_or_buf=args.fetch_out,sep="\t",index=False)
    ###########################
    ##########Finemap##########
    ###########################

    ###########################
    #######Annotate SNPs#######
    ###########################
    if (args.gnomad_exome_path == None) or (args.gnomad_genome_path == None) or (args.finngen_path==None):
        print("Annotation files missing, skipping gnomad & finngen annotation...")
        #args.compare_fname=args.annotate_fpath
        annotate_df = fetch_df
    else:
        print("Annotate SNPs")
        #annotate_df = annotate.annotate(fetch_df,args)
        annotate_df = annotate.annotate(df=fetch_df,gnomad_genome_path=args.gnomad_genome_path, gnomad_exome_path=args.gnomad_exome_path, batch_freq=args.batch_freq, finngen_path=args.finngen_path, fg_ann_version = args.fg_ann_version,
            functional_path=args.functional_path, prefix=args.prefix, column_labels=args.column_labels)
    annotate_df.to_csv(path_or_buf=args.annotate_out,sep="\t",index=False)
    ###########################
    ######Compare results######
    ###########################
    print("Compare results to previous findings")
    [report_df,ld_out_df] = compare.compare(annotate_df,compare_style=args.compare_style, summary_fpath=args.summary_fpath, endpoints=args.endpoints,ld_check=args.ld_check,
                                    plink_mem=args.plink_mem, ld_panel_path=args.ld_panel_path, prefix=args.prefix,
                                    gwascatalog_pval=args.gwascatalog_pval, gwascatalog_pad=args.gwascatalog_pad, gwascatalog_threads=args.gwascatalog_threads,
                                    ldstore_threads=args.ldstore_threads, ld_treshold=args.ld_treshold, cache_gwas=args.cache_gwas, column_labels=args.column_labels,
                                    localdb_path=args.localdb_path, database_choice=args.database_choice)
    if type(report_df) != type(None):
        report_df.to_csv(args.report_out,sep="\t",index=False)
        #create top report
        #top level df 
        columns=autoreporting_utils.columns_from_arguments(args.column_labels)
        top_df=compare.create_top_level_report(report_df,args.efo_traits,columns)
        top_df.to_csv(args.top_report_out,sep="\t",index=False)
    if type(ld_out_df) != type(None):
        ld_out_df.to_csv(args.ld_report_out,sep="\t")

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="FINNGEN automatic hit reporting tool")
    
    #gws_fetch
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed summary statistic file to be processed")
    parser.add_argument("--sign-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix")
    parser.add_argument("--fetch-out",dest="fetch_out",type=str,default="fetch_out.csv",help="GWS output filename, default is fetch_out.csv")
    parser.add_argument("--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("--alt-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    parser.add_argument("--overlap",dest="overlap",action="store_true",help="Are groups allowed to overlap")
    parser.add_argument("--ignore-region",dest="ignore_region",type=str,default="",help="Ignore the given region, e.g. HLA region, from analysis. Give in CHROM:BPSTART-BPEND format.")
    parser.add_argument("--credible-set-file",dest="cred_set_file",type=str,default="",help="bgzipped SuSiE credible set file.")
    parser.add_argument("--ld-api",dest="ld_api_choice",type=str,default="plink",help="LD interface to use. Valid options are 'plink' and 'online'.")
    #finemap
    
    #annotate
    parser.add_argument("--gnomad-genome-path",dest="gnomad_genome_path",type=str,help="Gnomad genome annotation file filepath")
    parser.add_argument("--gnomad-exome-path",dest="gnomad_exome_path",type=str,help="Gnomad exome annotation file filepath")
    parser.add_argument("--include-batch-freq",dest="batch_freq",action="store_true",help="Include batch frequencies from finngen annotations")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default="",help="Finngen annotation file filepath")
    parser.add_argument("--functional-path",dest="functional_path",type=str,default="",help="File path to functional annotations file")
    parser.add_argument("--annotate-out",dest="annotate_out",type=str,default="annotate_out.csv",help="Annotation output filename, default is annotate_out.csv")
    parser.add_argument("--finngen-annotation-version",dest="fg_ann_version",type=str,default="r3",help="Finngen annotation release version: 3 or under or 4 or higher? Allowed values: 'r3' and 'r4'. Default 'r3' ")
    
    #compare results
    parser.add_argument("--compare-style",type=str,default="gwascatalog",help="use 'file', 'gwascatalog' or 'both'")
    parser.add_argument("--summary-fpath",dest="summary_fpath",type=str,help="Summary listing file path.")
    parser.add_argument("--endpoint-fpath",dest="endpoints",type=str,help="Endpoint listing file path.")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--report-out",dest="report_out",type=str,default="report_out.csv",help="Comparison report output path")
    parser.add_argument("--ld-report-out",dest="ld_report_out",type=str,default="ld_report_out.csv",help="LD check report output path")
    parser.add_argument("--gwascatalog-pval",default=5e-8,help="P-value cutoff for GWASCatalog searches")
    parser.add_argument("--gwascatalog-width-kb",dest="gwascatalog_pad",type=int,default=25,help="gwascatalog range padding")
    parser.add_argument("--gwascatalog-threads",dest="gwascatalog_threads",type=int,default=4,help="Number of concurrent queries to GWAScatalog API. Default 4. Increase if the gwascatalog api takes too long.")
    parser.add_argument("--ldstore-threads",type=int,default=4,help="Number of threads to use with ldstore. Default 4")
    parser.add_argument("--ld-treshold",type=float,default=0.9,help="ld treshold for including ld associations in ld report")
    parser.add_argument("--cache-gwas",action="store_true",help="save gwascatalog results into gwas_out_mapping.csv and load them from there if it exists. Use only for testing.")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval"],help="Names for data file columns. Default is '#chrom pos ref alt pval'.")
    parser.add_argument("--top-report-out",dest="top_report_out",type=str,default="top_report.csv",help="Top level report filename.")
    parser.add_argument("--efo-codes",dest="efo_traits",type=str,nargs="+",default=[],help="Specific EFO codes to look for in the top level report")
    parser.add_argument("--local-gwascatalog",dest='localdb_path',type=str,help="Path to local GWAS Catalog DB.")
    parser.add_argument("--db",dest="database_choice",type=str,default="gwas",help="Database to use for comparison. use 'local','gwas' or 'summary_stats'.")
    args=parser.parse_args()
    if args.prefix!="":
        args.prefix=args.prefix+"."
    main(args)
