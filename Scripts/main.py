#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd 
import numpy as np
import gws_fetch, compare, annotate

def main(args):
    ###########################
    ###Filter and Group SNPs###
    ###########################
    print("filter & group SNPs")
    args.out_fname="temp_tsv.out"
    args.annotate_fpath=args.out_fname
    args.compare_fname=args.out_fname
    gws_fetch.fetch_gws(args)
    
    ###########################
    ##########Finemap##########
    ###########################

    ###########################
    #######Annotate SNPs#######
    ###########################
    print("Annotate SNPs")
    annotate.annotate(args)
    ###########################
    ######Compare results######
    ###########################
    print("Compare results to previous findings")
    compare.compare(args)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="FINNGEN automatic hit reporting tool")
    
    #gws_fetch
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed tsv")
    parser.add_argument("--signifigance-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    #parser.add_argument("--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    parser.add_argument("--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("--alternate-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    #finemap
    
    #annotate
    #parser.add_argument("annotate_fpath",type=str,help="Filepath of the results to be annotated")
    parser.add_argument("--gnomad-genome-path",dest="gnomad_genome_path",type=str,help="Gnomad genome annotation file filepath")
    parser.add_argument("--gnomad-exome-path",dest="gnomad_exome_path",type=str,help="Gnomad exome annotation file filepath")
    parser.add_argument("--include-batch-freq",dest="batch_freq",action="store_true",help="Include batch frequencies from finngen annotations")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default=None,help="Finngen annotation file filepath")
    #parser.add_argument("--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    
    #compare results
    #parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,help="use 'file' or 'gwascatalog'")
    parser.add_argument("--summary-fpath",dest="summary_files",metavar="FILE",nargs="+",help="comparison summary filepaths")
    parser.add_argument("--endpoints",type=str,nargs="+",help="biological endpoint, as many as summaries")
    parser.add_argument("--build-38",dest="build_38",action="store_true",help="Whether is in GRCh38")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    #parser.add_argument("--ld-panel-path",dest="ld_panel_path",help="The path for the LD panel to determine what samples are in LD with each other")
    parser.add_argument("--raport-out",dest="raport_out",type=str,default="raport_output.csv",help="Comparison raport output path")
    parser.add_argument("--gwascatalog-pval",default=5e-8,help="P-value cutoff for GWASCatalog searches")
    parser.add_argument("--gwascatalog-width-kb",dest="gwascatalog_pad",type=int,default=25,help="gwascatalog range padding")
    parser.add_argument("--ldstore-threads",type=int,default=4,help="Number of threads to use with ldstore")
    parser.add_argument("--ld-treshold",type=float,default=0.4,help="ld treshold")
    args=parser.parse_args()
    main(args)