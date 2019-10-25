#!/usr/bin/env python3

import argparse,shlex,subprocess,os
from subprocess import Popen, PIPE
import pandas as pd 
import numpy as np
import tabix
from autoreporting_utils import *
#TODO: make a system for making sure we can calculate all necessary fields,
#e.g by checking that the columns exist


def calculate_enrichment(gnomad_df,fi_af_col,count_nfe_lst,number_nfe_lst):
    """Calculate enrichment for finns vs the other group, 
    which is defined in the column name lists
    In: gnomad dataframe, finnish allele frequency column,
    other group allele count and number columns
    Out: enrichment column"""
    nfe_counts=gnomad_df.loc[:,count_nfe_lst].sum(axis=1,skipna=True)
    nfe_numbers=gnomad_df.loc[:,number_nfe_lst].sum(axis=1,skipna=True)
    finn_freq=gnomad_df.loc[:,fi_af_col]
    enrichment=nfe_numbers*finn_freq/nfe_counts
    # clip enrichment
    enrichment=enrichment.clip(lower=0.0,upper=1e6)
    return enrichment

def create_rename_dict(list_of_names, prefix):
    d={}
    for value in list_of_names:
        d[value]="{}{}".format(prefix,value)
    return d

def annotate(args):
    """
    Annotates variants with allele frequencies, 
    enrichment numbers, and most severe gene/consequence data
    Annotations from gnomad exome data, gnomad genome data,
    finngen annotation file 
    """
    columns={"chrom":args.column_labels[0],"pos":args.column_labels[1],"ref":args.column_labels[2],"alt":args.column_labels[3],"pval":args.column_labels[4]}
    #columns that we want to take from gnomad and finngen annotations
    gnomad_gen_cols=["AF_fin","AF_nfe","AF_nfe_est","AF_nfe_nwe","AF_nfe_onf","AF_nfe_seu","FI_enrichment_nfe","FI_enrichment_nfe_est"]
    gnomad_exo_cols=["AF_nfe_bgr","AF_fin","AF_nfe","AF_nfe_est","AF_nfe_swe","AF_nfe_nwe","AF_nfe_onf",\
        "AF_nfe_seu","FI_enrichment_nfe","FI_enrichment_nfe_est","FI_enrichment_nfe_swe","FI_enrichment_nfe_est_swe"]
    finngen_cols=["most_severe_gene","most_severe_consequence"]

    #load main file
    df=pd.read_csv(args.annotate_fpath,sep="\t")
    if df.empty:
        df.to_csv(path_or_buf="{}".format(args.annotate_out),sep="\t",index=False)
        return
    #load gnomad_genomes
    if not os.path.exists("{}.tbi".format(args.gnomad_genome_path)):
        raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(args.gnomad_genome_path))
    tb_g=tabix.open(args.gnomad_genome_path)
    gnomad_genomes=load_tb_df(df,tb_g,args.gnomad_genome_path,columns=columns)

    #load gnomad_exomes
    if not os.path.exists("{}.tbi".format(args.gnomad_exome_path)):
        raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(args.gnomad_exome_path))
    tb_e=tabix.open(args.gnomad_exome_path)
    gnomad_exomes=load_tb_df(df,tb_e,args.gnomad_exome_path,columns=columns)

    #load finngen annotations
    if not os.path.exists("{}.tbi".format(args.finngen_path)):
        raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(args.finngen_path))
    tb_f=tabix.open(args.finngen_path)
    fg_df=load_tb_df(df,tb_f,args.finngen_path,chrom_prefix="chr",na_value="NA",columns=columns)
    fg_df=fg_df.drop_duplicates(subset=["#variant"])

    #load functional annotations. 
    if args.functional_path == "":
        func_df = pd.DataFrame(columns=["chrom","pos","ref","alt","consequence"])
    else:
        if not os.path.exists("{}.tbi".format(args.functional_path)):
            raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(args.functional_path))
        tb_func=tabix.open(args.functional_path)
        func_df = load_tb_df(df,tb_func,args.functional_path,chrom_prefix="chr",na_value="NA",columns=columns)
        func_cols=["chrom","pos","ref","alt","consequence"]
        func_df = func_df[func_cols]

    if not gnomad_genomes.empty:
        gnomad_genomes=gnomad_genomes.drop_duplicates(subset=["#CHROM","POS","REF","ALT"]).rename(columns={"#CHROM":columns["chrom"],"POS":columns["pos"],"REF":columns["ref"],"ALT":columns["alt"]})
        gnomad_genomes.loc[:,"#variant"]=create_variant_column(gnomad_genomes,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        #calculate enrichment for gnomad genomes, nfe, nfe without est
        gn_gen_nfe_counts=["AC_nfe","AC_nfe_est","AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
        gn_gen_nfe_nums=["AN_nfe","AN_nfe_est","AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
        gn_gen_nfe_est_counts=["AC_nfe","AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
        gn_gen_nfe_est_nums=["AN_nfe","AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
        #calculate enrichment
        
        gnomad_genomes.loc[:,"FI_enrichment_nfe"]=calculate_enrichment(gnomad_genomes,"AF_fin",gn_gen_nfe_counts,gn_gen_nfe_nums)
        gnomad_genomes.loc[:,"FI_enrichment_nfe_est"]=calculate_enrichment(gnomad_genomes,"AF_fin",gn_gen_nfe_est_counts,gn_gen_nfe_est_nums)
    else:
        for val in ["FI_enrichment_nfe","FI_enrichment_nfe_est","#variant"]:
            gnomad_genomes[val]=None
    
    if not gnomad_exomes.empty:
        gnomad_exomes=gnomad_exomes.drop_duplicates(subset=["#CHROM","POS","REF","ALT"]).rename(columns={"#CHROM":columns["chrom"],"POS":columns["pos"],"REF":columns["ref"],"ALT":columns["alt"]})
        gnomad_exomes.loc[:,"#variant"]=create_variant_column(gnomad_exomes,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        #calculate enrichment for gnomax exomes, nfe, nfe without est, nfe without swe, nfe without est, swe?
        gn_exo_nfe_counts=["AC_nfe","AC_nfe_bgr","AC_nfe_est","AC_nfe_onf","AC_nfe_seu","AC_nfe_swe"]
        gn_exo_nfe_nums=["AN_nfe","AN_nfe_bgr","AN_nfe_est","AN_nfe_onf","AN_nfe_seu","AN_nfe_swe"]
        gn_exo_nfe_est_counts=["AC_nfe","AC_nfe_bgr","AC_nfe_onf","AC_nfe_seu","AC_nfe_swe"]
        gn_exo_nfe_est_nums=["AN_nfe","AN_nfe_bgr","AN_nfe_onf","AN_nfe_seu","AN_nfe_swe"]
        gn_exo_nfe_swe_counts=["AC_nfe","AC_nfe_bgr","AC_nfe_est","AC_nfe_onf","AC_nfe_seu"]
        gn_exo_nfe_swe_nums=["AN_nfe","AN_nfe_bgr","AN_nfe_est","AN_nfe_onf","AN_nfe_seu"]
        gn_exo_nfe_est_swe_counts=["AC_nfe","AC_nfe_bgr","AC_nfe_onf","AC_nfe_seu"]
        gn_exo_nfe_est_swe_nums=["AN_nfe","AN_nfe_bgr","AN_nfe_onf","AN_nfe_seu"]
    
        gnomad_exomes.loc[:,"FI_enrichment_nfe"]=calculate_enrichment(gnomad_exomes,"AF_fin",gn_exo_nfe_counts,gn_exo_nfe_nums)
        gnomad_exomes.loc[:,"FI_enrichment_nfe_est"]=calculate_enrichment(gnomad_exomes,"AF_fin",gn_exo_nfe_est_counts,gn_exo_nfe_est_nums)
        gnomad_exomes.loc[:,"FI_enrichment_nfe_swe"]=calculate_enrichment(gnomad_exomes,"AF_fin",gn_exo_nfe_swe_counts,gn_exo_nfe_swe_nums)
        gnomad_exomes.loc[:,"FI_enrichment_nfe_est_swe"]=calculate_enrichment(gnomad_exomes,"AF_fin",gn_exo_nfe_est_swe_counts,gn_exo_nfe_est_swe_nums)
    else:
        for val in ["FI_enrichment_nfe","FI_enrichment_nfe_est","FI_enrichment_nfe_swe","FI_enrichment_nfe_est_swe","#variant"]:
            gnomad_exomes[val]=None
    
    if not func_df.empty:
        #rename columns
        func_df = func_df.drop_duplicates(subset=["chrom","pos","ref","alt"]).rename(columns={"chrom":columns["chrom"],"pos":columns["pos"],"ref":columns["ref"],"alt":columns["alt"],"consequence":"functional_category"})
        #remove values that are not in the coding categories
        functional_categories = ["pLoF","LC","start_lost","stop_lost","inframe_indel","missense_variant"]
        func_df["functional_category"] = func_df["functional_category"].apply(lambda x: x if x in functional_categories else np.nan)
        func_df[columns["chrom"]] = func_df[columns["chrom"]].apply(lambda x:x.strip("chr"))
        func_df = func_df.dropna(axis="index")
        func_df.loc[:,"#variant"] = create_variant_column(func_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        func_df = func_df[["#variant","functional_category"]]
    else:
        for val in ["#variant","functional_category"]:
            func_df[val]=None
        func_df = func_df[["#variant","functional_category"]]

    #rename gnomad_exomes and gnomad_genomes
    gnomad_genomes=gnomad_genomes.loc[:,["#variant"]+gnomad_gen_cols]
    gn_gen_rename_d=create_rename_dict(gnomad_gen_cols,"GENOME_")
    gnomad_genomes=gnomad_genomes.rename(columns=gn_gen_rename_d)

    gnomad_exomes=gnomad_exomes.loc[:,["#variant"]+gnomad_exo_cols]
    gn_exo_rename_d=create_rename_dict(gnomad_exo_cols,"EXOME_")
    gnomad_exomes=gnomad_exomes.rename(columns=gn_exo_rename_d)
    
    fg_df=fg_df.rename(columns={"gene":"most_severe_gene","most_severe":"most_severe_consequence"})
    fg_cols=fg_df.columns.values.tolist()
    info=list(filter(lambda s: "INFO" in s,fg_cols))
    imp=list(filter(lambda s: "IMP" in s,fg_cols))
    af=list(filter(lambda s: "AFW" in s,fg_cols))

    fg_batch_lst=info+imp+af
    fg_batch_rename=create_rename_dict(fg_batch_lst,"FG_")
    fg_df=fg_df.rename(columns=fg_batch_rename)
    fg_batch_lst=["FG_{}".format(x) for x in fg_batch_lst]
    if args.batch_freq:
        finngen_cols=finngen_cols+fg_batch_lst
    fg_df=fg_df.loc[:,["#variant"]+finngen_cols]

    #merge the wanted columns into df
    df=df.merge(gnomad_genomes,how="left",on="#variant")
    df=df.merge(gnomad_exomes,how="left",on="#variant")
    df=df.merge(func_df,how="left",on="#variant")
    df=df.merge(fg_df,how="left",on="#variant")

    df.to_csv(path_or_buf=args.annotate_out,sep="\t",index=False)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Annotate results using gnoMAD and additional annotations")
    parser.add_argument("annotate_fpath",type=str,help="Filepath of the results to be annotated")
    parser.add_argument("--gnomad-genome-path",dest="gnomad_genome_path",type=str,help="Gnomad genome annotation file filepath")
    parser.add_argument("--gnomad-exome-path",dest="gnomad_exome_path",type=str,help="Gnomad exome annotation file filepath")
    parser.add_argument("--include-batch-freq",dest="batch_freq",action="store_true",help="Include batch frequencies from finngen annotations")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default="",help="Finngen annotation file filepath")
    parser.add_argument("--functional-path",dest="functional_path",type=str,default="",help="File path to functional annotations file")
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix. Default value is the base name (no path and no file extensions) of input file. ")
    parser.add_argument("--annotate-out",dest="annotate_out",type=str,default="annotate_out.csv",help="Output filename, default is out.csv")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval"],help="Names for data file columns. Default is '#chrom pos ref alt pval'.")
    args=parser.parse_args()
    if args.prefix!="":
        args.prefix=args.prefix+"."
    args.annotate_out = "{}{}".format(args.prefix,args.annotate_out)
    if (args.gnomad_exome_path == None) or (args.gnomad_genome_path == None) or (args.finngen_path==None):
        print("Annotation files missing, aborting...")
    else:    
        annotate(args)