#!/usr/bin/env python3

import argparse,shlex,subprocess,os
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
import tabix
from typing import Dict, Optional
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
    """
    Create a dictionary for renaming columns from different annotation files
    For eaxmple, given column names ["AF_1","AF_2","AF_3"] and a prefix "GNOMAD_",
    return rename_dict={"AF_1":"GNOMAD_AF1","AF_2":"GNOMAD_AF2","AF_3":"GNOMAD_AF3"}
    This can then be used as df.rename(columns=rename_dict)
    In: list of column names, prefix
    Out: dictionary with entries "oldname":"prefixoldname"
    """
    d={}
    for value in list_of_names:
        d[value]="{}{}".format(prefix,value)
    return d

def previous_release_annotate(fpath: Optional[str], df: pd.DataFrame, columns: Dict[str,str]) -> pd.DataFrame:
    """Create the previous release annotation
    Args:
        fpath (str): filepath to the previous release summary statistic
        df (pd.DataFrame): Variant dataframe
        columns (Dict[str,str]): column dictionary
    Returns:
        (pd.DataFrame): Dataframe with columns [#variant, beta_previous_release, pval_previous_release]
    """
    previous_cols = [columns["chrom"], columns["pos"], columns["ref"], columns["alt"], "beta_previous_release", "pval_previous_release"]
    if fpath:
        if not os.path.exists("{}.tbi".format(fpath)):
            raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(fpath))
         
        previous_df = load_tb_df(df,fpath, chrom_prefix="", na_value="", columns=columns)
        previous_df = previous_df.rename(columns={"beta":"beta_previous_release","pval":"pval_previous_release"})
        previous_df = previous_df[previous_cols]
    else:
        return pd.DataFrame(columns = ["#variant","beta_previous_release", "pval_previous_release"])
    
    if not previous_df.empty:
        previous_df = previous_df.drop_duplicates(subset=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"]])
        previous_df = df_replace_value(previous_df,columns["chrom"],"X","23")
        previous_df["#variant"] = create_variant_column(previous_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
    else:
        previous_df["#variant"] = None
    previous_df = previous_df[["#variant", "beta_previous_release", "pval_previous_release"]]
    return previous_df


def annotate(df: pd.DataFrame, gnomad_genome_path: str, gnomad_exome_path: str, batch_freq: bool, finngen_path: str, functional_path: str, previous_release_path: str ,prefix: str, columns: Dict[str, str]) -> pd.DataFrame :
    """
    Annotates variants with allele frequencies, enrichment numbers, and most severe gene/consequence data
    Annotations from gnomad exome data, gnomad genome data, finngen annotation file, functional annotation file.
    Args:
        df (pd.DataFrame): Variant dataframe
        gnomad_genome_path (str): gnomad genome annotation file path
        gnomad_exome_path (str): gnomad exome annotation file path
        batch_freq (bool): flag for whether to include batch-specific frequencies from finngen annotation
        finngen_path (str): finngen annotation file path
        functional_path (str): functional annotation file path
        previous_release_path (str): filepath for the previous release
        prefix (str): prefix for analysis files
        columns (Dict[str, str]): column dictionary
    Returns:
        (pd.DataFrame): Annotated dataframe
    Out: Annotated dataframe
    """
    #columns that we want to take from gnomad and finngen annotations
    gnomad_gen_cols=["AF_fin","AF_nfe","AF_nfe_est","AF_nfe_nwe","AF_nfe_onf","AF_nfe_seu","FI_enrichment_nfe","FI_enrichment_nfe_est"]
    gnomad_exo_cols=["AF_nfe_bgr","AF_fin","AF_nfe","AF_nfe_est","AF_nfe_swe","AF_nfe_nwe","AF_nfe_onf",\
        "AF_nfe_seu","FI_enrichment_nfe","FI_enrichment_nfe_est","FI_enrichment_nfe_swe","FI_enrichment_nfe_est_swe"]
    finngen_cols=["most_severe_gene","most_severe_consequence"]

    if df.empty:
        return df

    #create chr 23->X calling df
    #needs chromosome 23 as X
    call_df_x = df.copy()
    call_df_x[columns["chrom"]]=call_df_x[columns["chrom"]].astype(str)
    call_df_x = df_replace_value(call_df_x,columns["chrom"],"23","X") #TODO: IF/WHEN GNOMAD RESOURCES USE CHR 23, THIS NEEDS TO BE REMOVED
    #load gnomad_genomes
    if not os.path.exists("{}.tbi".format(gnomad_genome_path)):
        raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(gnomad_genome_path))
    gnomad_genomes=load_tb_df(call_df_x,gnomad_genome_path,columns=columns)

    #load gnomad_exomes
    if not os.path.exists("{}.tbi".format(gnomad_exome_path)):
        raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(gnomad_exome_path))
    gnomad_exomes=load_tb_df(call_df_x,gnomad_exome_path,columns=columns)
    #replace X with 23
    gnomad_genomes = df_replace_value(gnomad_genomes,"#CHROM","X","23")
    gnomad_exomes = df_replace_value(gnomad_exomes,"#CHROM","X","23")
     
    #load finngen annotations
    if not os.path.exists("{}.tbi".format(finngen_path)):
        raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(finngen_path))

    fg_df=load_tb_df(df,finngen_path,chrom_prefix="",na_value="NA",columns=columns)
    fg_df=fg_df.drop(labels="#variant",axis="columns")
    fg_df["#variant"]=create_variant_column(fg_df,chrom="chr",pos="pos",ref="ref",alt="alt")
    fg_df = fg_df.drop_duplicates(subset=["#variant"]).rename(columns={"chrom":columns["chrom"],"pos":columns["pos"],"ref":columns["ref"],"alt":columns["alt"]})
    #sanity check: if number of variants is >0 and FG annotations are smaller, emit a warning message.
    if df.shape[0]>0 and all(fg_df["INFO"].isna()):
        print("Warning: FG annotation does not have any hits but the input data has. Check that you are using a recent version of the finngen annotation file (R3_v1 or above)")

    #load functional annotations. 
    if functional_path == "":
        func_df = pd.DataFrame(columns=["chrom","pos","ref","alt","consequence"])
    else:
        if not os.path.exists("{}.tbi".format(functional_path)):
            raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(functional_path))
        func_df = load_tb_df(call_df_x,functional_path,chrom_prefix="chr",na_value="NA",columns=columns)
        func_df["chrom"] = func_df["chrom"].apply(lambda x:x.strip("chr"))
        func_df=df_replace_value(func_df,"chrom","X","23")
        func_cols=["chrom","pos","ref","alt","consequence"]
        func_df = func_df[func_cols]


    #load previous release annotations
    previous_df = previous_release_annotate(previous_release_path,call_df_x,columns)

    if not gnomad_genomes.empty:
        gnomad_genomes=gnomad_genomes.drop_duplicates(subset=["#CHROM","POS","REF","ALT"]).rename(columns={"#CHROM":columns["chrom"],"POS":columns["pos"],"REF":columns["ref"],"ALT":columns["alt"]})
        gnomad_genomes["#variant"]=create_variant_column(gnomad_genomes,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        #calculate enrichment for gnomad genomes, nfe, nfe without est
        gn_gen_nfe_counts=["AC_nfe_est","AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
        gn_gen_nfe_nums=["AN_nfe_est","AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
        gn_gen_nfe_est_counts=["AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
        gn_gen_nfe_est_nums=["AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
        #calculate enrichment
        
        gnomad_genomes.loc[:,"FI_enrichment_nfe"]=calculate_enrichment(gnomad_genomes,"AF_fin",gn_gen_nfe_counts,gn_gen_nfe_nums)
        gnomad_genomes.loc[:,"FI_enrichment_nfe_est"]=calculate_enrichment(gnomad_genomes,"AF_fin",gn_gen_nfe_est_counts,gn_gen_nfe_est_nums)
    else:
        for val in ["FI_enrichment_nfe","FI_enrichment_nfe_est","#variant"]:
            gnomad_genomes[val]=None
    
    if not gnomad_exomes.empty:
        gnomad_exomes=gnomad_exomes.drop_duplicates(subset=["#CHROM","POS","REF","ALT"]).rename(columns={"#CHROM":columns["chrom"],"POS":columns["pos"],"REF":columns["ref"],"ALT":columns["alt"]})
        gnomad_exomes["#variant"]=create_variant_column(gnomad_exomes,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        #calculate enrichment for gnomax exomes, nfe, nfe without est, nfe without swe, nfe without est, swe?
        gn_exo_nfe_counts=["AC_nfe_bgr","AC_nfe_est","AC_nfe_onf","AC_nfe_seu","AC_nfe_swe"]
        gn_exo_nfe_nums=["AN_nfe_bgr","AN_nfe_est","AN_nfe_onf","AN_nfe_seu","AN_nfe_swe"]
        gn_exo_nfe_est_counts=["AC_nfe_bgr","AC_nfe_onf","AC_nfe_seu","AC_nfe_swe"]
        gn_exo_nfe_est_nums=["AN_nfe_bgr","AN_nfe_onf","AN_nfe_seu","AN_nfe_swe"]
        gn_exo_nfe_swe_counts=["AC_nfe_bgr","AC_nfe_est","AC_nfe_onf","AC_nfe_seu"]
        gn_exo_nfe_swe_nums=["AN_nfe_bgr","AN_nfe_est","AN_nfe_onf","AN_nfe_seu"]
        gn_exo_nfe_est_swe_counts=["AC_nfe_bgr","AC_nfe_onf","AC_nfe_seu"]
        gn_exo_nfe_est_swe_nums=["AN_nfe_bgr","AN_nfe_onf","AN_nfe_seu"]
    
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
        functional_categories = ["pLoF","LC","start_lost","stop_lost","stop_gained","inframe_indel","missense_variant"]
        func_df["functional_category"] = func_df["functional_category"].apply(lambda x: x if x in functional_categories else np.nan)
        func_df = func_df.dropna(axis="index")
        func_df["#variant"] = create_variant_column(func_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
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
    if batch_freq:
        finngen_cols=finngen_cols+fg_batch_lst
    fg_df=fg_df.loc[:,["#variant"]+finngen_cols]

    #merge the wanted columns into df
    df=df.merge(gnomad_genomes,how="left",on="#variant")
    df=df.merge(gnomad_exomes,how="left",on="#variant")
    df=df.merge(func_df,how="left",on="#variant")
    df=df.merge(fg_df,how="left",on="#variant")
    df=df.merge(previous_df,how="left",on="#variant")

    return df

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Annotate results using gnoMAD and additional annotations")
    parser.add_argument("annotate_fpath",type=str,help="Filepath of the results to be annotated")
    parser.add_argument("--gnomad-genome-path",dest="gnomad_genome_path",type=str,help="Gnomad genome annotation file filepath")
    parser.add_argument("--gnomad-exome-path",dest="gnomad_exome_path",type=str,help="Gnomad exome annotation file filepath")
    parser.add_argument("--include-batch-freq",dest="batch_freq",action="store_true",help="Include batch frequencies from finngen annotations")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default="",help="Finngen annotation file filepath")
    parser.add_argument("--functional-path",dest="functional_path",type=str,default="",help="File path to functional annotations file")
    parser.add_argument("--previous-release-path",dest="previous_release_path",type=str,help="File path to previous release summary statistic file")
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix. Default value is the base name (no path and no file extensions) of input file. ")
    parser.add_argument("--annotate-out",dest="annotate_out",type=str,default="annotate_out.tsv",help="Output filename, default is out.tsv")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL","BETA","AF","AF_CASE","AF_CONTROL"),nargs=9,default=["#chrom","pos","ref","alt","pval","beta","maf","maf_cases","maf_controls"],help="Names for data file columns. Default is '#chrom pos ref alt pval beta maf maf_cases maf_controls'.")
    args=parser.parse_args()
    columns=columns_from_arguments(args.column_labels)
    if args.prefix!="":
        args.prefix=args.prefix+"."
    args.annotate_out = "{}{}".format(args.prefix,args.annotate_out)
    if (args.gnomad_exome_path == None) or (args.gnomad_genome_path == None) or (args.finngen_path==None):
        print("Annotation files missing, aborting...")
    else:    
        input_df = pd.read_csv(args.annotate_fpath,sep="\t")
        df = annotate(df=input_df,gnomad_genome_path=args.gnomad_genome_path, gnomad_exome_path=args.gnomad_exome_path, batch_freq=args.batch_freq, finngen_path=args.finngen_path,
        functional_path=args.functional_path, previous_release_path=args.previous_release_path, prefix=args.prefix, columns=columns)
        df.fillna("NA").replace("","NA").to_csv(path_or_buf=args.annotate_out,sep="\t",index=False,float_format="%.3g")
