#!/usr/bin/env python3

import argparse,shlex,subprocess, glob
from subprocess import Popen, PIPE
import sys,os,io
import pandas as pd, numpy as np
#import tabix
from autoreporting_utils import *

def parse_region(region):
    chrom=region.split(":")[0]
    start=region.split(":")[1].split("-")[0]
    end=region.split(":")[1].split("-")[1]
    return {"chrom":str(chrom),"start":int(start),"end":int(end)}


def solve_groups(all_variants,group_data):
    """
    Returns a dataframe containing the grouped variants, and only those.
    In: All variants considered for grouping (basically the variants with p-value < args.sig_tresh_2), parsed group data
    Out: Dataframe with same columns as all_variants, containing the rows in group_data
    """
    retval=pd.DataFrame(columns=all_variants.columns)
    for t in group_data.itertuples():
        group=[t.SNP]
        if t.TOTAL>0:
            sp2_split=t.SP2.split(",")
            strip_= lambda x,y: x[:-len(y)] if x.endswith(y) else x
            sp2_split=[strip_(x.strip(),"(1)" ) for x in sp2_split]
            group=group+sp2_split
        group_df=all_variants[ all_variants["#variant"].isin(group) ].copy()
        group_df.loc[:,"locus_id"]=t.SNP
        retval=pd.concat([retval,group_df],axis="index",sort=True)
    return retval

def simple_grouping(df_p1,df_p2,r,overlap,columns):
    """
    Simple grouping function
    Groups the variants based on p-values
    Group 1: lead variants
    Group 2: lead variants + variants that can not be lead variants. group 2 contains group 1.
    In: group 1 (df_p1, df), group 2 (df_p2, df), group width (r, int), do we overlap (overlap, bool), columns (columns, dict{key:str})
    Out: grouped dataframe
    """
    #simple grouping
    df=df_p1.copy()
    group_df=df_p2.copy()
    new_df=pd.DataFrame(columns=df.columns)
    while not df.empty:
        ms_snp=df.loc[df[columns["pval"] ].idxmin(),:]#get most sig SNP
        rowidx=(group_df[ columns["pos"] ]<=(ms_snp[columns["pos"]]+r) )&(group_df[ columns["pos"] ]>= (ms_snp[columns["pos"]]-r) )&(group_df[ columns["chrom"] ]==ms_snp[columns["chrom"]])#get group indexes
        tmp=group_df.loc[rowidx,:].copy()
        tmp.loc[:,"locus_id"]=ms_snp["#variant"]#necessary
        tmp.loc[:,"pos_rmin"]=tmp[columns["pos"]].min()#get group pos_rmin and pos_rmax from location of the SNPs instead of just a 
        tmp.loc[:,"pos_rmax"]=tmp[columns["pos"]].max()
        new_df=pd.concat([new_df,tmp],ignore_index=True,axis=0,join='inner')
        #convergence: remove the indexes from result_dframe. The rowidx =\= dropidx, as they come from different dataframes (df, group_df)
        dropidx=(df[ columns["pos"] ]<=(ms_snp[columns["pos"]]+r) )&(df[ columns["pos"] ]>= (ms_snp[columns["pos"]]-r) )&(df[ columns["chrom"] ]==ms_snp[columns["chrom"]])
        df=df.loc[~dropidx,:]
        #if not overlap, we also delete these from the group_df, so they can't be included again in other groups.
        if not overlap:
            t_dropidx=(group_df[ columns["pos"] ]<=(ms_snp[columns["pos"]]+r) )&(group_df[ columns["pos"] ]>=(ms_snp[columns["pos"]]-r) )&(group_df[ columns["chrom"] ]==ms_snp[columns["chrom"]])
            group_df=group_df.loc[~t_dropidx,:]
    return new_df

def ld_grouping(df_p1,df_p2, sig_treshold , sig_treshold_2, locus_width, ld_treshold,ld_panel_path,plink_memory,overlap, prefix, columns):
    """
    LD Clumping function
    Groups the variants based on PLINK's ld-clumping
    In: variant group 1, variant group 2, p-value threshold 1, p-value threshold 2, group width, ld threshold, prefix,  columns
    Out: grouped dataframe
    """
    #1:create PLINK variant list
    temp_variants="{}clump_variants.csv".format(prefix)
    df_p2.loc[:,["#variant",columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"] ]].to_csv(path_or_buf=temp_variants,index=False,sep="\t")
    plink_fname="{}plink_clump".format(prefix)
    #set up overlap flag for PLINK
    allow_overlap=""
    if overlap==True:
        allow_overlap="--clump-allow-overlap"
    plink_command="plink --allow-extra-chr --bfile {} --clump {} --clump-field {} --clump-snp-field '{}'  --clump-r2 {}"\
        " --clump-kb {} --clump-p1 {} --clump-p2 {} --out {} --memory {} {}".format(
        ld_panel_path,
        temp_variants,
        columns["pval"],
        "#variant",
        ld_treshold,
        locus_width,
        sig_treshold,
        sig_treshold_2,
        plink_fname,
        plink_memory,
        allow_overlap)
    #run PLINK
    pr = subprocess.Popen(shlex.split(plink_command), stdout=PIPE,stderr=subprocess.STDOUT,encoding='ASCII' )
    pr.wait()
    #get plink log
    plink_log=pr.stdout.readlines()
    if pr.returncode!=0:
        print("PLINK FAILURE. Error code {}".format(pr.returncode)  )
        [print(l) for l in plink_log]
        raise ValueError("Plink clumping returned code {}".format(pr.returncode))
    #Check if PLINK returned something or not
    no_sig_res_string="Warning: No significant --clump results.  Skipping."
    #if there is data
    if os.path.exists("{}.clumped".format(plink_fname)):
        group_data=pd.read_csv("{}.clumped".format(plink_fname),sep="\s+")
        group_data=group_data.loc[:,["SNP","TOTAL","SP2"]]
        new_df=solve_groups(df_p2.copy(),group_data)
        for var in new_df["locus_id"].unique():
            new_df.loc[new_df["locus_id"]==var,"pos_rmin"]=new_df.loc[new_df["locus_id"]==var,"pos"].min()#r["min"]
            new_df.loc[new_df["locus_id"]==var,"pos_rmax"]=new_df.loc[new_df["locus_id"]==var,"pos"].max()
        p1_group_leads = df_p1["#variant"].isin(new_df["#variant"])
        p1_singletons = ~p1_group_leads
        new_df=pd.concat([new_df,df_p1.loc[p1_singletons,:]],axis="index",sort=True).sort_values(by=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],"#variant"])
        new_df.loc[:,"pos_rmin"]=new_df.loc[:,"pos_rmin"].astype(np.int32)
        new_df.loc[:,"pos_rmax"]=new_df.loc[:,"pos_rmax"].astype(np.int32)
    #if there is not data
    else:
        if any([no_sig_res_string in string for string in plink_log]):#no significant results
            print("No significant results with PLINK clumping. All groups are singletons.")
            new_df=df_p1.copy()
        else:
            print("Plink .clumped file not found. Check the logs for information:")
            [print(l) for l in plink_log]
            raise FileNotFoundError("Plink .clumped file not found.")
    #cleanup plink files
    plink_files=glob.glob("{}.*".format(plink_fname))
    subprocess.call(["rm",temp_variants]+plink_files,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    return new_df

def fetch_gws(args):
    #column names
    columns={"chrom":args.column_labels[0],"pos":args.column_labels[1],"ref":args.column_labels[2],"alt":args.column_labels[3],"pval":args.column_labels[4]}

    sig_tresh=max(args.sig_treshold,args.sig_treshold_2)
    c_size=100000
    r=args.loc_width*1000#range for location width, originally in kb
    dtype={columns["chrom"]:str,
                columns["pos"]:np.int32,
                columns["ref"]:str,
                columns["alt"]:str,
                columns["pval"]:np.float64}
    temp_df=pd.DataFrame()
    for dframe in pd.read_csv(args.gws_fpath,compression="gzip",sep="\t",dtype=dtype,engine="c",chunksize=c_size):
        #filter gws snps
        temp_df=pd.concat([temp_df,dframe.loc[dframe[columns["pval"]]<=sig_tresh,:]],axis="index",ignore_index=True)
    #remove ignored region if there is one
    if args.ignore_region:
        ignore_region=parse_region(args.ignore_region)
        ign_idx=( ( temp_df[columns["chrom"]]==ignore_region["chrom"] ) & ( temp_df[columns["pos"]]<=ignore_region["end"] )&( temp_df[columns["pos"]]>=ignore_region["start"] ) )
        temp_df=temp_df.loc[~ign_idx,:]
    
    if temp_df.empty:
        print("The input file {} contains no gws-significant hits with signifigance treshold of {}. Aborting.".format(args.gws_fpath,args.sig_treshold))
        return 1
    temp_df=temp_df.reset_index(drop=True)
    temp_df.loc[:,"#variant"]=create_variant_column(temp_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
    temp_df.loc[:,"locus_id"]=temp_df.loc[:,"#variant"]
    temp_df.loc[:,"pos_rmax"]=temp_df.loc[:,columns["pos"]]
    temp_df.loc[:,"pos_rmin"]=temp_df.loc[:,columns["pos"]]
    df_p1=temp_df.loc[temp_df[columns["pval"]] <= args.sig_treshold,: ].copy()
    df_p2=temp_df.loc[temp_df[columns["pval"]] <= args.sig_treshold_2,: ].copy()
    if args.grouping and not df_p1.empty:
        if args.grouping_method=="ld":
            new_df=ld_grouping(df_p1,df_p2,args.sig_treshold,args.sig_treshold_2,args.loc_width,args.ld_r2,args.ld_panel_path,args.plink_mem,args.overlap,args.prefix,columns)
        else:
            new_df=simple_grouping(df_p1=df_p1,df_p2=df_p2,r=r,overlap=args.overlap,columns=columns)
        new_df=new_df.sort_values(["locus_id","#variant"])
        new_df.to_csv(path_or_buf=args.fetch_out,sep="\t",index=False)
    else:
        #take only gws hits, no groups. Therefore, use df_p1
        df_p1.sort_values(["locus_id","#variant"]).to_csv(path_or_buf=args.fetch_out,sep="\t",index=False)
    return 0
    
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Fetch and group genome-wide significant variants from summary statistic")
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed summary statistic")
    parser.add_argument("--sign-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix. Default value is the base name (no path and no file extensions) of input file. ")
    parser.add_argument("--fetch-out",dest="fetch_out",type=str,default="fetch_out.csv",help="GWS output filename, default is fetch_out.csv")
    parser.add_argument("--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("--alt-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    parser.add_argument("--overlap",dest="overlap",action="store_true",help="Are groups allowed to overlap")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval"],help="Names for data file columns. Default is '#chrom pos ref alt pval'.")
    parser.add_argument("--ignore-region",dest="ignore_region",type=str,default="",help="Ignore the given region, e.g. HLA region, from analysis. Give in CHROM:BPSTART-BPEND format.")
    args=parser.parse_args()
    if args.prefix!="":
        args.prefix=args.prefix+"."
    args.fetch_out = "{}{}".format(args.prefix,args.fetch_out)
    fetch_gws(args)