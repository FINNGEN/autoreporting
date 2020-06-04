#!/usr/bin/env python3

import argparse,shlex,subprocess, glob
from subprocess import Popen, PIPE
import sys,os,io
import pandas as pd, numpy as np
import tabix
from typing import Dict, List
from autoreporting_utils import *
from data_access.linkage import PlinkLD, OnlineLD
from data_access.db import LDAccess

def parse_region(region):
    chrom=region.split(":")[0]
    start=region.split(":")[1].split("-")[0]
    end=region.split(":")[1].split("-")[1]
    return {"chrom":str(chrom),"start":int(start),"end":int(end)}

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

def load_credsets(fname,columns):
    """
    Load SuSiE credible sets from one bgzipped file.
    In: filename of the SuSie credible set file
    Out: A Dataframe containing all credible sets for this phenotype.
    """
    input_data = pd.read_csv(fname,sep="\t",compression="gzip")
    if input_data.empty:
        return pd.DataFrame(columns = columns.values()+["cs_prob","cs_id"])
    input_data = input_data[input_data["cs"]!=-1]#filter to credible sets
    input_data["credsetid"]=input_data[["region","cs"]].apply(lambda x: "".join([str(y) for y in x]),axis=1)
    data=input_data.rename(columns={"chromosome":columns["chrom"], "position": columns["pos"], "allele1": columns["ref"], "allele2": columns["alt"], "prob": "cs_prob"}).copy()
    data[columns["chrom"]] = data[columns["chrom"]].str.strip("chr")
    data["cs_id"]=np.NaN
    for i in data["credsetid"].unique():
        cs=data[data["credsetid"]==i]
        rsid=cs.loc[cs["cs_prob"].idxmax(),"rsid"]
        idx=cs.loc[cs["cs_prob"].idxmax(),"cs"]
        data.loc[data["credsetid"]==i,"cs_id"] = "{}_{}".format(rsid,idx)
    cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"], "cs_prob", "cs_id" ]
    data=data.loc[:,cols]
    return data

def ld_grouping(df_p1,df_p2, sig_treshold_2,locus_width,ld_treshold, overlap,prefix, ld_api, columns):
    """
    Create groups based on the LD between variants.
    In: df filtered with p1, df filtered with p2, 
    """
    all_variants=df_p2.copy()
    group_leads = df_p1.copy()
    ld_ranges = group_leads[ [columns["chrom"], columns["pos"], columns["ref"], columns["alt"], "#variant"] ].rename(columns={ columns["chrom"]:"chr", columns["pos"]:"pos", columns["ref"]:"ref", columns["alt"]:"alt" })
    ld_data = ld_api.get_ranges(ld_ranges,locus_width*1000,ld_treshold)
    ld_df = pd.merge(all_variants[ ["#variant",columns["pval"] ] ],ld_data,how="inner",left_on="#variant",right_on="variant_2" )
    ld_df=ld_df.drop(columns=["chrom_2","pos_2","variant_2"])
    ld_df = ld_df[ld_df[columns["pval"]] <= sig_treshold_2 ]
    out_df = pd.DataFrame(columns=list(df_p2.columns)+["r2_to_lead"])
    #Grouping: greedily group those variants that have not been yet grouped into the most significant variants.
    #Range has been taken care of in the LD fetching, as well as r^2 threshold, and p-value was taken care of in filtering the ld_df.
    while not group_leads.empty:
        lead_variant = group_leads.loc[group_leads[columns["pval"]].idxmin(),"#variant" ]
        ld_data = ld_df.loc[ld_df["variant_1"] == lead_variant,["#variant","r2"] ].copy().rename(columns={"r2":"r2_to_lead"})
        group = all_variants[all_variants["#variant"].isin(ld_data["#variant"]) ].copy()
        group= pd.merge(group,ld_data,on="#variant",how="left")
        grouplead=all_variants[all_variants["#variant"]==lead_variant].copy()
        grouplead["r2_to_lead"]=1.0#the group lead is, of course, in perfect LD with the group lead
        group=pd.concat([group,grouplead],ignore_index=True,axis=0,join='inner').drop_duplicates(subset=["#variant"])
        group["locus_id"]=lead_variant
        group["pos_rmin"]=group[columns["pos"]].min()
        group["pos_rmax"]=group[columns["pos"]].max()
        out_df=pd.concat([out_df,group],ignore_index=True,axis=0,join='inner')
        #remove all of the variants with p<sig_tresh from lead_variants, since those in groups can not become leads
        group_leads=group_leads[ ~group_leads["#variant"].isin( group["#variant"].unique() ) ]
        #overlap
        if not overlap:
            all_variants=all_variants[~all_variants["#variant"].isin( group["#variant"].unique() )]
    return out_df


def credible_set_grouping(data: pd.DataFrame, ld_threshold: float, locus_range: int, overlap: bool, ld_api: LDAccess, columns: Dict[str, str]) -> pd.DataFrame:

    """Group variants using credible sets
    Create groups using credible set most probable variants as the lead variants, and rest of the data as the additional variants
    Args:
        data (pd.DataFrame): Input data
        ld_threshold (float): LD threshold for including a variant in the group
        locus_range (int): Range around the locus on which LD is calculated, in kb
        overlap (bool): Whether groups can overlap or not
        ld_api (LDAccess): ld api object
        columns (Dict[str,str]): column dictionary
    Returns:
        (pd.DataFrame): Grouped variants in a pandas dataframe
    """
    df=data.copy()
    lead_vars=[]
    #determine group leads. Group leads are 'the variants with largest cs_prob in that credible set'.
    for credible_set in df.loc[~df["cs_id"].isna(),"cs_id"].unique():
        group_lead =  df.loc[ df.loc[df["cs_id"]==credible_set,"cs_prob"].idxmax(),:]
        lead_vars.append(group_lead["#variant"])
    if len(lead_vars) == 0:
        return pd.DataFrame(columns=df.columns)
    lead_df = df.loc[df["#variant"].isin(lead_vars)].copy()
    ld_ranges=lead_df[ [columns["chrom"], columns["pos"], columns["ref"], columns["alt"], "#variant"] ].rename(columns={ columns["chrom"]:"chr", columns["pos"]:"pos", columns["ref"]:"ref", columns["alt"]:"alt" })
    all_lead_ld_data=ld_api.get_ranges(ld_ranges,locus_range*1000)
    #join
    ld_df = pd.merge(df[["#variant",columns["chrom"],columns["pos"],columns["pval"]]],all_lead_ld_data, how="inner",left_on="#variant",right_on="variant_2") #does include all of the lead variants as well
    ld_df=ld_df.drop(columns=["chrom_2","pos_2","variant_2"])
    #filter by p-value
    out_df = pd.DataFrame(columns=list(data.columns)+["r2_to_lead"])
    #create df with only lead variants
    leads = df[df["#variant"].isin(lead_vars)].loc[:,["#variant",columns["pval"]]].copy()
    while not leads.empty:
        #get lead variant and ld data
        lead_variant = leads.loc[leads[columns["pval"]].idxmin(),"#variant"] #choose the lead variant with smallest p-value
        ld_data = ld_df.loc[ld_df["variant_1"] == lead_variant,["#variant","r2"] ].copy().rename(columns={"r2":"r2_to_lead"})
        
        #create credible set group and ld partner group
        credible_id = data.loc[data["#variant"]==lead_variant,"cs_id"].values[0]
        credible_set= data.loc[data["cs_id"] == credible_id,:].copy()
        cred_group = credible_set.merge(ld_data,on="#variant",how="left")#contains all variants of this CS, even though LD might be smaller than ld threshold
        ld_data_filtered = ld_data[ld_data["r2_to_lead"]>=ld_threshold].copy() #filter ld
        group = df[df["#variant"].isin(ld_data_filtered["#variant"] )].copy()
        group = pd.merge(group,ld_data_filtered,on="#variant",how="left")
        #concat the two, remove duplicate entries. Entries with cs_id are preferred over entries without cs_id.
        #Though it shouldn't be possible for there to be variants that are both in the credible set and out of it. 
        group=pd.concat([group,cred_group],ignore_index=True,sort=False).sort_values(by=["cs_id","#variant","r2_to_lead"]).drop_duplicates(subset=["#variant"],keep="first")
        group["locus_id"]=lead_variant
        group["pos_rmin"]=group[columns["pos"]].min()
        group["pos_rmax"]=group[columns["pos"]].max()
        
        #add the group to output
        out_df=pd.concat([out_df,group],ignore_index=True,axis=0,join='inner')
        
        #convergence: remove lead_variant, remove group from df if overlap is not true
        leads=leads[~ (leads["#variant"] == lead_variant) ]
        if not overlap:
            df=df[~df["#variant"].isin(ld_data_filtered["#variant"])]
            df=df[~(df["cs_id"]==credible_id)]
    return out_df

def extract_cols(df: pd.DataFrame, cols: List[str])-> pd.DataFrame:
    """Extract columns from a dataframe
    Args:
    Returns:
        (pd.DataFrame): The dataframe with only those columns
    """
    try:
        df=df[ cols ]
    except KeyError:
        raise KeyError("DataFrame did not contain all of the required columns. Missing columns:{} Supplied columns:{}  ".format([a for a in cols if a not in df.columns],
            cols))
    except:
        raise
    return df

def get_gws_variants(fname, sign_treshold=5e-8,dtype=None,columns={},extra_cols=[],compression="gzip"):
    """
    Get genome-wide significant variants from a summary statistic file.
    In: filename, significance threshold, dtype,columns,compression
    Out: dataframe containing the significant variants. No additional columns will be added.
    """
    chunksize=100000
    if not dtype:
        dtype={columns["chrom"]:str,
                columns["pos"]:np.int32,
                columns["ref"]:str,
                columns["alt"]:str,
                columns["pval"]:np.float64}
    retval=pd.DataFrame()
    for df in pd.read_csv(fname,compression=compression,sep="\t",dtype=dtype,engine="c",chunksize=chunksize):
        retval=pd.concat( [retval,df.loc[df[columns["pval"] ] <=sign_treshold,: ] ], axis="index", ignore_index=True,sort=False )
    extracted_cols=list(columns.values())+extra_cols
    retval=extract_cols(retval,extracted_cols)
    return retval

def merge_credset(gws_df,cs_df,fname,columns):
    """
    Merge credible set to the genome-wide significant variants. 
    In case variants in the credible set are not included in the gws variants, 
    the rows corresponding to them are fetched using tabix.
    In: Dataframe containing gws variants, dataframe containing credible sets, filename for summary statistic.
    Out: Dataframe containing the gws variants + any credible set variants that are not gws. Columns 'cs_id','cs_prob' added to the dataframe. 
    """
    join_cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"]]
    # fetch rows using tabix
    cred_row_df = load_tb_df(cs_df,fname,columns=columns)
    cols = list(gws_df.columns)
    cred_row_df = cred_row_df[ cols ].drop_duplicates(keep="first")
    cred_row_df=cred_row_df.astype(dtype={columns["chrom"]:str,columns["pos"]:np.int64,columns["ref"]:str,columns["alt"]:str})
    cs_df=cs_df.astype(dtype={columns["chrom"]:str,columns["pos"]:np.int64,columns["ref"]:str,columns["alt"]:str})
    # ensure only the credible sets were included
    cred_row_df = pd.merge(cred_row_df,cs_df[join_cols],how="right",on=join_cols)
    df = pd.concat( [gws_df,cred_row_df], axis="index", ignore_index=True, sort=False)\
        .astype(dtype={columns["chrom"]:str,columns["pos"]:np.int64,columns["ref"]:str,columns["alt"]:str})\
        .drop_duplicates(subset=list( join_cols ) )
    # merge the credible set
    merged = pd.merge(df,cs_df,how="left",on=join_cols)
    return merged

def fetch_gws(gws_fpath: str, sig_tresh_1: float, prefix: str, group: bool, grouping_method: str, locus_width: int, sig_tresh_2: float,
                ld_r2: float, overlap: bool,columns: Dict[str,str], ignore_region: str, cred_set_file: str, ld_api: LDAccess, extra_cols: List[str]):
    """Filter and group variants.
    Args:
        gws_fpath (str): summary statistic filename
        sig_tresh_1 (float): primary significance threshold (used in simple and ld grouping)
        prefix (str): prefix for all of the files
        group (bool): whether to group the variants or not
        grouping_method (str): grouping method
        locus_width (int): grouping width in kb
        sig_tresh_2 (float): alternate significance threshold
        ld_r2 (float): LD correlation threshold for including variants in groups (for cred and ld grouping only)
        overlap (bool): Whether groups are allowed to overlap, i.e. same variant can be in multiple groups or not
        columns (Dict[str,str]): column name dictionary
        ignore_region (str): Region to be ignored in the analysis
        cred_set_file (str): Credible set filename
        ld_api (LDAccess): ld api object 
        extra_cols (List[str]): Extra columns to include in 
    Returns:
        (pd.DataFrame): Filtered and grouped variants
    """
    if ignore_region:
        ignore_region_=parse_region(ignore_region)
    if group and grouping_method == "cred":
        #load credsets
        cs_df=load_credsets(cred_set_file,columns)
        #remove ignored region if there is one
        if ignore_region:
            ign_idx=( ( cs_df[columns["chrom"]]==ignore_region_["chrom"] ) & ( cs_df[columns["pos"]]<=ignore_region_["end"] )&( cs_df[columns["pos"]]>=ignore_region_["start"] ) )
            cs_df=cs_df.loc[~ign_idx,:]

        #combine with summary stats to have all rows

        join_cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"]]
        cs_var_dtypes = {columns["chrom"]:str, columns["pos"]:int, columns["ref"]:str, columns["alt"]:str,columns["pval"]:float}
        cs_variants = load_groups_from_tabix(cs_df, gws_fpath, columns, "cs_id", columns,extract_cols = list(columns.values()) + extra_cols,dtypes=cs_var_dtypes)
        cs_variants = cs_df.merge(cs_variants,how="inner",on=join_cols) #join to add cs_df columns to cs_variants
        cs_variants.loc[:,"#variant"]=create_variant_column(cs_variants,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])

        #get ranges around credible set leads
        cs_lead_index = cs_df[["cs_id","cs_prob"]].reset_index().groupby("cs_id").max()["index"]
        cs_ranges = cs_df.loc[cs_lead_index, [columns["chrom"],columns["pos"]]]
        cs_ranges["min"] = cs_ranges[columns["pos"]]-locus_width*1000
        cs_ranges["max"] = cs_ranges[columns["pos"]]+locus_width*1000
        cs_ranges=cs_ranges.rename(columns={columns["chrom"]:"chrom"}).drop(columns=columns["pos"])
        
        #load summary stats around credsets, add columns for data
        summ_stat_variants = load_tb_ranges(cs_ranges,gws_fpath,"",".").drop_duplicates().astype(cs_var_dtypes)
        summ_stat_variants = extract_cols(summ_stat_variants,list(columns.values())+extra_cols)
        summ_stat_variants = summ_stat_variants.astype({columns["pos"]:int,columns["pval"]:float})
        summ_stat_variants.loc[:,"#variant"]=create_variant_column(summ_stat_variants,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        #remove variants that are already in the cs_variants
        summ_stat_variants = summ_stat_variants.loc[~summ_stat_variants["#variant"].isin(cs_variants["#variant"]),:]
        #concatenate the two 
        not_grouped_data = pd.concat([cs_variants, summ_stat_variants], axis="index", ignore_index=True, sort=False)\
            .sort_values(axis="index",by=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],"cs_id"],na_position="last")
        not_grouped_data=not_grouped_data.drop_duplicates(keep="first",subset=["#variant"])

        not_grouped_data=not_grouped_data.reset_index(drop=True)
        #not_grouped_data.loc[:,"#variant"]=create_variant_column(not_grouped_data,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        not_grouped_data.loc[:,"locus_id"]=not_grouped_data.loc[:,"#variant"]
        not_grouped_data.loc[:,"pos_rmax"]=not_grouped_data.loc[:,columns["pos"]]
        not_grouped_data.loc[:,"pos_rmin"]=not_grouped_data.loc[:,columns["pos"]]
        #group, sort data
        grouped_data = credible_set_grouping(data=not_grouped_data,ld_threshold=ld_r2,locus_range=locus_width,
                overlap=overlap,ld_api=ld_api,columns=columns)
        retval = grouped_data.sort_values(["locus_id","#variant"])
        
    else:
        sig_tresh_2=max(sig_tresh_1,sig_tresh_2)
        r=locus_width*1000#range for location width, originally in kb
        dtype={columns["chrom"]:str,
                    columns["pos"]:np.int32,
                    columns["ref"]:str,
                    columns["alt"]:str,
                    columns["pval"]:np.float64}

        #data input: get genome-wide significant variants.
        temp_df=get_gws_variants(gws_fpath,sign_treshold=sig_tresh_2,dtype=dtype,columns=columns,compression="gzip",extra_cols=extra_cols)
        #remove ignored region if there is one
        if ignore_region:
            ign_idx=( ( temp_df[columns["chrom"]]==ignore_region_["chrom"] ) & ( temp_df[columns["pos"]]<=ignore_region_["end"] )&( temp_df[columns["pos"]]>=ignore_region_["start"] ) )
            temp_df=temp_df.loc[~ign_idx,:]
        
        if temp_df.empty:
            print("The input file {} contains no gws-significant hits with signifigance treshold of {}. Aborting.".format(gws_fpath,sig_tresh_1))
            return None
        
        #data input: get credible set variants
        join_cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"]]
        if cred_set_file != "":
            cs_df=load_credsets(cred_set_file,columns)
        else:
            cs_df=pd.DataFrame(columns=join_cols+["cs_prob","cs_id"])
        #merge with gws_df, by using chrom,pos,ref,alt
        temp_df = merge_credset(temp_df,cs_df,gws_fpath,columns)
        #create necessary columns for the data
        temp_df=temp_df.reset_index(drop=True)
        temp_df.loc[:,"#variant"]=create_variant_column(temp_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        temp_df.loc[:,"locus_id"]=temp_df.loc[:,"#variant"]
        temp_df.loc[:,"pos_rmax"]=temp_df.loc[:,columns["pos"]]
        temp_df.loc[:,"pos_rmin"]=temp_df.loc[:,columns["pos"]]
        df_p1=temp_df.loc[temp_df[columns["pval"]] <= sig_tresh_1,: ].copy()
        df_p2=temp_df.loc[temp_df[columns["pval"]] <= sig_tresh_2,: ].copy()
        #grouping
        if group:
            if grouping_method=="ld":
                new_df=ld_grouping(df_p1=df_p1,df_p2=df_p2,
                sig_treshold_2=sig_tresh_2,locus_width=locus_width,ld_treshold=ld_r2,
                overlap=overlap,
                prefix=prefix,ld_api=ld_api,columns=columns)
            else :
                new_df=simple_grouping(df_p1=df_p1,df_p2=df_p2,r=r,overlap=overlap,columns=columns)
                new_df["r2_to_lead"]=np.nan
            new_df=new_df.sort_values(["locus_id","#variant"])
            retval = new_df
        else:
            #take only gws hits, no groups. Therefore, use df_p1
            retval = df_p1.sort_values(["locus_id","#variant"])
            retval["r2_to_lead"]=np.nan
        
    ## NOTE:Change X-chromosome to 23 in chrom, cs_id, #variant and locus_id
    retval = df_replace_value(retval,columns["chrom"],"X","23")
    retval = df_replace_value(retval,"cs_id",r'^chrX(.*)',r'chr23\1',regex=True)
    retval = df_replace_value(retval,"#variant",r'^chrX(.*)',r'chr23\1',regex=True)# NOTE: change 23 back in variant_id
    retval = df_replace_value(retval,"locus_id",r'^chrX(.*)',r'chr23\1',regex=True)
    return retval
    
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Fetch and group genome-wide significant variants from summary statistic")
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed summary statistic")
    parser.add_argument("--sign-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix. Default value is the base name (no path and no file extensions) of input file. ")
    parser.add_argument("--fetch-out",dest="fetch_out",type=str,default="fetch_out.tsv",help="GWS output filename, default is fetch_out.tsv")
    parser.add_argument("--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, options ['ld','simple','cred']")
    parser.add_argument("--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("--alt-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    parser.add_argument("--overlap",dest="overlap",action="store_true",help="Are groups allowed to overlap")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval","beta","maf","maf_cases","maf_controls"],help="Names for data file columns. Default is '#chrom pos ref alt pval beta maf maf_cases maf_controls'.")
    parser.add_argument("--extra-cols",dest="extra_cols",nargs="*",default=[],help="extra columns in the summary statistic you want to add to the results")
    parser.add_argument("--ignore-region",dest="ignore_region",type=str,default="",help="Ignore the given region, e.g. HLA region, from analysis. Give in CHROM:BPSTART-BPEND format.")
    parser.add_argument("--credible-set-file",dest="cred_set_file",type=str,default="",help="bgzipped SuSiE credible set file.")
    parser.add_argument("--ld-api",dest="ld_api_choice",type=str,default="plink",help="LD interface to use. Valid options are 'plink' and 'online'.")
    args=parser.parse_args()
    columns=columns_from_arguments(args.column_labels)
    if args.prefix!="":
        args.prefix=args.prefix+"."
    args.fetch_out = "{}{}".format(args.prefix,args.fetch_out)
    ld_api=None
    if args.ld_api_choice == "plink":
        ld_api = PlinkLD(args.ld_panel_path,args.plink_mem)
    elif args.ld_api_choice == "online":
        ld_api = OnlineLD("http://api.finngen.fi/api/ld")
    else:
        raise ValueError("Wrong argument for --ld-api:{}".format(args.ld_api_choice)) 
    fetch_df = fetch_gws(gws_fpath=args.gws_fpath, sig_tresh_1=args.sig_treshold, prefix=args.prefix, group=args.grouping, grouping_method=args.grouping_method, locus_width=args.loc_width,
        sig_tresh_2=args.sig_treshold_2, ld_r2=args.ld_r2, overlap=args.overlap, columns=columns,
        ignore_region=args.ignore_region, cred_set_file=args.cred_set_file,ld_api=ld_api, extra_cols=args.extra_cols)
    fetch_df.fillna("NA").replace("","NA").to_csv(path_or_buf=args.fetch_out,sep="\t",index=False,float_format="%.3g")