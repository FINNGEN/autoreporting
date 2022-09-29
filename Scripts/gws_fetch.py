#!/usr/bin/env python3

import argparse,shlex,subprocess, glob
from subprocess import Popen, PIPE
import sys,os,io
import pandas as pd, numpy as np
import scipy.stats as stats
from typing import Dict, List
from autoreporting_utils import *
from data_access.linkage import PlinkLD, OnlineLD, Variant, LDData
from data_access.db import LDAccess
from data_access.db import CSAccess, CS, CSVariant
from data_access.cs import cs_to_df
from data_access.csfactory import csfactory

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



def ld_grouping(
    df_p1: pd.DataFrame,
    df_p2: pd.DataFrame,
    locus_range: int,
    dynamic_r2:bool,
    ld_threshold:float,
    overlap: bool,
    ld_api: LDAccess,
    columns: Dict[str,str]):
    """Group variants using LD
    Create groups using most significant variant as group lead, and LD partners with r2>ld_threshold as LD partners. Choose leads greedily with remaining min pval.
    Args:
        df_p1 (pd.DataFrame): Dataframe with variants that have pval>sig_threshold. These variants can be lead variants.
        df_p2 (pd.DataFrame): Dataframe with variants that have pval>sig_threshold_2. These variants can be LD partners.
        locus_range (int): Range around the locus on which LD is calculated, in bp
        dynamic_r2 (bool): Whether to adjust LD threshold by lead variant pval or not.
        ld_threshold (float): LD threshold for including a variant in the group. If dynamic r2 is False, it is the r^2 threshold. If dynamic r2 is True, it is ld threshold in ld_threshold/stats.chi2.isf(lead_pval,df=1)
        overlap (bool): Whether groups can overlap or not. If overlapping is set to True, variants grouped in a group are not removed from df_p2 after grouping, so other groups can claim them as well.
        ld_api (LDAccess): ld api object
        columns (Dict[str,str]): column dictionary
    Returns:
        (pd.DataFrame): Grouped variants in a pandas dataframe
    """
    leads = df_p1.copy()
    all_variants = df_p2.copy()
    #if r2 to lead data was gotten from cs, it's dropped as obsolete.
    if "r2_to_lead" in all_variants.columns:
        leads = leads.drop(columns=["r2_to_lead"])
        all_variants = all_variants.drop(columns=["r2_to_lead"])
    out_df = pd.DataFrame(columns=list(all_variants.columns)+["r2_to_lead"])
    iteration=0
    while not leads.empty:
        #get min pval variant
        lead_var_row = leads.loc[leads[columns["pval"]].idxmin(),: ]
        lead_var_id = lead_var_row["#variant"]
        #get ld threshold
        group_ld_threshold = ld_threshold
        if dynamic_r2:
            lead_pval = lead_var_row[columns["pval"]]
            group_ld_threshold = min(ld_threshold/stats.chi2.isf(lead_pval,df=1),1.0)
        #get LD neighbourhood
        ld_data = ld_api.get_range(Variant(lead_var_row[columns["chrom"]], lead_var_row[columns["pos"]], lead_var_row[columns["ref"]], lead_var_row[columns["alt"]]), locus_range, group_ld_threshold)
        flat_ld = [a.to_flat() for a in ld_data]
        ld_df = pd.DataFrame(flat_ld, columns=['chrom1','pos1','ref1','alt1','chrom2','pos2','ref2','alt2','r2'])
        ld_df = ld_df.rename(columns={"r2":"r2_to_lead"})
        if not ld_df.empty:
            ld_df['variant1']=create_variant_column(ld_df,'chrom1','pos1','ref1','alt1')
            ld_df['#variant']=create_variant_column(ld_df,'chrom2','pos2','ref2','alt2')
            ld_df = ld_df.drop(columns=["chrom1","pos1","chrom2","pos2","ref1","ref2","alt1","alt2"])
            ld_df = ld_df[ld_df["variant1"]==lead_var_id]
            merged_df = pd.merge(all_variants, ld_df[["#variant","r2_to_lead"]], how="inner",on="#variant")
            group_vars = merged_df.loc[:,out_df.columns]
            grouplead = all_variants[all_variants["#variant"]==lead_var_id].copy()
            grouplead["r2_to_lead"]=1.0
            group = pd.concat([group_vars,grouplead],ignore_index=True,axis=0,join='inner').drop_duplicates(subset=["#variant"])
        else:
            group = leads.loc[leads["#variant"]==lead_var_id,:]
            group["r2_to_lead"]=1.0
    
        group["locus_id"]=lead_var_id
        group["pos_rmin"]=group[columns["pos"]].min()
        group["pos_rmax"]=group[columns["pos"]].max()
        out_df = pd.concat([out_df,group],ignore_index=True,axis=0,join="inner")

        #remove all group variants from leads
        leads = leads[~leads["#variant"].isin(group["#variant"].unique())]
        #if overlap false, we remove grouped variants from all variants
        if not overlap:
            all_variants = all_variants[~all_variants["#variant"].isin(group["#variant"].unique())]
        iteration+=1
    return out_df


def credible_grouping(data: pd.DataFrame, dynamic_r2: bool, ld_threshold: float, locus_range: int, overlap: bool, ld_api: LDAccess, columns: Dict[str, str]) -> pd.DataFrame:
    """Group variants using credible sets
    Create groups using credible set most probable variants as the lead variants, and rest of the data as the additional variants
    Args:
        data (pd.DataFrame): Input data
        dynamic_r2 (bool): Whether to adjust LD threshold by lead variant pval or not.
        ld_threshold (float): LD threshold for including a variant in the group. If dynamic r2 is False, it is the r^2 threshold. If dynamic r2 is True, it is ld threshold in ld_threshold/stats.chi2.isf(lead_pval,df=1)
        locus_range (int): Range around the locus on which LD is calculated, in bp
        overlap (bool): Whether groups can overlap or not. If overlapping is set to True, variants grouped in a group are not removed from df_p2 after grouping, so other groups can claim them as well.
        ld_api (LDAccess): ld api object
        columns (Dict[str,str]): column dictionary
    Returns:
        (pd.DataFrame): Grouped variants in a pandas dataframe
    """
    df = data.copy()
    lead_vars = []
    for name, group in df.groupby(["cs_id"]):
        loc_id = "_".join(name.split("_")[:-1])#remove cs_number from cs_id
        group_lead =  group.loc[group["#variant"]==loc_id,"#variant"].iat[0]
        lead_vars.append(group_lead)
    if len(lead_vars) == 0:
        return pd.DataFrame(columns=df.columns)
    leads = df.loc[df["#variant"].isin(lead_vars)].copy()
    out_df = pd.DataFrame(columns=list(data.columns))

    while not leads.empty:
        lead_var_row = leads.loc[leads[columns["pval"]].idxmin(),: ]
        lead_variant = lead_var_row["#variant"] #choose the lead variant with smallest p-value
        cs_id = lead_var_row["cs_id"]
        #get ld threshold
        group_ld_threshold = ld_threshold
        if dynamic_r2:
            lead_pval = lead_var_row[columns["pval"]]
            group_ld_threshold = min(ld_threshold/stats.chi2.isf(lead_pval,df=1),1.0)
        #get LD data to cs lead
        ld_data = ld_api.get_range(Variant(lead_var_row[columns["chrom"]], lead_var_row[columns["pos"]], lead_var_row[columns["ref"]], lead_var_row[columns["alt"]]), locus_range, group_ld_threshold)
        flat_ld = [a.to_flat() for a in ld_data]
        ld_df = pd.DataFrame(flat_ld, columns=['chrom1','pos1','ref1','alt1','chrom2','pos2','ref2','alt2','r2'])
        ld_df = ld_df.rename(columns={"r2":"r2_to_lead"})
        if not ld_df.empty:
            ld_df['variant1']=create_variant_column(ld_df,'chrom1','pos1','ref1','alt1')
            ld_df['#variant']=create_variant_column(ld_df,'chrom2','pos2','ref2','alt2')
            ld_df = ld_df.drop(columns=["chrom1","pos1","chrom2","pos2","ref1","ref2","alt1","alt2"])
            ld_df = ld_df[ld_df["variant1"]==lead_variant]
            #merged_df = pd.merge(df, ld_df[["#variant","r2_to_lead"]], how="inner",on="#variant")
        else:
            #empty df, 
            ld_df = pd.DataFrame(columns=["#variant","variant1","r2_to_lead"])
        #separate credible set. It is deliberately taken from the not mutated 'data'-dataframe, so that even if those variants were grouped somewhere before, they are still included.
        cs = data.loc[data["cs_id"]==cs_id,:].copy()
        #fill r2 from the ld data if it's not all in the cs data
        cs_group = cs.merge(ld_df,on="#variant",how="left",suffixes = ("","_right"))#contains all variants of this CS, even though LD might be smaller than ld threshold
        cs_group["r2_to_lead"]=cs_group["r2_to_lead"].fillna(cs_group["r2_to_lead_right"])
        cs_group = cs_group.drop(columns=["r2_to_lead_right"])

        #get LD partners
        ld_partners = df.merge(ld_df[["#variant","r2_to_lead"]], on="#variant",how="inner",suffixes = ("_old",""))
        ld_partners = ld_partners.loc[ld_partners["cs_id"]!= cs_id,:]
        #fill with right
        ld_partners = ld_partners.drop(columns=["r2_to_lead_old"])
        #filter out the variants in the current credible set
        #remove credset data from LD partner variants
        wipe_credset_data = ["cs_prob",
            "cs_min_r2",
            "cs_log10bf",
            "good_cs",
            "cs_region",
            "cs_size",
            "cs_id"]
        for col in wipe_credset_data:
            ld_partners[col] = np.nan
        group=pd.concat([cs_group,ld_partners],ignore_index=True,sort=False).sort_values(by=["cs_id","#variant","r2_to_lead"]).drop_duplicates(subset=["#variant"],keep="first")
        group["locus_id"]=lead_variant
        group["pos_rmin"]=group[columns["pos"]].min()
        group["pos_rmax"]=group[columns["pos"]].max()
        out_df=pd.concat([out_df,group],ignore_index=True,axis=0,join='inner')
        
        #convergence: remove lead_variant, remove group from df if overlap is not true
        leads=leads[~ (leads["#variant"] == lead_variant) ]
        if not overlap:
            df=df[~df["#variant"].isin(ld_partners["#variant"])]
            df=df[~(df["cs_id"]==cs_id)]
        
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
    the rows corresponding to them are fetched using pysam.
    In: Dataframe containing gws variants, dataframe containing credible sets, filename for summary statistic.
    Out: Dataframe containing the gws variants + any credible set variants that are not gws. Columns 'cs_id','cs_prob' added to the dataframe. 
    """
    cols = list(gws_df.columns)
    join_cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"]]
    # fetch rows using pysam
    #create x and chr23 versions
    cs_f_x = df_replace_value(cs_df[join_cols].copy(),columns["chrom"],"X","23")
    cs_f_23 = df_replace_value(cs_df[join_cols].copy(),columns["chrom"],"23","X")
    #load both versions
    cred_df_x = load_pysam_df(cs_f_x,fname,columns,"",".")
    cred_df_23 = load_pysam_df(cs_f_23,fname,columns,"",".")
    #take the one which is longer
    cred_row_df = cred_df_23 if cred_df_23.shape[0] > cred_df_x.shape[0] else cred_df_x
    #turn chrX into chr23
    cred_row_df = df_replace_value(cred_row_df,columns["chrom"],"X","23")
    cred_row_df = cred_row_df[ cols ].drop_duplicates(keep="first")

    cred_row_df=cred_row_df.astype(dtype={columns["chrom"]:str,columns["pos"]:np.int64,columns["ref"]:str,columns["alt"]:str,columns["pval"]:float})
    cs_df=cs_df.astype(dtype={columns["chrom"]:str,columns["pos"]:np.int64,columns["ref"]:str,columns["alt"]:str})
    # ensure only the credible sets were included
    cred_row_df = pd.merge(cred_row_df,cs_df[join_cols],how="right",on=join_cols)
    #now, by concatting these two, we should have all of the cs variants summstat data present.
    df = pd.concat( [gws_df,cred_row_df], axis="index", ignore_index=True, sort=False)\
        .astype(dtype={columns["chrom"]:str,columns["pos"]:np.int64,columns["ref"]:str,columns["alt"]:str,columns["pval"]:float})\
        .drop_duplicates(subset=list( join_cols ) )
    # merge the credible set data to the summstat data
    merged = pd.merge(df,cs_df,how="left",on=join_cols)
    return merged

def fetch_gws(gws_fpath: str, 
                sig_tresh_1: float,
                prefix: str,
                group: bool,
                grouping_method: str,
                locus_width: int,
                sig_tresh_2: float,
                dynamic_r2: bool,
                ld_r2: float,
                overlap: bool,
                columns: Dict[str,str],
                ignore_region: str,
                cred_set_data: CSAccess,
                ld_api: LDAccess, 
                extra_cols: List[str],
                pheno_name: str,
                pheno_data_file :str):
    """Filter and group variants.
    Args:
        gws_fpath (str): summary statistic filename
        sig_tresh_1 (float): primary significance threshold (used in simple and ld grouping)
        prefix (str): prefix for all of the files
        group (bool): whether to group the variants or not
        grouping_method (str): grouping method
        locus_width (int): grouping width in kb
        sig_tresh_2 (float): alternate significance threshold
        ld_r2 (float): LD correlation threshold for including variants in groups (for cred   and ld grouping only)
        overlap (bool): Whether groups are allowed to overlap, i.e. same variant can be in multiple groups or not
        columns (Dict[str,str]): column name dictionary
        ignore_region (str): Region to be ignored in the analysis
        cred_set_file (str): Credible set filename
        ld_api (LDAccess): ld api object 
        extra_cols (List[str]): Extra columns to include in results
        pheno_name (str): Phenotype name
        pheno_data_file (str): Phenotype info file
    Returns:
        (pd.DataFrame): Filtered and grouped variants
    """
    locus_width_bp = locus_width*1000
    if ignore_region:
        ignore_region_=parse_region(ignore_region)
    if group and grouping_method == "cred":
        if not cred_set_data:
            raise Exception("--credible-set-file not specified")
        cs_data = cred_set_data.get_cs()
        cs_df = cs_to_df(cs_data,columns)
        #cs df X->23
        cs_df = df_replace_value(cs_df,columns["chrom"],"X","23")
        cs_df = df_replace_value(cs_df,"cs_id",r'^chrX(.*)',r'chr23\1',regex=True)
        cs_df = df_replace_value(cs_df,"cs_region",r'^chrX(.*)',r'chr23\1',regex=True)
        #remove ignored region if there is one
        if ignore_region:
            ign_idx=( ( cs_df[columns["chrom"]]==ignore_region_["chrom"] ) & ( cs_df[columns["pos"]]<=ignore_region_["end"] )&( cs_df[columns["pos"]]>=ignore_region_["start"] ) )
            cs_df=cs_df.loc[~ign_idx,:]
        if cs_df.empty:
            print("The input file {} contains no credible sets. Aborting.".format(gws_fpath))
            return None
        cs_leads = cs_df.loc[cs_df[["cs_id","cs_prob"]].reset_index().groupby("cs_id").max()["index"],:]
        cs_ranges = cs_leads[[columns["chrom"],columns["pos"]]].copy()
        cs_ranges["min"] = cs_ranges[columns["pos"]]-locus_width_bp
        cs_ranges["max"] = cs_ranges[columns["pos"]]+locus_width_bp
        cs_ranges=cs_ranges.rename(columns={columns["chrom"]:"chrom"}).drop(columns=columns["pos"])
        #make an X and 23 version of ranges
        cs_ranges_23 = df_replace_value(cs_ranges.copy(),"chrom","X","23")
        cs_ranges_x = df_replace_value(cs_ranges.copy(),"chrom","23","X")
        #load summary stats around credsets, add columns for data
        #use both cs_ranges objects, and choose the one that is longer
        summ_stat_variants_x = load_pysam_ranges(cs_ranges_x,gws_fpath,"",".")
        summ_stat_variants_23 = load_pysam_ranges(cs_ranges_23,gws_fpath,"",".")
        summ_stat_variants = summ_stat_variants_x if summ_stat_variants_x.shape[0] > summ_stat_variants_23.shape[0] else summ_stat_variants_23
        summ_stat_variants = df_replace_value(summ_stat_variants,columns["chrom"],"X","23")#finally in chrom 23 
        
        #filter them by p-value threshold 2. Merge credset will take care of cs variants that are filtered out.
        summ_stat_variants= summ_stat_variants.loc[summ_stat_variants[columns["pval"]]<=sig_tresh_2,:]
        
        summ_stat_variants = extract_cols(summ_stat_variants,list(columns.values())+extra_cols)

        #check that every CS variant is in summ stat variants
        #check that all cs vars are in cred_row_df
        cs_df_x = df_replace_value(cs_df.copy(),columns["chrom"],"23","X")
        cs_df_23 = df_replace_value(cs_df.copy(),columns["chrom"],"X","23")
        cs_df_ss_x = load_pysam_df(cs_df_x,gws_fpath,columns,"",".")
        cs_df_ss_23 = load_pysam_df(cs_df_23,gws_fpath,columns,"",".")
        cs_df_ss = cs_df_ss_x if cs_df_ss_x.shape[0] > cs_df_ss_23.shape[0] else cs_df_ss_23
        cs_df_ss = df_replace_value(cs_df_ss,columns["chrom"],"X","23")
        summstat_var_col = create_variant_column(cs_df_ss,columns["chrom"],columns["pos"],columns["ref"],columns["alt"])
        cs_var_col = create_variant_column(cs_df,columns["chrom"],columns["pos"],columns["ref"],columns["alt"])
        if not all(cs_var_col.isin(summstat_var_col)):
            print("ERROR: not all CS variants were in summary statistic")
            print([cs_var for cs_var in cs_var_col if cs_var not in summstat_var_col])
            raise Exception("Not all cs variants were in summary statistic file. Grouping can not continue.")

        not_grouped_data = merge_credset(summ_stat_variants,cs_df,gws_fpath,columns)\
            .sort_values(axis="index",by=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],"cs_id"],na_position="last")
        not_grouped_data=not_grouped_data.reset_index(drop=True)
        not_grouped_data.loc[:,"#variant"]=create_variant_column(not_grouped_data,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        not_grouped_data.loc[:,"locus_id"]=not_grouped_data.loc[:,"#variant"]
        not_grouped_data.loc[:,"pos_rmax"]=not_grouped_data.loc[:,columns["pos"]]
        not_grouped_data.loc[:,"pos_rmin"]=not_grouped_data.loc[:,columns["pos"]]
        #group, sort data
        grouped_data = credible_grouping(
            data=not_grouped_data,
            dynamic_r2=dynamic_r2,
            ld_threshold=ld_r2,
            locus_range=locus_width_bp,
            overlap=overlap,
            ld_api=ld_api,
            columns=columns)
        retval = grouped_data.sort_values(["locus_id","#variant"])
        
    else:
        sig_tresh_2=max(sig_tresh_1,sig_tresh_2)
        dtype={columns["chrom"]:str,
                    columns["pos"]:np.int32,
                    columns["ref"]:str,
                    columns["alt"]:str,
                    columns["pval"]:np.float64}

        #data input: get genome-wide significant variants.
        temp_df=get_gws_variants(gws_fpath,sign_treshold=sig_tresh_2,dtype=dtype,columns=columns,compression="gzip",extra_cols=extra_cols)
        temp_df = df_replace_value(temp_df,columns["chrom"],"X","23")#summ stat data to 23 if not there already
        #remove ignored region if there is one
        if ignore_region:
            ign_idx=( ( temp_df[columns["chrom"]]==ignore_region_["chrom"] ) & ( temp_df[columns["pos"]]<=ignore_region_["end"] )&( temp_df[columns["pos"]]>=ignore_region_["start"] ) )
            temp_df=temp_df.loc[~ign_idx,:]
        
        if temp_df.empty:
            print("The input file {} contains no gws-significant hits with signifigance treshold of {}. Aborting.".format(gws_fpath,sig_tresh_1))
            return None
        
        #data input: get credible set variants
        join_cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"]]
        if cred_set_data:
            cs_data = cred_set_data.get_cs()
            cs_df = cs_to_df(cs_data,columns)
        else:
            cs_df=pd.DataFrame(columns=join_cols+["cs_prob","cs_id","cs_region"])
        cs_df = df_replace_value(cs_df,columns["chrom"],"X","23")
        cs_df = df_replace_value(cs_df,"cs_id",r'^chrX(.*)',r'chr23\1',regex=True)
        cs_df = df_replace_value(cs_df,"cs_region",r'^chrX(.*)',r'chr23\1',regex=True)
        #merge with gws_df, by using chrom,pos,ref,alt
        temp_df = merge_credset(temp_df,cs_df,gws_fpath,columns)
        #create necessary columns for the data
        temp_df=temp_df.reset_index(drop=True)
        temp_df.loc[:,"#variant"]=create_variant_column(temp_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        temp_df.loc[:,"locus_id"]=temp_df.loc[:,"#variant"]
        temp_df.loc[:,"pos_rmax"]=temp_df.loc[:,columns["pos"]]
        temp_df.loc[:,"pos_rmin"]=temp_df.loc[:,columns["pos"]]
        if "r2_to_lead" in temp_df.columns: #cs r2 to lead is not useful since we don't necessarily group around cs variants
            temp_df=temp_df.drop(columns=["r2_to_lead"])
        df_p1=temp_df.loc[temp_df[columns["pval"]] <= sig_tresh_1,: ].copy()
        df_p2=temp_df.loc[temp_df[columns["pval"]] <= sig_tresh_2,: ].copy()
        if df_p1.empty:
            print("The input file {} contains no gws-significant hits with signifigance treshold of {}. Aborting.".format(gws_fpath,sig_tresh_1))
            return None
        #grouping
        if group:
            if grouping_method=="ld":
                new_df=ld_grouping(
                    df_p1=df_p1,
                    df_p2=df_p2,
                    locus_range=locus_width_bp,
                    dynamic_r2 = dynamic_r2,
                    ld_threshold=ld_r2,
                    overlap=overlap,
                    ld_api=ld_api,
                    columns=columns
                )
            else :
                new_df=simple_grouping(df_p1=df_p1,df_p2=df_p2,r=locus_width_bp,overlap=overlap,columns=columns)
                new_df["r2_to_lead"]=np.nan
            new_df=new_df.sort_values(["locus_id","#variant"])
            retval = new_df
        else:
            #take only gws hits, no groups. Therefore, use df_p1
            retval = df_p1.sort_values(["locus_id","#variant"])
            retval["r2_to_lead"]=np.nan

    #add phenotype name
    retval["phenotype"] = pheno_name
    #add phenotype data
    #load phenotype datafile
    try:
        pheno_data = pd.read_csv(pheno_data_file,sep="\t")
        pheno_row = pheno_data[pheno_data["phenocode"] == pheno_name].iloc[0]
        retval["longname"] = pheno_row["name"]
        retval["category"] = pheno_row["category"]
        retval["n_cases"] = pheno_row["num_cases"]
        retval["n_controls"] = pheno_row["num_controls"]
    except:
        print(f"Error in reading phenotype {pheno_name} data from phenotype datafile {pheno_data_file}. Either the file did not exist, or the endpoint was not present in the file. Filling the values with NA")
        retval["longname"] = np.nan
        retval["category"] = np.nan
        retval["n_cases"] = np.nan
        retval["n_controls"] = np.nan

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
    #r2 static bound or dynamic per peak
    r2_group = parser.add_mutually_exclusive_group()
    r2_group.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    r2_group.add_argument("--dynamic-r2-chisq",type=float,nargs="?",const=5.0,default=None,help="If flag is passed, r2 threshold is set per peak so that leadvar_chisq*r2=value (default 5).")
    
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    parser.add_argument("--overlap",dest="overlap",action="store_true",help="Are groups allowed to overlap")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval","beta","maf","maf_cases","maf_controls"],help="Names for data file columns. Default is '#chrom pos ref alt pval beta maf maf_cases maf_controls'.")
    parser.add_argument("--extra-cols",dest="extra_cols",nargs="*",default=[],help="extra columns in the summary statistic you want to add to the results")
    parser.add_argument("--ignore-region",dest="ignore_region",type=str,default="",help="Ignore the given region, e.g. HLA region, from analysis. Give in CHROM:BPSTART-BPEND format.")
    parser.add_argument("--credible-set-file",dest="cred_set_file",type=str,default="",help="bgzipped SuSiE credible set file.")
    parser.add_argument("--ld-api",dest="ld_api_choice",type=str,default="plink",help="LD interface to use. Valid options are 'plink' and 'online'.")
    parser.add_argument("--pheno-name",dest="pheno_name",type=str,default="",help="Phenotype name")
    parser.add_argument("--pheno-info-file",dest="pheno_info_file",type=str,default="",help="Phenotype information file path")
    args=parser.parse_args()
    columns=columns_from_arguments(args.column_labels)
    if args.prefix!="":
        args.prefix=args.prefix+"."
    args.fetch_out = "{}{}".format(args.prefix,args.fetch_out)
    #ld api loading
    ld_api=None
    if args.ld_api_choice == "plink":
        ld_api = PlinkLD(args.ld_panel_path,args.plink_mem)
    elif args.ld_api_choice == "online":
        ld_api = OnlineLD("http://api.finngen.fi/api/ld")
    else:
        raise ValueError("Wrong argument for --ld-api:{}".format(args.ld_api_choice)) 
    #cs access object loading
    if args.cred_set_file:
        cs_access = csfactory(args.cred_set_file)
    else:
        cs_access=None

    #parse r2 
    if args.dynamic_r2_chisq != None:
        dynamic_r2 = True
        r2_thresh =  args.dynamic_r2_chisq
    else:
        dynamic_r2 = False
        r2_thresh = args.ld_r2

    fetch_df = fetch_gws(
        gws_fpath=args.gws_fpath,
        sig_tresh_1=args.sig_treshold,
        prefix=args.prefix,
        group=args.grouping,
        grouping_method=args.grouping_method,
        locus_width=args.loc_width,
        sig_tresh_2=args.sig_treshold_2,
        dynamic_r2=dynamic_r2,
        ld_r2=r2,
        overlap=args.overlap,
        columns=columns,
        ignore_region=args.ignore_region,
        cred_set_data=cs_access,
        ld_api=ld_api,
        extra_cols=args.extra_cols,
        pheno_name=args.pheno_name,
        pheno_data_file =args.pheno_info_file
    )
    fetch_df.fillna("NA").replace("","NA").to_csv(path_or_buf=args.fetch_out,sep="\t",index=False,float_format="%.3g")