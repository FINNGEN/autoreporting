import argparse,shlex,subprocess, os
from subprocess import Popen, PIPE
import pandas as pd, numpy as np #typing: ignore
import tabix
import pysam
from typing import List, Dict
"""
Utility functions that are used in the scripts, put here for keeping the code clearer
"""

def filebasename(s):
    if s != "":
        return os.path.basename(s).split(".")[0]
    return ""

def df_replace_value(df,column,value,replace_with,regex=False):
    """Replace value on column with values

    Args:
        df (pd.DataFrame): Pandas Dataframe with at least column column
        column (str): Dictionary with enytries for the column names of dataframe
        value (Any): value to replace
        replace_with (Any): value to replace value with
    Returns:
        (pd.DataFrame):Dataframe with value replaced
    """
    df[column] = df[column].replace(value,replace_with,regex=regex)
    return df


def create_variant_column(df,chrom="#chrom",pos="pos",ref="ref",alt="alt"):
    """Create 'chr$#chrom_$pos_$ref_$alt' column
    In: dataframe, with potentially defined column names
    Out: Variant column as pd.Series
    """
    if df.empty:
        return None
    return df.apply( lambda x: "chr{}_{}_{}_{}".format(x[chrom],x[pos],x[ref],x[alt]) ,axis=1)

def load_annotation_df(df: pd.DataFrame, fpath: str, columns: Dict[str,str], resource_columns: Dict[str,str], chrom_prefix: str="", na_value: str=".") -> pd.DataFrame:
    """Load annotation data by reading the whole annotation file
    This is slower than load_pysam_df for phenotypes with little results, but massively faster for phenotypes with a lot of results (>40k rows)
    Also, the time is not that dependent on input size, which is a nice bonus. 
    """
    df_colsubset = [columns["chrom"],columns["pos"],columns["ref"],columns["alt"]]
    data_to_load = df.drop_duplicates(subset=df_colsubset)
    data_to_load = data_to_load[df_colsubset]
    data_to_load[columns["chrom"]]=data_to_load[columns["chrom"]].apply(lambda x:chrom_prefix+x)
    dtype = {resource_columns["chrom"]:str,
        resource_columns["pos"]:np.int32,
        resource_columns["ref"]:str,
        resource_columns["alt"]:str
    }
    chunksize = 1000*1000
    chr_set = set(data_to_load[columns["chrom"]])
    pos_set = set(data_to_load[columns["pos"]])
    ref_set = set(data_to_load[columns["ref"]])
    alt_set = set(data_to_load[columns["alt"]])
    out=pd.DataFrame()
    for partial_df in pd.read_csv(fpath, compression="gzip",sep="\t",engine="c",dtype=dtype,chunksize=chunksize,na_values=na_value):
        bool_condition = partial_df[resource_columns["chrom"]].isin(chr_set) &\
            partial_df[resource_columns["pos"]].isin(pos_set) &\
            partial_df[resource_columns["ref"]].isin(ref_set) &\
            partial_df[resource_columns["alt"]].isin(alt_set)
        out = pd.concat([out,partial_df.loc[bool_condition,:]],axis="index",ignore_index=True,sort=False)
    out[out.columns]=out[out.columns].apply(pd.to_numeric,errors="ignore")
    return out

def load_pysam_df(df,fpath,columns,chrom_prefix="",na_value=".") -> pd.DataFrame:
    """Load variants using pysam from tabix-indexed file
    Args:

    Returns:
        (pd.DataFrame): DataFrame with same columns as the annotation file
    """
    tb = pysam.TabixFile(fpath)
    #generate chrom-position df
    chrompos_df = df[[columns["chrom"], columns["pos"]]].copy()
    chrompos_df = chrompos_df.rename(columns = {columns["chrom"]:"chrom",columns["pos"]:"pos"}).drop_duplicates(keep="first")
    tbxlst = []
    for t in chrompos_df.itertuples():
        #get rows
        try:
            rows = tb.fetch("{}{}".format(chrom_prefix,t.chrom), int(t.pos)-1,int(t.pos))
        except:
            rows = []
        data = [a.strip('\n').split('\t') for a in rows]
        tbxlst.extend(data)
    header = tb.header[0].split('\t')
    out_df = pd.DataFrame(tbxlst, columns = header)
    out_df=out_df.replace(na_value,np.nan)
    out_df[out_df.columns]=out_df[out_df.columns].apply(pd.to_numeric,errors="ignore")
    tb.close()
    return out_df

def load_pysam_ranges(df: pd.DataFrame, fpath: str, chrom_prefix: str = "", na_value: str = ".") -> pd.DataFrame:
    tb = pysam.TabixFile(fpath)
    tbxlst=[]
    for _,row in df.iterrows():
        try:
            rows = tb.fetch("{}{}".format(chrom_prefix,row["chrom"]),int(row["min"])-1,int(row["max"]))
        except:
            rows = []
        data = [a.strip('\n').split('\t') for a in rows]
        tbxlst.extend(data)
    header = tb.header[0].split('\t')
    out_df = pd.DataFrame(tbxlst, columns = header)
    out_df=out_df.replace(na_value,np.nan)
    out_df[out_df.columns]=out_df[out_df.columns].apply(pd.to_numeric,errors="ignore")
    tb.close()
    return out_df

def prune_regions(df):
    """Prune overlapping tabix regions so that no duplicate calls are made
    In: dataframe with the regions
    Out:dataframe with pruned regions
    """
    regions=[]
    for t in df.itertuples():
        if regions:
            found=False
            for region in regions:
                if (( int(t.pos_rmax)<region["min"] ) or ( int(t.pos_rmin) > region["max"] ) ) or (str(t._1)!=region[ "chrom" ] ):
                    continue
                elif int(t.pos_rmin)>=region["min"] and int(t.pos_rmax)<=region["max"]:
                    found=True
                    break
                else: 
                    if int(t.pos_rmin)<region["min"] and int(t.pos_rmax)>region["min"]:
                        region["min"]=int(t.pos_rmin)
                    if int(t.pos_rmax)>region["max"] and int(t.pos_rmin)<region["max"]:
                        region["max"]=int(t.pos_rmax)
                    found=True
                    break
            if not found:
                regions.append({"chrom":str(t._1),"min":int(t.pos_rmin),"max":int(t.pos_rmax)})
        else:
            regions.append({"chrom":str(t._1),"min":int(t.pos_rmin),"max":int(t.pos_rmax)})
    return pd.DataFrame(regions)

def columns_from_arguments(column_labels):
    """
    Return a dict of columns (used pervasively throughout the script) from the argument column_labels
    In: column labels, as a list
    Out: Dictionary with the members 'chrom','pos','ref','alt','pval'
    """
    return {
        "chrom":column_labels[0],
        "pos":column_labels[1],
        "ref":column_labels[2],
        "alt":column_labels[3],
        "pval":column_labels[4]}