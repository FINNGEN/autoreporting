import argparse,shlex,subprocess, os
from subprocess import Popen, PIPE
import pandas as pd, numpy as np #typing: ignore
import tabix

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

def pytabix(tb,chrom,start,end):
    """Get genomic region from tabixed file
    In: pytabix handle, chromosome, start of region, end of region
    Out: list of variants in region 
    """
    try:
        retval=tb.querys("{}:{}-{}".format(chrom,start,end))
        return list(retval)
    except tabix.TabixError:
        return []

def create_variant_column(df,chrom="#chrom",pos="pos",ref="ref",alt="alt"):
    """Create 'chr$#chrom_$pos_$ref_$alt' column
    In: dataframe, with potentially defined column names
    Out: Variant column as pd.Series
    """
    if df.empty:
        return None
    return df.apply( lambda x: "chr{}_{}_{}_{}".format(x[chrom],x[pos],x[ref],x[alt]) ,axis=1)

def get_gzip_header(fname):
    """"Returns header for gzipped tsvs, as that is not currently possible using pytabix
    In: file path of gzipped tsv
    Out: header of tsv as a list of column names"""
    #create gzip call
    gzip_call=shlex.split("gzip -cd {}".format(fname))
    head_call=shlex.split("head -n 1")
    out=[]
    with Popen(gzip_call,stdout=PIPE) as gz:
        with Popen(head_call,stdin=gz.stdout,stdout=PIPE) as hd:
            for line in hd.stdout:
                out.append(line.decode().strip().split("\t"))
    return out[0]

def load_tb_df(df,fpath,chrom_prefix="",na_value=".",columns={"chrom":"#chrom"}):
    tb=tabix.open(fpath)
    tbxlst=[]
    for _,row in df.iterrows():
        tbxlst=tbxlst+pytabix( tb,"{}{}".format(chrom_prefix,row[columns["chrom"] ]),int(row[columns["pos"] ]),int(row[columns["pos"]]) )
    header=get_gzip_header(fpath)
    out_df=pd.DataFrame(tbxlst,columns=header )
    out_df=out_df.replace(na_value,np.nan)
    out_df[out_df.columns]=out_df[out_df.columns].apply(pd.to_numeric,errors="ignore")
    return out_df

def load_tb_ranges(df: pd.DataFrame, fpath: str, chrom_prefix: str = "", na_value: str = ".") -> pd.DataFrame:
    tb=tabix.open(fpath)
    tbxlst=[]
    for _,row in df.iterrows():
        tbxlst=tbxlst+pytabix( tb,"{}{}".format(chrom_prefix,row["chrom"]),int(row["min"]),int(row["max"]) )
    header=get_gzip_header(fpath)
    out_df=pd.DataFrame(tbxlst,columns=header )
    out_df=out_df.replace(na_value,np.nan)
    out_df[out_df.columns]=out_df[out_df.columns].apply(pd.to_numeric,errors="ignore")
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
    Out: Dictionary with the members 'chrom','pos','ref','alt','pval', 'beta', 'af' , 'af_case', 'af_control'
    """
    return {
        "chrom":column_labels[0],
        "pos":column_labels[1],
        "ref":column_labels[2],
        "alt":column_labels[3],
        "pval":column_labels[4],
        "beta":column_labels[5],
        "af":column_labels[6],
        "af_cases":column_labels[7],
        "af_controls":column_labels[8]}