import argparse,shlex,subprocess, os
from subprocess import Popen, PIPE
import pandas as pd, numpy as np #typing: ignore
import tabix
from typing import List, Dict, Any, Optional
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

def load_tb_df(df,fpath,columns,chrom_prefix="",na_value="."):
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
    """Load variant ranges from a file using tabix
    Args:
        df (pd.DataFrame): range dataframe with columns chrom, min, max
        fpath (str): tabixed file path
        chrom_prefix (str): Whether the tabix sequence has a prefix, e.g. 'chr10' instead of '10'
        na_value (str): NA representation that gets filled with NaN
    Returns:
        (pd.DataFrame): Dataframe with variants from the specified ranges
    """
    tb=tabix.open(fpath)
    tbxlst=[]
    for _,row in df.iterrows():
        tbxlst=tbxlst+pytabix( tb,"{}{}".format(chrom_prefix,row["chrom"]),int(row["min"]),int(row["max"]) )
    header=get_gzip_header(fpath)
    out_df=pd.DataFrame(tbxlst,columns=header )
    out_df=out_df.replace(na_value,np.nan)
    return out_df

def load_tb_ranges_spec_cols(df:pd.DataFrame, fpath: str, extract_cols: List[str], chrom_prefix: str = "", na_value: str = ".") -> pd.DataFrame:
    """Load variant ranges from a file using tabix, but limit the columns to just the specified ones to save memory
    Args:
        df (pd.DataFrame): range dataframe with columns chrom, min, max
        fpath (str): tabixed file path
        extract_cols (List[str]): List of columns to extract from the file
        chrom_prefix (str): Whether the tabix sequence has a prefix, e.g. 'chr10' instead of '10'
        na_value (str): NA representation that gets filled with NaN
    Returns:
        (pd.DataFrame): Dataframe with variants from the specified ranges, with the specified columns
    """
    tb=tabix.open(fpath)
    header=get_gzip_header(fpath)
    dflist=[]
    for _, row in df.iterrows():
        region_df = pytabix( tb,"{}{}".format(chrom_prefix,row["chrom"]),int(row["min"]),int(row["max"]) )
        region_df=pd.DataFrame(region_df,columns=header)
        region_df=region_df.loc[:,extract_cols]
        dflist.append(region_df)
    output_df = pd.concat(dflist,ignore_index=True).replace(na_value, np.nan)
    return output_df

def load_groups_from_tabix(df: pd.DataFrame, fpath: str, columns: Dict[str, str], group_column: str, data_columns: Optional[Dict[str, str]] = None, chrom_prefix: Optional[str] = "", na_value: Optional[str] = ".",extract_cols: Optional[List[str]] = None) -> pd.DataFrame:
    """Load dataframe's groups from another tabixed file
    Args:
        df (pd.DataFrame): input dataframe
        fpath (str): data filepath
        column (Dict[str, str]): column name dictionary
        group_column (str): group column name
        data_columns (Optional[Dict[str, str]]): Optional data column names, if different from the ones in columns
        chrom_prefix (Optional[str]): chromosome prefix, default ""
        na_value (Optional[str]): na representation in data file, default "."
        extract_cols (Optional[List[str]]): columns to extract from dataframe.
    Returns:
        (pd.DataFrame): Dataframe with the groups from the input dataframe. Contains extra variants.
    """
    if not data_columns:
        data_columns = columns
    groups = df[[group_column,columns["chrom"],columns["pos"]]]\
        .rename(columns={columns["chrom"]:"chrom",
                         columns["pos"]:"pos"})\
        .groupby(group_column)\
        .aggregate(chrom = pd.NamedAgg(column="chrom",aggfunc=pd.Series.mode),
                   max = pd.NamedAgg(column="pos", aggfunc = max),
                   min = pd.NamedAgg(column="pos", aggfunc = min))
    groups = prune_df_regions(groups)
    groups.to_csv("debug_ranges.tsv",sep="\t",index=False)
    if not extract_cols:
        output = load_tb_ranges(df=groups, fpath=fpath,chrom_prefix=chrom_prefix, na_value=na_value)
    else:
        output = load_tb_ranges_spec_cols(df=groups, fpath=fpath, extract_cols=extract_cols, chrom_prefix=chrom_prefix, na_value=na_value)
    output = output.rename(columns={data_columns["chrom"]:columns["chrom"],
                            data_columns["pos"]:columns["pos"],
                            data_columns["ref"]:columns["ref"],
                            data_columns["alt"]:columns["alt"]})
    if chrom_prefix != "":
        output[columns["chrom"]] = output[columns["chrom"]].apply(lambda x: x.strip(chrom_prefix))
    output = output.drop_duplicates()
    join_cols = [columns["chrom"], columns["pos"], columns["ref"], columns["alt"]]
    join_types = {columns["chrom"]:str, columns["pos"]:int, columns["ref"]:str, columns["alt"]:str}
    output = output.astype(join_types)
    df = df.astype(join_types)
    output = output.merge(df[join_cols],how="right",on=join_cols)
    return output

def prune_df_regions(df: pd.DataFrame, cols: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Prune dataframe regions
    Args:
        df (pd.DataFrame): dataframe with columns for sequence, min, max
        cols (Optional[Dict[str, str]]): Optional column names
    Returns:
        (pd.DataFrame): dataframe with same columns, with overlapping regions merged
    """
    if not cols:
        cols = {"chrom":"chrom","min":"min","max":"max"}
    df = df.copy()
    df=df.rename(columns={v: k for k,v in cols.items()})#rename to chrom, min, max so tuples naming limitations do not apply
    regions = []
    for t in df.itertuples():
        if regions:
            found=False
            for region in regions:
                if (( int(t.max)<region[cols["min"]] ) or ( int(t.min) > region[cols["max"]] ) ) or (str(t.chrom)!=region[ cols["chrom"] ] ):
                    continue
                elif int(t.min)>=region[cols["min"]] and int(t.max)<=region[cols["max"]]:
                    found=True
                    break
                else: 
                    if int(t.min)<region[cols["min"]] and int(t.max)>region[cols["min"]]:
                        region[cols["min"]]=int(t.min)
                    if int(t.max)>region[cols["max"]] and int(t.min)<region[cols["max"]]:
                        region["max"]=int(t.max)
                    found=True
                    break
            if not found:
                regions.append({"chrom":str(t.chrom),"min":int(t.min),"max":int(t.max)})
        else:
            regions.append({"chrom":str(t.chrom),"min":int(t.min),"max":int(t.max)})
    return pd.DataFrame(regions).rename(columns=cols) #rename back

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