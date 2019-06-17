import argparse,shlex,subprocess
from subprocess import Popen, PIPE
import pandas as pd, numpy as np
import tabix

"""
Utility functions that are used in the scripts, put here for keeping the code clearer
"""

def pytabix(tb,chrom,start,end):
    try:
        retval=tb.querys("{}:{}-{}".format(chrom,start,end))
        return retval
    except tabix.TabixError:
        return []

def create_variant_column(df,chrom="#chrom",pos="pos",ref="ref",alt="alt"):
    """Create 'chr$#chrom_$pos_$ref_$alt' column
    In: dataframe, with potentially defined column names
    Out: Variant column as pd.Series
    """
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
                if ((t.pos_rmax<region["min"]) or (t.pos_rmin>region["max"]) ) or (t._1!=region["#chrom"]):
                    continue
                elif t.pos_rmin>=region["min"] and t.pos_rmax<=region["max"]:
                    found=True
                    break
                else: 
                    if t.pos_rmin<region["min"] and t.pos_rmax>region["min"]:
                        region["min"]=t.pos_rmin
                    if t.pos_rmax>region["max"] and t.pos_rmin<region["max"]:
                        region["max"]=t.pos_rmax
                    found=True
                    break
            if not found:
                regions.append({"#chrom":t._1,"min":t.pos_rmin,"max":t.pos_rmax})
        else:
            regions.append({"#chrom":t._1,"min":t.pos_rmin,"max":t.pos_rmax})
    return pd.DataFrame(regions)