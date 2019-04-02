#! /usr/bin/python3

import sys, argparse
import pandas as pd 
import numpy as np
#read dataframe, it's basically a tab delimited csv

def fetch_gws(args):
    fname=args.fpath
    sig_tresh=args.sig_treshold
    #get data
    c_size=10000
    result_dframe=None
    for dframe in pd.read_csv(fname,compression="gzip",sep="\t",dtype={"#chrom":str},chunksize=c_size):
        #filter gws snps
        if type(result_dframe)==type(None):
            result_dframe=dframe[dframe["pval"]<sig_tresh]
        else:
            result_dframe=result_dframe.append(dframe[dframe["pval"]<sig_tresh])
        
    #save gws snps to new vcf
    result_dframe.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Fetching gws SNPs from phenotype data")
    parser.add_argument("fpath",type=str,help="Filepath of the compressed tsv")
    parser.add_argument("-s","--signifigance-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("-o","--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    args=parser.parse_args()
    fetch_gws(args)