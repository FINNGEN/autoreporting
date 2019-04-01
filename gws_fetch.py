#! /usr/bin/python3

import sys
import pandas as pd 
import numpy as np
#read dataframe, it's basically a tab delimited csv

def main():
    fname=sys.argv[1]
    #get data
    dframe=pd.read_csv(fname,compression="gzip",sep="\t",dtype={"#chrom":str})
    #filter gws snps
    gws_snp=dframe[dframe["pval"]<5e-8]
    #save gws snps to new vcf
    gws_snp.to_csv(path_or_buf="gws_snps.csv",sep="\t",index=False)


if __name__=="__main__":
    main()