#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd 
import numpy as np
#read dataframe, it's basically a tab delimited csv

def fetch_gws(args):
    fname=args.fpath
    sig_tresh=args.sig_treshold
    #get data
    c_size=100000
    result_dframe=None
    r=args.loc_width#range for location width
    for dframe in pd.read_csv(fname,compression="gzip",sep="\t",dtype={"#chrom":str,
        "pos":np.int32,
        "ref":str,
        "alt":str,
        "rsids":str,
        "nearest_genes":str,
        "pval":np.float64,
        "beta":np.float64,
        "sebeta":np.float64,
        "maf":np.float64,
        "maf_cases":np.float64,
        "maf_controls":np.float64
        },engine="c",chunksize=c_size):
        #filter gws snps
        temp_dframe=dframe[dframe["pval"]<sig_tresh].copy()
        if not temp_dframe.empty:
            temp_dframe.loc[:,"zero"]=0
            temp_dframe.loc[:,"pos_rmin"]=temp_dframe.loc[:,("pos","zero") ].max(axis=1)
            temp_dframe.loc[:,"pos_rmax"]=temp_dframe.loc[:,"pos"]+r
            #drop zero
            temp_dframe=temp_dframe.drop(axis=1,columns=["zero"])
            temp_dframe.loc[:,"locid"]="c_"+temp_dframe["#chrom"].map(str)+"_p_"+temp_dframe["pos"].map(str)
            if type(result_dframe) == type(None):
                result_dframe=temp_dframe.copy()
            else:
                result_dframe.append(temp_dframe)
    print("filtered gws values")
    new_df=None
    if args.grouping:
        #PROBLEM
        #what if we search overlapping areas? setting the group will set them erratically, and give wrong results
        #if we think of the tabix-fetched rows as one list of rows (worse, we will have duplicate rows, and how do we dowrk with those?
        # filter them out?) Oh well, most likely we won't come across this in the beginning, as the rows are filtered based on the secondary
        #pval, but it might happen.
        #First, make tabix calls and write them to a file, then fetch all the relevant rows
        #so, we need the chromosome and the min and max locations
        tcall="tabix "+fname+" -h "
        for _idx,row in result_dframe.iterrows():
            tcall =tcall+" chr "+row["#chrom"]+":"+str(row["pos_rmin"])+"-"+str(row["pos_rmax"])
        #print(tcall) 
        call=shlex.split(tcall)
        #execute tabix call
        with open("temp.out","w") as out: 
            tbx=subprocess.run(call,stdout=out)
        #load ouput file
        print("Tabix ran")
        tabixdf=pd.read_csv("temp.out",sep="\t",dtype={"#chrom":str,
        "pos":np.int32,
        "ref":str,
        "alt":str,
        "rsids":str,
        "nearest_genes":str,
        "pval":np.float64,
        "beta":np.float64,
        "sebeta":np.float64,
        "maf":np.float64,
        "maf_cases":np.float64,
        "maf_controls":np.float64
        },engine="c")

        #then, we need to 
        # 1)filter them according to p-value
        tabixdf=tabixdf[tabixdf["pval"]<args.sig_treshold_2]
        #remove duplicates from tabixdf
        print("drop duplicates")
        tabixdf=tabixdf.drop_duplicates(subset=["#chrom","pos"])
        # 2) attach the correct groups, I guess
        #this might be harder to do.
        new_df=pd.DataFrame(columns=result_dframe.columns).drop(["pos_rmin","pos_rmax"],axis=1)
        for t in result_dframe.itertuples():
            rowidx=(tabixdf["pos"]<=t.pos_rmax)&(tabixdf["pos"]>=t.pos_rmin)
            tmp=tabixdf.loc[rowidx,:].copy()
            tmp.loc[:,"locid"]=t.locid
            #print(tmp)
            new_df=pd.concat([new_df,tmp],ignore_index=True,axis=0,join='inner',copy=False)
            #tabixdf.loc[rowidx,"locid"]=t.locid
        # it makes no sense to call tabix many times just because it would be easier to attach a group id that way
        #instead, we might iterate over the resulting rows? It might take time but shouldn't be too slow, since the 
        #dataframe isn't that big anymore (maybe a few million rows max)
    if args.grouping:
        new_df.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
    else:
        result_dframe.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
    #remove temp.out
    subprocess.run(shlex.split("rm temp.out"))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Fetching gws SNPs from phenotype data")
    parser.add_argument("fpath",type=str,help="Filepath of the compressed tsv")
    parser.add_argument("-s","--signifigance-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("-o","--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    parser.add_argument("-w","--locus-width",dest="loc_width",type=int,default=250000,help="location width to include for each SNP")
    parser.add_argument("-g" "--group", dest="grouping",action='store_true',help="Whether to include p-values that are within location width and have pval<sig_treshold_2")
    parser.add_argument("-s2","--alternate-sign-treshold",dest="sig_treshold_2",type=float, default=5e-7,help="optional group treshold")
    args=parser.parse_args()
    fetch_gws(args)