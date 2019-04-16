#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd 
import numpy as np

def fetch_gws(args):
    fname=args.fpath
    sig_tresh=args.sig_treshold

    c_size=100000
    result_dframe=None
    r=args.loc_width#range for location width
    #get data
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
    #print("filtered gws values")
    new_df=None
    if args.grouping:
        #basic idea:
        #we have our gws as result list, and tabixed list
        #while SNPs in result list:
        #   we find the most significant SNP in result list
        #   Then, wefind the SNPs from tabix file that are in its locus width area
        #   Remove those that are in result list, add tabixed ones to group, so groups can overlap.
        #   Repeat until no SNPs in result list.  
        tcall="tabix "+fname+" -h "
        for _idx,row in result_dframe.iterrows():
            tcall =tcall+" chr "+row["#chrom"]+":"+str(row["pos_rmin"])+"-"+str(row["pos_rmax"])
        #print(tcall) 
        call=shlex.split(tcall)
        #execute tabix call
        with open("temp.out","w") as out: 
            tbx=subprocess.run(call,stdout=out)
        #print("Tabix ran")
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
        tabixdf=tabixdf[tabixdf["pval"]<args.sig_treshold_2]
        #print("drop duplicates")
        tabixdf=tabixdf.drop_duplicates(subset=["#chrom","pos"])
        new_df=pd.DataFrame(columns=result_dframe.columns).drop(["pos_rmin","pos_rmax"],axis=1)
        #find groups
        i=1
        total=result_dframe.shape[0]
        while not result_dframe.empty:
            #find most significant SNP
            ms_snp=result_dframe.loc[result_dframe["pval"].idxmin(),:]
            #take the grop from tabixdf
            rowidx=(tabixdf["pos"]<=ms_snp["pos_rmax"])&(tabixdf["pos"]>=ms_snp["pos_rmin"])
            tmp=tabixdf.loc[rowidx,:].copy()
            tmp.loc[:,"locid"]=ms_snp["locid"]
            new_df=pd.concat([new_df,tmp],ignore_index=True,axis=0,join='inner')

            #convergence: remove the indexes from result_dframe
            dropidx=(result_dframe["pos"]<=ms_snp["pos_rmax"])&(result_dframe["pos"]>=ms_snp["pos_rmin"])
            result_dframe=result_dframe.loc[~dropidx,:]
            if i%10==0:
                print("iter: {}, SNPs remaining:{}/{}".format(i,result_dframe.shape[0],total))
            i+=1
        subprocess.run(shlex.split("rm temp.out"))
        new_df.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
    else:
        result_dframe.loc[:,["#chrom","pos","ref","alt","rsids",
                    "nearest_genes","pval","beta","sebeta","maf",
                    "maf_cases","maf_controls","locid"]].to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
    

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Fetching gws SNPs from phenotype data")
    parser.add_argument("fpath",type=str,help="Filepath of the compressed tsv")
    parser.add_argument("-s","--signifigance-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("-o","--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    parser.add_argument("-w","--locus-width",dest="loc_width",type=int,default=250000,help="location width to include for each SNP")
    parser.add_argument("-g" "--group", dest="grouping",action='store_true',help="Whether to include p-values that are within location width and have pval<sig_treshold_2")
    parser.add_argument("-s2","--alternate-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    args=parser.parse_args()
    fetch_gws(args)