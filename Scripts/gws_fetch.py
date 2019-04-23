#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd 
import numpy as np

def fetch_gws(args):
    fname=args.fpath
    sig_tresh=args.sig_treshold

    c_size=100000
    r=args.loc_width#range for location width
    df=None
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
        temp_dframe=dframe.loc[dframe["pval"]<sig_tresh,:]
        if not temp_dframe.empty:
            if type(df) == type(None):
                df=temp_dframe.copy()
            else:
                df=df.append(temp_dframe)
    df=df.reset_index(drop=True)
    df.to_csv("df_debug.csv",sep="\t",index=False)
    df.loc[:,"pos_rmin"]=df.loc[:,"pos"]-r
    df.loc[:,"pos_rmax"]=df.loc[:,"pos"]+r
    df.loc[:,"pos_rmin"]=df.loc[:,"pos_rmin"].clip(lower=0)
    df.loc[:,"locid"]="c_"+df.loc[:,"#chrom"].map(str)+"_p_"+df.loc[:,"pos"].map(str)
    result_dframe=df
    new_df=None
    if args.grouping:
        #basic idea:
        #we have our gws as result list, and tabixed list
        #while SNPs in result list:
        #   we find the most significant SNP in result list
        #   Then, wefind the SNPs from tabix file that are in its locus width area
        #   Remove those that are in result list, add tabixed ones to group, so groups can overlap.
        #   Repeat until no SNPs in result list.  

        #simplify tabix call by removing overlapping tabix calls
        regions=[]
        for t in df.loc[:,["#chrom", "pos_rmin","pos_rmax"]].itertuples():
            if regions:
                #do region stuff
                found=False
                for region in regions:
                    if (t.pos_rmax<region["min"]) or (t.pos_rmin>region["max"]):
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
        print("amount of tabix regions: {}".format(df.shape[0]))
        print("amount of pruned regions: {}".format(len(regions)))
        tcall="tabix "+fname+" -h "
        for region in regions:
            tcall =tcall+" chr "+region["#chrom"]+":"+str(region["min"])+"-"+str(region["max"])
        call=shlex.split(tcall)
        #execute tabix call
        with open("temp.out","w") as out: 
            tbx=subprocess.run(call,stdout=out)
        tabixdf=None
        for t_df in pd.read_csv("temp.out",sep="\t",dtype={"#chrom":str,
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
            },chunksize=c_size,engine="c"):
            t_df=t_df[t_df["pval"]<args.sig_treshold_2]
            if not t_df.empty:
                if type(tabixdf) == type(None):
                    tabixdf=t_df.copy()
                else:
                    tabixdf=tabixdf.append(t_df)
        #tabixdf=tabixdf[tabixdf["pval"]<args.sig_treshold_2]
        tabixdf=tabixdf.drop_duplicates(subset=["#chrom","pos"],keep="first")
        new_df=pd.DataFrame(columns=result_dframe.columns).drop(["pos_rmin","pos_rmax"],axis=1)
        #find groups
        i=1
        total=result_dframe.shape[0]
        while not result_dframe.empty:
            #find most significant SNP
            ms_snp=result_dframe.loc[result_dframe["pval"].idxmin(),:]
            #take the group from tabixdf
            rowidx=(tabixdf["pos"]<=ms_snp["pos_rmax"])&(tabixdf["pos"]>=ms_snp["pos_rmin"])
            tmp=tabixdf.loc[rowidx,:].copy()
            tmp.loc[:,"locid"]=ms_snp["locid"]
            new_df=pd.concat([new_df,tmp],ignore_index=True,axis=0,join='inner')
            #convergence: remove the indexes from result_dframe
            dropidx=(result_dframe["pos"]<=ms_snp["pos_rmax"])&(result_dframe["pos"]>=ms_snp["pos_rmin"])
            result_dframe=result_dframe.loc[~dropidx,:]
            if i%100==0:
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
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("-w","--locus-width",dest="loc_width",type=int,default=250000,help="location width to include for each SNP")
    parser.add_argument("-g", "--group", dest="grouping",action='store_true',help="Whether to include p-values that are within location width and have pval<sig_treshold_2")
    parser.add_argument("-s2","--alternate-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    args=parser.parse_args()
    fetch_gws(args)