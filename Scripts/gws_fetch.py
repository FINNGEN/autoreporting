#!/usr/bin/env python3

import argparse,shlex,subprocess, glob
from subprocess import Popen, PIPE
import sys,os,io
import pandas as pd, numpy as np
#import tabix
from autoreporting_utils import *

def parse_region(region):
    chrom=region.split(":")[0]
    start=region.split(":")[1].split("-")[0]
    end=region.split(":")[1].split("-")[1]
    return {"chrom":str(chrom),"start":int(start),"end":int(end)}

def parse_plink_output(df,columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt"}):
    """Parse plink --clump output, which is a multiple-space separated dataframe with the groups
    separated by commas. Is used to get tabixdf
    In: plink --clump output dataframe
    Out: dataframe that has chrom,pos,ref,alt of all of the grouped and group variants"""
    chrom_pos_lst=[]
    for t in df.itertuples():
        chrom_pos_lst.append(t.SNP)
        if t.TOTAL != 0:
            sp2_split=t.SP2.split(",")
            sp2_split=[x.strip().strip("(1)") for x in sp2_split]
            chrom_pos_lst=chrom_pos_lst+sp2_split
    #separate to chrom and pos
    res={columns["chrom"]:[],columns["pos"]:[],columns["ref"]:[],columns["alt"]:[]}
    for row in chrom_pos_lst:
        tmp=row.strip("chr").split("_")
        res[columns["chrom"] ].append(tmp[0])
        res[columns["pos"] ].append(tmp[1])
        res[columns["ref"] ].append(tmp[2])
        res[columns["alt"] ].append(tmp[3])
    res=pd.DataFrame(res)
    res.loc[:,"#variant"]=create_variant_column(res,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
    return res

def solve_groups(result_dframe,group_data,tabixdf):
    for t in group_data.itertuples():
        group=[t.SNP]
        if t.TOTAL>0:
            sp2_split=t.SP2.split(",")
            sp2_split=[x.strip().strip("(1)") for x in sp2_split]
            group=group+sp2_split
        tmp=tabixdf.loc[tabixdf["#variant"].isin(group),].copy(deep=True)
        tmp.loc[:,"locus_id"]=t.SNP
        result_dframe=pd.concat([result_dframe,tmp],axis=0)
    return result_dframe

def get_group_range(dframe,group_variant,columns={"pos":"pos"}):
    #find min and max position from group, return them
    temp_df=dframe.loc[dframe["locus_id"]==group_variant]
    if temp_df.empty:
        return None
    min_=np.min(temp_df[columns["pos"]])
    max_=np.max(temp_df[columns["pos"]])
    return {"min":min_,"max":max_}

def fetch_gws(args):
    #column names
    columns={"chrom":args.column_labels[0],"pos":args.column_labels[1],"ref":args.column_labels[2],"alt":args.column_labels[3],"pval":args.column_labels[4]}
    sig_tresh=max(args.sig_treshold,args.sig_treshold_2)
    c_size=100000
    r=args.loc_width*1000#range for location width, originally in kb
    dtype={columns["chrom"]:str,
                columns["pos"]:np.int32,
                columns["ref"]:str,
                columns["alt"]:str,
                columns["pval"]:np.float64}
    temp_df=pd.DataFrame()
    for dframe in pd.read_csv(args.gws_fpath,compression="gzip",sep="\t",dtype=dtype,engine="c",chunksize=c_size):
        #filter gws snps
        temp_df=pd.concat([temp_df,dframe.loc[dframe[columns["pval"]]<=sig_tresh,:]],axis="index",ignore_index=True)
    if temp_df.empty:
        print("The input file {} contains no gws-significant hits with signifigance treshold of {}. Aborting.".format(args.gws_fpath,args.sig_treshold))
        return 1
    #remove ignored region if there is one
    if args.ignore_region:
        ignore_region=parse_region(args.ignore_region)
        ign_idx=( ( temp_df[columns["chrom"]]==ignore_region["chrom"] ) & ( temp_df[columns["pos"]]<=ignore_region["end"] )&( temp_df[columns["pos"]]>=ignore_region["start"] ) )
        temp_df=temp_df.loc[~ign_idx,:]
    temp_df=temp_df.reset_index(drop=True)
    temp_df.loc[:,"#variant"]=create_variant_column(temp_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
    temp_df.loc[:,"locus_id"]=temp_df.loc[:,"#variant"]
    temp_df.loc[:,"pos_rmax"]=temp_df.loc[:,columns["pos"]]
    temp_df.loc[:,"pos_rmin"]=temp_df.loc[:,columns["pos"]]
    df_p1=temp_df.loc[temp_df[columns["pval"]] <= args.sig_treshold,: ].copy()
    df_p2=temp_df.loc[temp_df[columns["pval"]] <= args.sig_treshold_2,: ].copy()
    if args.grouping and not df_p1.empty:
        if args.grouping_method=="ld":
            #write current SNPs to file
            temp_variants="{}clump_variants.csv".format(args.prefix)
            df_p2.loc[:,["#variant",columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"] ]].to_csv(path_or_buf=temp_variants,index=False,sep="\t")
            plink_fname="{}plink_clump".format(args.prefix)
            allow_overlap=""
            if args.overlap==True:
                allow_overlap="--clump-allow-overlap"
            plink_command="plink --allow-extra-chr --bfile {} --clump {} --clump-field {} --clump-snp-field '{}'  --clump-r2 {}"\
                " --clump-kb {} --clump-p1 {} --clump-p2 {} --out {} --memory {} {}".format(
                args.ld_panel_path,
                temp_variants,
                columns["pval"],
                "#variant",
                args.ld_r2,
                args.loc_width,
                args.sig_treshold,
                args.sig_treshold_2,
                plink_fname,
                args.plink_mem,
                allow_overlap)
            #call plink
            pr = subprocess.Popen(shlex.split(plink_command), stdout=PIPE,stderr=subprocess.STDOUT,encoding='ASCII' )
            pr.wait()
            plink_log=pr.stdout.readlines()
            if pr.returncode!=0:
                print("PLINK FAILURE. Error code {}".format(pr.returncode)  )
                [print(l) for l in plink_log]
                raise ValueError("Plink clumping returned code {}".format(pr.returncode))
            #parse output file, find locus width
            with open("{}plink_log.log".format(args.prefix),"w") as f:
                f.writelines(plink_log)
            try:
                group_data=pd.read_csv("{}.clumped".format(plink_fname),sep="\s+")
            except:
                print("Plink .clumped file not found. Plink logs:")
                [print(l) for l in plink_log]
                raise 
            group_data=group_data.loc[:,["SNP","TOTAL","SP2"]]
            res=parse_plink_output(group_data,columns=columns)
            new_df=pd.DataFrame(columns=df_p2.columns)
            new_df=solve_groups(new_df,group_data,df_p2)
            for var in new_df["locus_id"].unique():
                r=get_group_range(new_df,var,columns=columns)
                new_df.loc[new_df["locus_id"]==var,"pos_rmin"]=r["min"]
                new_df.loc[new_df["locus_id"]==var,"pos_rmax"]=r["max"]
            new_df.loc[:,"pos_rmin"]=new_df.loc[:,"pos_rmin"].astype(np.int32)
            new_df.loc[:,"pos_rmax"]=new_df.loc[:,"pos_rmax"].astype(np.int32)
            #cleanup plink files
            plink_files=glob.glob("{}.*".format(plink_fname))
            subprocess.call(["rm",temp_variants]+plink_files,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        else:
            #simple grouping
            df=df_p1.copy()
            df.loc[:,"pos_rmin"]+= -r
            df.loc[:,"pos_rmax"]+= r
            df.loc[:,"pos_rmin"]=df.loc[:,"pos_rmin"].clip(lower=0)
            group_df=df_p2.copy()
            new_df=pd.DataFrame(columns=df.columns)
            i=1
            total=df.shape[0]
            while not df.empty:
                ms_snp=df.loc[df[columns["pval"] ].idxmin(),:]
                rowidx=(group_df[ columns["pos"] ]<=ms_snp["pos_rmax"])&(group_df[ columns["pos"] ]>=ms_snp["pos_rmin"])&(group_df[ columns["chrom"] ]==ms_snp[columns["chrom"]])
                tmp=group_df.loc[rowidx,:].copy()
                tmp.loc[:,"locus_id"]=ms_snp["#variant"]
                tmp.loc[:,"pos_rmin"]=ms_snp["pos_rmin"]
                tmp.loc[:,"pos_rmax"]=ms_snp["pos_rmax"]
                new_df=pd.concat([new_df,tmp],ignore_index=True,axis=0,join='inner')
                #convergence: remove the indexes from result_dframe
                dropidx=(df[ columns["pos"] ]<=ms_snp["pos_rmax"])&(df[ columns["pos"] ]>=ms_snp["pos_rmin"])&(df[ columns["chrom"] ]==ms_snp[columns["chrom"]])
                df=df.loc[~dropidx,:]
                if not args.overlap:
                    t_dropidx=(group_df[ columns["pos"] ]<=ms_snp["pos_rmax"])&(group_df[ columns["pos"] ]>=ms_snp["pos_rmin"])&(group_df[ columns["chrom"] ]==ms_snp[columns["chrom"]])
                    group_df=group_df.loc[~t_dropidx,:]
                if i%100==0:
                    print("iter: {}, lead SNPs remaining:{}/{}".format(i,df.shape[0],total))
                i+=1
        new_df=new_df.sort_values(["locus_id","#variant"])
        new_df.to_csv(path_or_buf=args.fetch_out,sep="\t",index=False)
    else:
        #take only gws hits, no groups. Therefore, use df_p1
        df_p1.sort_values(["locus_id","#variant"]).to_csv(path_or_buf=args.fetch_out,sep="\t",index=False)
    return 0
    
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Fetch and group genome-wide significant variants from summary statistic")
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed summary statistic")
    parser.add_argument("--sign-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix. Default value is the base name (no path and no file extensions) of input file. ")
    parser.add_argument("--fetch-out",dest="fetch_out",type=str,default="fetch_out.csv",help="GWS output filename, default is fetch_out.csv")
    parser.add_argument("--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("--alt-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    parser.add_argument("--overlap",dest="overlap",action="store_true",help="Are groups allowed to overlap")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval"],help="Names for data file columns. Default is '#chrom pos ref alt pval'.")
    parser.add_argument("--ignore-region",dest="ignore_region",type=str,default="",help="Ignore the given region, e.g. HLA region, from analysis. Give in CHROM:BPSTART-BPEND format.")
    args=parser.parse_args()
    if args.prefix!="":
        args.prefix=args.prefix+"."
    args.fetch_out = "{}{}".format(args.prefix,args.fetch_out)
    fetch_gws(args)