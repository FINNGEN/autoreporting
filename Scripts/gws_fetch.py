#! /usr/bin/python3

import argparse,shlex,subprocess
from subprocess import Popen, PIPE
import pandas as pd 
import numpy as np
import tabix

def tabix_command(df,path,min_="pos_rmin",max_="pos_rmax",prefix=""):
    tcall="tabix "+path+" -h " 
    for idx,r in df.iterrows():
        #add correct row
        tcall=tcall+str(prefix)+" "+str(r["#chrom"])+":"+str(r[min_])+"-"+str(r[max_])+" "
    return tcall

def pytabix(tb,chrom,start,end):
    return tb.querys("{}:{}-{}".format(chrom,start,end))

def create_variant_column(df,chrom="#chrom",pos="pos",ref="ref",alt="alt"):
    return df.apply( lambda x: "chr{}_{}_{}_{}".format(x[chrom],x[pos],x[ref],x[alt]) ,axis=1)

def prune_regions(df):
    regions=[]
    for t in df.itertuples():
        if regions:
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
    return pd.DataFrame(regions)

def get_gzip_header(fname):
    """"Returns header for gzipped tsvs, as that is not currently possible using pytabix
    In: file path of gzipped tsv
    Out: header of tsv as a list of column names"""
    #create gzip call
    gzip_call=shlex.split("gzip -cd {}".format(fname))
    head_call=shlex.split("head -n 1")
    gzip_process=Popen(gzip_call,stdout=PIPE)
    head_proc=Popen(head_call,stdin=gzip_process.stdout,stdout=PIPE)
    out=[]
    for line in head_proc.stdout:
        out.append(line.decode().strip().split("\t"))
    return out[0]

def parse_plink_output(df):
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
    res={"#chrom":[],"pos":[],"ref":[],"alt":[]}
    for row in chrom_pos_lst:
        tmp=row.strip("chr").split("_")
        res["#chrom"].append(tmp[0])
        res["pos"].append(tmp[1])
        res["ref"].append(tmp[2])
        res["alt"].append(tmp[3])
    res=pd.DataFrame(res)
    res.loc[:,"variant"]=create_variant_column(res)
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

def fetch_gws(args):
    fname=args.gws_fpath
    sig_tresh=args.sig_treshold
    c_size=100000
    r=args.loc_width*1000#range for location width, originally in kb
    df=None
    dtype={"#chrom":str,
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
                }
    for dframe in pd.read_csv(fname,compression="gzip",sep="\t",dtype=dtype,engine="c",chunksize=c_size):
        #filter gws snps
        temp_dframe=dframe.loc[dframe["pval"]<sig_tresh,:]
        if not temp_dframe.empty:
            if type(df) == type(None):
                df=temp_dframe.copy()
            else:
                df=df.append(temp_dframe)
    df=df.reset_index(drop=True)
    df.loc[:,"#variant"]=create_variant_column(df)
    df.loc[:,"locus_id"]=df.loc[:,"#variant"]
    result_dframe=df
    new_df=None

    if args.grouping:
        if args.grouping_method=="ld":
            #write current SNPs to file
            temp_variants="ld_variants.csv"
            df.loc[:,["#variant","#chrom","pos","ref","alt","pval"]].to_csv(path_or_buf=temp_variants,index=False,sep="\t")
            plink_fname="temp_plink"
            plink_command="plink --allow-extra-chr --bfile {} --clump {} --clump-field {} --clump-snp-field '{}'  --clump-r2 {}"\
                " --clump-kb {} --clump-p1 {} --clump-p2 {} --out {} --memory {} --clump-allow-overlap".format(
                args.ld_panel_path,
                temp_variants,
                "pval",
                "#variant",
                args.ld_r2,
                args.loc_width,
                args.sig_treshold,
                args.sig_treshold_2,
                plink_fname,
                args.plink_mem)
            #call plink
            subprocess.call(shlex.split(plink_command), stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL )
            #parse output file, find locus width
            group_data=pd.read_csv(plink_fname+".clumped",sep="\s+")
            group_data=group_data.loc[:,["SNP","TOTAL","SP2"]]
            res=parse_plink_output(group_data)
            #fetch tabix
            tbxlst=[]
            tb=tabix.open(fname)
            for _,row in res.iterrows():
                tbxlst=tbxlst+list(pytabix(tb,row["#chrom"],int(row["pos"]),int(row["pos"]) ) )
            tbxheader=get_gzip_header(fname)
            #transform into dataframe
            tabixdf=pd.DataFrame(tbxlst,columns=tbxheader )
            tabixdf=tabixdf.astype(dtype=dtype)
            #filter, create variant column
            tabixdf=tabixdf.drop_duplicates(subset=["#chrom","pos","ref","alt"],keep="first")
            tabixdf=tabixdf[tabixdf["pval"]<args.sig_treshold_2]
            tabixdf.loc[:,"#variant"]=create_variant_column(tabixdf)
            #add tabix info to groups and write to file
            df=pd.DataFrame(columns=tbxheader+["#variant","locus_id"])
            df=solve_groups(df,group_data,tabixdf)
            df.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)

        else:
            #simple grouping
            df.loc[:,"pos_rmin"]=df.loc[:,"pos"]-r
            df.loc[:,"pos_rmax"]=df.loc[:,"pos"]+r
            df.loc[:,"pos_rmin"]=df.loc[:,"pos_rmin"].clip(lower=0)
            reg_df=prune_regions(df.loc[:,["#chrom", "pos_rmin","pos_rmax"]])
            print("amount of tabix regions: {}".format(df.shape[0]))
            print("amount of pruned regions: {}".format(reg_df.shape[0]))
            tbxlst=[]
            tb=tabix.open(fname)
            for _,row in reg_df.iterrows():
                tbxlst=tbxlst+list(pytabix(tb,row["#chrom"],int(row["min"]),int(row["max"]) ))
            tbxheader=get_gzip_header(fname)
            tabixdf=pd.DataFrame(tbxlst,columns=tbxheader)
            tabixdf=tabixdf.astype(dtype=dtype)

            tabixdf=tabixdf[tabixdf["pval"]<args.sig_treshold_2]
            tabixdf=tabixdf.drop_duplicates(subset=["#chrom","pos","ref","alt"],keep="first")
            new_df=pd.DataFrame(columns=df.columns).drop(["pos_rmin","pos_rmax"],axis=1)
            tabixdf.to_csv("tbx_debug.csv",sep="\t",index=False)
            i=1
            total=df.shape[0]
            while not df.empty:
                ms_snp=df.loc[df["pval"].idxmin(),:]
                rowidx=(tabixdf["pos"]<=ms_snp["pos_rmax"])&(tabixdf["pos"]>=ms_snp["pos_rmin"])
                tmp=tabixdf.loc[rowidx,:].copy()
                tmp.loc[:,"locus_id"]=ms_snp["#variant"]
                tmp.loc[:,"#variant"]=create_variant_column(tmp)
                new_df=pd.concat([new_df,tmp],ignore_index=True,axis=0,join='inner')
                #convergence: remove the indexes from result_dframe
                dropidx=(df["pos"]<=ms_snp["pos_rmax"])&(df["pos"]>=ms_snp["pos_rmin"])
                df=df.loc[~dropidx,:]
                if i%100==0:
                    print("iter: {}, SNPs remaining:{}/{}".format(i,result_dframe.shape[0],total))
                i+=1
            new_df.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
    else:
        df.loc[:,["#chrom","pos","ref","alt","rsids",
                    "nearest_genes","pval","beta","sebeta","maf",
                    "maf_cases","maf_controls","#variant","locus_id"]].to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
    
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Fetch and group gws SNPs from summary statistic")
    parser.add_argument("gws_fpath",type=str,help="Filepath of the compressed tsv")
    parser.add_argument("-s","--signifigance-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("-o","--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    parser.add_argument("-g", "--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("-w","--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("-s2","--alternate-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    
    args=parser.parse_args()
    fetch_gws(args)