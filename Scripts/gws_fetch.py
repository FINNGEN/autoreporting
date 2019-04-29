#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd 
import numpy as np

def tabix_command(df,path,min_="pos_rmin",max_="pos_rmax",prefix=""):
    tcall="tabix "+path+" -h " 
    for idx,r in df.iterrows():
        #add correct row
        tcall=tcall+str(prefix)+" "+str(r["#chrom"])+":"+str(r[min_])+"-"+str(r[max_])+" "
    return tcall

def fetch_gws(args):
    fname=args.fpath
    sig_tresh=args.sig_treshold

    c_size=100000
    r=args.loc_width*1000#range for location width, originally in kb
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
    df.loc[:,"#variant"]="chr"+df.loc[:,"#chrom"].map(str)+"_"+df.loc[:,"pos"].map(str)+"_"+df.loc[:,"ref"].map(str)+"_"+df.loc[:,"alt"].map(str)
    df.loc[:,"locus_id"]=df.loc[:,"#variant"]
    result_dframe=df
    new_df=None
    if args.grouping:
        if args.grouping_method=="ld":
            #write current SNPs to file
            
            df.loc[:,["#variant","#chrom","pos","ref","alt","pval"]].to_csv(path_or_buf="ld_variants.csv",index=False,sep="\t")
            #raise NotImplementedError("ld grouping not implemented yet")
            plink_fname="temp_plink"
            plink_command="plink --allow-extra-chr --bfile "+args.ld_panel_path+" --clump ld_variants.csv --clump-field pval --clump-snp-field '#variant' --clump-r2 "\
            +str(args.ld_r2)+" --clump-kb "+str(args.loc_width) +" --clump-p1 "+str(args.sig_treshold)+" --clump-p2 "+str(args.sig_treshold_2) +" --out "+ plink_fname+ " --memory 12000 --clump-allow-overlap"
            #call plink
            subprocess.call(shlex.split(plink_command), stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL )
            #parse output file, find locus width
            #The output file is organized into the groups, for which they have all the members in a separate column, as comma-separated, with (#) added.
            #So, it seems I'll have to take the rows, and unroll them so that the groups have a group id, and then the id can be added to those SNPs
            group_data=pd.read_csv(plink_fname+".clumped",sep="\s+")
            group_data=group_data.loc[:,["SNP","TOTAL","SP2"]]
            #fetch the correct SNPs from our data, maybe easiest to get using tabix
            #TODO: do this, first make a df with the correct positions and then fetch it with tabix
            chrom_pos_lst=[]
            for t in group_data.itertuples():
                chrom_pos_lst.append(t.SNP)
                if t.TOTAL != 0:
                    sp2_split=t.SP2.split(",")
                    sp2_split=[x.strip().strip("(1)") for x in sp2_split]
                    chrom_pos_lst=chrom_pos_lst+sp2_split
            #separate to chrom and pos
            res={}
            res["#chrom"]=[]
            res["pos"]=[]

            for row in chrom_pos_lst:
                tmp=row.strip("chr").split("_")
                res["#chrom"].append(tmp[0])
                res["pos"].append(tmp[1])
            #fetch tabix
            res=pd.DataFrame(res)
            tcall=tabix_command(res,fname,"pos","pos")
            call=shlex.split(tcall)
            with open("temp.out","w") as out: 
                tbx=subprocess.run(call,stdout=out)
            #add groups to data
            tabixdf=pd.read_csv("temp.out",sep="\t")
            tabixdf=tabixdf.drop_duplicates(subset=["#chrom","pos","ref","alt"],keep="first")
            tabixdf.loc[:,"#variant"]="chr"+tabixdf.loc[:,"#chrom"].map(str)+"_"+tabixdf.loc[:,"pos"].map(str)+"_"+tabixdf.loc[:,"ref"].map(str)+"_"+tabixdf.loc[:,"alt"].map(str)
            tabixdf.loc[:,"locus_id"]=tabixdf.loc[:,"#variant"]
            for t in group_data.itertuples():
                #split SP2 to different values
                if t.TOTAL>0:
                    sp2_split=t.SP2.split(",")
                    sp2_split=[x.strip().strip("(1)") for x in sp2_split]
                    #add the correct locus id
                    tabixdf.loc[tabixdf["#variant"].isin(sp2_split),"locus_id"]=t.SNP
            tabixdf.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
            subprocess.run(shlex.split("rm temp.out"))


        else:
            #simple grouping
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
            reg_df=pd.DataFrame(regions)
            tcall=tabix_command(reg_df,fname,"min","max")
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
            tabixdf=tabixdf.drop_duplicates(subset=["#chrom","pos"],keep="first")
            new_df=pd.DataFrame(columns=result_dframe.columns).drop(["pos_rmin","pos_rmax"],axis=1)
            i=1
            total=result_dframe.shape[0]
            while not result_dframe.empty:
                ms_snp=result_dframe.loc[result_dframe["pval"].idxmin(),:]
                rowidx=(tabixdf["pos"]<=ms_snp["pos_rmax"])&(tabixdf["pos"]>=ms_snp["pos_rmin"])
                tmp=tabixdf.loc[rowidx,:].copy()
                tmp.loc[:,"#variant"]=ms_snp["#variant"]
                new_df=pd.concat([new_df,tmp],ignore_index=True,axis=0,join='inner')
                #convergence: remove the indexes from result_dframe
                dropidx=(result_dframe["pos"]<=ms_snp["pos_rmax"])&(result_dframe["pos"]>=ms_snp["pos_rmin"])
                result_dframe=result_dframe.loc[~dropidx,:]
                if i%100==0:
                    print("iter: {}, SNPs remaining:{}/{}".format(i,result_dframe.shape[0],total))
                i+=1
            #subprocess.run(shlex.split("rm temp.out"))
            new_df.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
    else:
        result_dframe.loc[:,["#chrom","pos","ref","alt","rsids",
                    "nearest_genes","pval","beta","sebeta","maf",
                    "maf_cases","maf_controls","#variant","locus_id"]].to_csv(path_or_buf=args.out_fname,sep="\t",index=False)
    

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Fetching gws SNPs from phenotype data")
    parser.add_argument("fpath",type=str,help="Filepath of the compressed tsv")
    parser.add_argument("-s","--signifigance-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("-o","--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    parser.add_argument("-g", "--group", dest="grouping",action='store_true',help="Whether to group SNPs")
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple or ld, default simple")
    parser.add_argument("-w","--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("-s2","--alternate-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    
    args=parser.parse_args()
    fetch_gws(args)