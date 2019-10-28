#!/usr/bin/env python3

import argparse,shlex,subprocess, glob
from subprocess import Popen, PIPE
import sys,os,io
import pandas as pd, numpy as np
import tabix
from autoreporting_utils import *

def parse_region(region):
    chrom=region.split(":")[0]
    start=region.split(":")[1].split("-")[0]
    end=region.split(":")[1].split("-")[1]
    return {"chrom":str(chrom),"start":int(start),"end":int(end)}


def solve_groups(all_variants,group_data):
    """
    Returns a dataframe containing the grouped variants, and only those.
    In: All variants considered for grouping (basically the variants with p-value < args.sig_tresh_2), parsed group data
    Out: Dataframe with same columns as all_variants, containing the rows in group_data
    """
    retval=pd.DataFrame(columns=all_variants.columns)
    for t in group_data.itertuples():
        group=[t.SNP]
        if t.TOTAL>0:
            sp2_split=t.SP2.split(",")
            strip_= lambda x,y: x[:-len(y)] if x.endswith(y) else x
            sp2_split=[strip_(x.strip(),"(1)" ) for x in sp2_split]
            group=group+sp2_split
        group_df=all_variants[ all_variants["#variant"].isin(group) ].copy()
        group_df.loc[:,"locus_id"]=t.SNP
        retval=pd.concat([retval,group_df],axis="index",sort=True)
    return retval

def simple_grouping(df_p1,df_p2,r,overlap,columns):
    """
    Simple grouping function
    Groups the variants based on p-values
    Group 1: lead variants
    Group 2: lead variants + variants that can not be lead variants. group 2 contains group 1.
    In: group 1 (df_p1, df), group 2 (df_p2, df), group width (r, int), do we overlap (overlap, bool), columns (columns, dict{key:str})
    Out: grouped dataframe
    """
    #simple grouping
    df=df_p1.copy()
    group_df=df_p2.copy()
    new_df=pd.DataFrame(columns=df.columns)
    while not df.empty:
        ms_snp=df.loc[df[columns["pval"] ].idxmin(),:]#get most sig SNP
        rowidx=(group_df[ columns["pos"] ]<=(ms_snp[columns["pos"]]+r) )&(group_df[ columns["pos"] ]>= (ms_snp[columns["pos"]]-r) )&(group_df[ columns["chrom"] ]==ms_snp[columns["chrom"]])#get group indexes
        tmp=group_df.loc[rowidx,:].copy()
        tmp.loc[:,"locus_id"]=ms_snp["#variant"]#necessary
        tmp.loc[:,"pos_rmin"]=tmp[columns["pos"]].min()#get group pos_rmin and pos_rmax from location of the SNPs instead of just a 
        tmp.loc[:,"pos_rmax"]=tmp[columns["pos"]].max()
        new_df=pd.concat([new_df,tmp],ignore_index=True,axis=0,join='inner')
        #convergence: remove the indexes from result_dframe. The rowidx =\= dropidx, as they come from different dataframes (df, group_df)
        dropidx=(df[ columns["pos"] ]<=(ms_snp[columns["pos"]]+r) )&(df[ columns["pos"] ]>= (ms_snp[columns["pos"]]-r) )&(df[ columns["chrom"] ]==ms_snp[columns["chrom"]])
        df=df.loc[~dropidx,:]
        #if not overlap, we also delete these from the group_df, so they can't be included again in other groups.
        if not overlap:
            t_dropidx=(group_df[ columns["pos"] ]<=(ms_snp[columns["pos"]]+r) )&(group_df[ columns["pos"] ]>=(ms_snp[columns["pos"]]-r) )&(group_df[ columns["chrom"] ]==ms_snp[columns["chrom"]])
            group_df=group_df.loc[~t_dropidx,:]
    return new_df

def load_credsets(fname,columns):
    """
    Load SuSiE credible sets from one bgzipped file.
    In: filename of the SuSie credible set file
    Out: A Dataframe containing all credible sets for this phenotype.
    """
    input_data = pd.read_csv(fname,sep="\t",compression="gzip")
    if input_data.empty:
        return pd.DataFrame(columns = columns.values()+["cs_prob","cs_id"])
    input_data = input_data[input_data["cs"]!=-1]#filter to credible sets
    input_data["credsetid"]=input_data[["region","cs"]].apply(lambda x: "".join([str(y) for y in x]),axis=1)
    data=input_data.rename(columns={"chromosome":columns["chrom"], "position": columns["pos"], "allele1": columns["ref"], "allele2": columns["alt"], "prob": "cs_prob"}).copy()
    data[columns["chrom"]] = data[columns["chrom"]].str.strip("chr")
    data["cs_id"]=np.NaN
    for i in data["credsetid"].unique():
        cs=data[data["credsetid"]==i]
        rsid=cs.loc[cs["cs_prob"].idxmax(),"rsid"]
        idx=cs.loc[cs["cs_prob"].idxmax(),"cs"]
        data.loc[data["credsetid"]==i,"cs_id"] = "{}_{}".format(rsid,idx)
    cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"], "cs_prob", "cs_id" ]
    data=data.loc[:,cols]
    return data

def ld_grouping(df_p1,df_p2, sig_treshold , sig_treshold_2, locus_width, ld_treshold,ld_panel_path,plink_memory,overlap, prefix, columns):
    """
    LD Clumping function
    Groups the variants based on PLINK's ld-clumping
    In: variant group 1, variant group 2, p-value threshold 1, p-value threshold 2, group width, ld threshold, prefix,  columns
    Out: grouped dataframe
    """
    #1:create PLINK variant list
    temp_variants="{}clump_variants.csv".format(prefix)
    df_p2.loc[:,["#variant",columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"] ]].to_csv(path_or_buf=temp_variants,index=False,sep="\t")
    plink_fname="{}plink_clump".format(prefix)
    #set up overlap flag for PLINK
    allow_overlap=""
    if overlap==True:
        allow_overlap="--clump-allow-overlap"
    plink_command="plink --allow-extra-chr --bfile {} --clump {} --clump-field {} --clump-snp-field '{}'  --clump-r2 {}"\
        " --clump-kb {} --clump-p1 {} --clump-p2 {} --out {} --memory {} {}".format(
        ld_panel_path,
        temp_variants,
        columns["pval"],
        "#variant",
        ld_treshold,
        locus_width,
        sig_treshold,
        sig_treshold_2,
        plink_fname,
        plink_memory,
        allow_overlap)
    #run PLINK
    pr = subprocess.Popen(shlex.split(plink_command), stdout=PIPE,stderr=subprocess.STDOUT,encoding='ASCII' )
    pr.wait()
    #get plink log
    plink_log=pr.stdout.readlines()
    if pr.returncode != 0:
        print("PLINK FAILURE. Error code {}".format(pr.returncode)  )
        [print(l) for l in plink_log]
        raise ValueError("Plink clumping returned code {}".format(pr.returncode))
    #Check if PLINK returned something or not
    no_sig_res_string="Warning: No significant --clump results.  Skipping."
    #if there is data
    if os.path.exists("{}.clumped".format(plink_fname)):
        group_data=pd.read_csv("{}.clumped".format(plink_fname),sep="\s+")
        group_data=group_data.loc[:,["SNP","TOTAL","SP2"]]
        new_df=solve_groups(df_p2.copy(),group_data)
        for var in new_df["locus_id"].unique():
            new_df.loc[new_df["locus_id"]==var,"pos_rmin"]=new_df.loc[new_df["locus_id"]==var,"pos"].min()#r["min"]
            new_df.loc[new_df["locus_id"]==var,"pos_rmax"]=new_df.loc[new_df["locus_id"]==var,"pos"].max()
        p1_group_leads = df_p1["#variant"].isin(new_df["#variant"])
        p1_singletons = ~p1_group_leads
        new_df=pd.concat([new_df,df_p1.loc[p1_singletons,:]],axis="index",sort=True).sort_values(by=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],"#variant"])
        new_df.loc[:,"pos_rmin"]=new_df.loc[:,"pos_rmin"].astype(np.int32)
        new_df.loc[:,"pos_rmax"]=new_df.loc[:,"pos_rmax"].astype(np.int32)
    #if there is not data
    else:
        if any([no_sig_res_string in string for string in plink_log]):#no significant results
            print("No significant results with PLINK clumping. All groups are singletons.")
            new_df=df_p1.copy()
        else:
            print("Plink .clumped file not found. Check the logs for information:")
            [print(l) for l in plink_log]
            raise FileNotFoundError("Plink .clumped file not found.")
    #cleanup plink files
    plink_files=glob.glob("{}.*".format(plink_fname))
    subprocess.call(["rm",temp_variants]+plink_files,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    return new_df

def credible_set_grouping(data,alt_sign_treshold,ld_panel_path,ld_treshold, locus_range,overlap,columns,prefix=""):
    """
    Create groups using credible set most probable variants as the lead variants, and rest of the data as the additional variants
    In: df(containing the credible set information), alternate significance threshold, ld panel path, columns
    Out: grouped df
    """
    df=data.copy()
    ld_window=1000
    lead_vars=[]
    #determine group leads. Group leads are 'the variants with largest cs_prob in that credible set'.
    for credible_set in df.loc[~df["cs_id"].isna(),"cs_id"].unique():
        group_lead =  df.loc[ df.loc[df["cs_id"]==credible_set,"cs_prob"].idxmax(),:]
        lead_vars.append(group_lead["#variant"])
    if len(lead_vars) == 0:
        return pd.DataFrame(columns=df.columns)
    #write lead_vars to a file, one per row
    fname="{}plink_ld.variants".format(prefix)
    output="{}plink_ld".format(prefix)
    with open(fname,"w") as f:
        for var in lead_vars:
            f.write("{}\n".format(var) )
    #perform plink computation
    plink_cmd = "plink --allow-extra-chr --bfile {} --r2 --ld-snp-list {} --ld-window-r2 {} --ld-window-kb {} --ld-window {} --out {}".format(ld_panel_path,
        fname,
        ld_treshold,
        locus_range,
        ld_window,
        output)
    pr = subprocess.Popen(shlex.split(plink_cmd),stdout=PIPE,stderr=subprocess.STDOUT,encoding='ASCII')
    pr.wait()
    plink_log = pr.stdout.readlines()
    if pr.returncode != 0:
        print("PLINK FAILURE. Error code {}".format(pr.returncode)  )
        [print(l) for l in plink_log]
        raise ValueError("Plink r2 calculation returned code {}".format(pr.returncode))
    #read in the variants
    ld_data=pd.read_csv("{}.ld".format(output),sep="\s+")
    #join
    df["index"]=df.index
    ld_df = pd.merge(df[["#variant",columns["chrom"],columns["pos"],columns["pval"],"index"]],ld_data, how="inner",left_on="#variant",right_on="SNP_B") #does include all of the lead variants as well
    ld_df=ld_df.set_index("index")
    ld_df.index.name=None
    #filter by p-value
    ld_df = ld_df[ld_df[columns["pval"]] <= alt_sign_treshold ]
    #Now should have columns variant, chrom, pos, pval, CHR_A,BP_A,SNP_A,CHR_B,BP_B, SNP_B, R2, with index being the same as in df.
    out_df = pd.DataFrame(columns=data.columns)
    #create df with only lead variants
    leads = df[df["#variant"].isin(lead_vars)].loc[:,["#variant",columns["pval"]]].copy()
    while not leads.empty:
        #all of the variants are in sufficient LD with the lead variants if there is a row where SNP_A is variant, and SNP_B (or #variant) is the variant in LD with the lead variant.
        lead_variant = leads.loc[leads[columns["pval"]].idxmin(),"#variant"]
        group_idx = ld_df[ld_df["SNP_A"] == lead_variant].index
        group_idx = group_idx[group_idx.isin(df.index)]#make sure that the index is present in both dataframes
        group = df.loc[group_idx,:].copy()
        #also add the variants that are in the credible set. Just in case they might not be included in the 
        credible_id = data.loc[data["#variant"]==lead_variant,"cs_id"].values[0]
        credible_set= data.loc[data["cs_id"] == credible_id,:].copy()
        #concat the two, remove duplicate entries
        group=pd.concat([group,credible_set],ignore_index=True,sort=False).drop_duplicates(subset=["#variant","cs_id"])
        group["locus_id"]=lead_variant
        group["pos_rmin"]=group["pos"].min()
        group["pos_rmax"]=group["pos"].max()
        out_df=pd.concat([out_df,group],ignore_index=True,axis=0,join='inner')
        #convergence: remove lead_variant, remove group from df in case overlap is not done
        leads=leads[~ (leads["#variant"] == lead_variant) ] 
        if not overlap:
            df=df[~df.index.isin(group_idx)]
            df=df[~(df["cs_id"]==credible_id)]
    #cleanup: delete variant file, plink files
    cleanup_cmd= "rm {}".format(fname)
    plink_files = glob.glob( "{}.*".format(output) )
    subprocess.call(shlex.split(cleanup_cmd)+plink_files, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return out_df

def get_gws_variants(fname, sign_treshold=5e-8,dtype=None,columns={"chrom":"#chrom","pos":"pos","ref":"ref","alt":"alt","pval":"pval"},compression="gzip"):
    """
    Get genome-wide significant variants from a summary statistic file.
    In: filename, significance threshold, dtype,columns,compression
    Out: dataframe containing the significant variants. No additional columns will be added.
    """
    chunksize=100000
    if not dtype:
        dtype={columns["chrom"]:str,
                columns["pos"]:np.int32,
                columns["ref"]:str,
                columns["alt"]:str,
                columns["pval"]:np.float64}
    retval=pd.DataFrame()
    for df in pd.read_csv(fname,compression=compression,sep="\t",dtype=dtype,engine="c",chunksize=chunksize):
        retval=pd.concat( [retval,df.loc[df[columns["pval"] ] <=sign_treshold,: ] ], axis="index", ignore_index=True )
    return retval

def merge_credset(gws_df,cs_df,fname,columns):
    """
    Merge credible set to the genome-wide significant variants. 
    In case variants in the credible set are not included in the gws variants, 
    the rows corresponding to them are fetched using tabix.
    In: Dataframe containing gws variants, dataframe containing credible sets, filename for summary statistic.
    Out: Dataframe containing the gws variants + any credible set variants that are not gws. Columns 'cs_id','cs_prob' added to the dataframe. 
    """
    join_cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"]]
    # fetch rows using tabix
    ## make tabix file constructor
    tb=tabix.open(fname)
    cred_row_df = load_tb_df(cs_df,tb,fname,columns=columns)
    # ensure only the credible sets were included
    cred_row_df = pd.merge(cred_row_df,cs_df[join_cols],how="right",on=join_cols)
    # concat to gws_df, remove duplicate rows using drop_duplicates
    df=gws_df.copy()
    df = pd.concat( [gws_df,cred_row_df], axis="index", ignore_index=True, sort=False).drop_duplicates(subset=list( join_cols ) )
    # merge the credible set
    merged = pd.merge(df,cs_df,how="left",on=join_cols)
    return merged

def fetch_gws(args):
    #column names
    columns={"chrom":args.column_labels[0],"pos":args.column_labels[1],"ref":args.column_labels[2],"alt":args.column_labels[3],"pval":args.column_labels[4]}
    sig_tresh=max(args.sig_treshold,args.sig_treshold_2)
    r=args.loc_width*1000#range for location width, originally in kb
    dtype={columns["chrom"]:str,
                columns["pos"]:np.int32,
                columns["ref"]:str,
                columns["alt"]:str,
                columns["pval"]:np.float64}

    #data input: get genome-wide significant variants.
    temp_df=get_gws_variants(args.gws_fpath,sign_treshold=sig_tresh,dtype=dtype,columns=columns,compression="gzip")
    #remove ignored region if there is one
    if args.ignore_region:
        ignore_region=parse_region(args.ignore_region)
        ign_idx=( ( temp_df[columns["chrom"]]==ignore_region["chrom"] ) & ( temp_df[columns["pos"]]<=ignore_region["end"] )&( temp_df[columns["pos"]]>=ignore_region["start"] ) )
        temp_df=temp_df.loc[~ign_idx,:]
    
    if temp_df.empty:
        print("The input file {} contains no gws-significant hits with signifigance treshold of {}. Aborting.".format(args.gws_fpath,args.sig_treshold))
        return 1
    
    #data input: get credible set variants
    join_cols=[columns["chrom"], columns["pos"], columns["ref"], columns["alt"]]
    if args.cred_set_file != "":
        cs_df=load_credsets(args.cred_set_file,columns)
    else:
        cs_df=pd.DataFrame(columns=join_cols+["cs_prob","cs_id"])
    #merge with gws_df, by using chrom,pos,ref,alt
    #temp_df = pd.merge(temp_df,cs_df,how="left",on=join_cols)
    temp_df = merge_credset(temp_df,cs_df,args.gws_fpath,columns)

    #create necessary columns for the data
    temp_df=temp_df.reset_index(drop=True)
    temp_df.loc[:,"#variant"]=create_variant_column(temp_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
    temp_df.loc[:,"locus_id"]=temp_df.loc[:,"#variant"]
    temp_df.loc[:,"pos_rmax"]=temp_df.loc[:,columns["pos"]]
    temp_df.loc[:,"pos_rmin"]=temp_df.loc[:,columns["pos"]]
    df_p1=temp_df.loc[temp_df[columns["pval"]] <= args.sig_treshold,: ].copy()
    df_p2=temp_df.loc[temp_df[columns["pval"]] <= args.sig_treshold_2,: ].copy()
    #grouping
    if args.grouping and not df_p1.empty:
        if args.grouping_method=="ld":
            new_df=ld_grouping(df_p1,df_p2,args.sig_treshold,args.sig_treshold_2,args.loc_width,args.ld_r2,args.ld_panel_path,args.plink_mem,args.overlap,args.prefix,columns)
        elif args.grouping_method=="cred":
            new_df = credible_set_grouping(df_p2,args.sig_treshold_2,args.ld_panel_path,args.ld_r2,args.loc_width,args.overlap,columns,args.prefix)
        else :
            new_df=simple_grouping(df_p1=df_p1,df_p2=df_p2,r=r,overlap=args.overlap,columns=columns)
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
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, options ['ld','simple','cred']")
    parser.add_argument("--locus-width-kb",dest="loc_width",type=int,default=250,help="locus width to include for each SNP, in kb")
    parser.add_argument("--alt-sign-treshold",dest="sig_treshold_2",type=float, default=5e-8,help="optional group treshold")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--ld-r2", dest="ld_r2", type=float, default=0.4, help="r2 cutoff for ld clumping")
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    parser.add_argument("--overlap",dest="overlap",action="store_true",help="Are groups allowed to overlap")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval"],help="Names for data file columns. Default is '#chrom pos ref alt pval'.")
    parser.add_argument("--ignore-region",dest="ignore_region",type=str,default="",help="Ignore the given region, e.g. HLA region, from analysis. Give in CHROM:BPSTART-BPEND format.")
    parser.add_argument("--credible-set-file",dest="cred_set_file",type=str,default="",help="bgzipped SuSiE credible set file.")
    args=parser.parse_args()
    if args.prefix!="":
        args.prefix=args.prefix+"."
    args.fetch_out = "{}{}".format(args.prefix,args.fetch_out)
    fetch_gws(args)