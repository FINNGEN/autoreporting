#! /usr/bin/python3

import argparse,shlex,subprocess
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
from autoreporting_utils import *
import gwcatalog_api
import os

def map_alleles(a1,a2):
    """
    Flips alleles to the A strand if neccessary and orders them lexicogaphically
    Author: Pietro?
    """
    allele_dict={}
    for n1,n2 in [("T","A"),("C","G"),("G","C")]:
        allele_dict[n1]=n2
     # check if the A variant is present
    if 'A' not in a1 + a2 and 'a' not in a1+a2:
        # for both/ref and alt map them to the A strand and order each one lexicographically
        a1 = ''.join([allele_dict[elem.upper()] for elem in a1])
        a2 = ''.join([allele_dict[elem.upper()] for elem in a2])
    # further sorting
    return sorted([a1,a2])
    

def compare(args):
    
    #load original file
    df=pd.read_csv(args.compare_fname,sep="\t")
    
    #if using summary file
    if args.compare_style=="file":
        #load summary file
        summary_df=pd.DataFrame()
        for idx in range(0,len(args.summary_files) ):
            s_df=pd.read_csv(args.summary_files[idx],sep="\t")
            #make sure the variant column exists
            summary_cols=s_df.columns
            if "#variant" not in summary_cols:
                s_df.loc[:, "summary_#variant"]=create_variant_column(summary_df)
            
            s_df.loc[:,"trait"] = args.endpoints[idx]
            #s_df=s_df.rename(columns=summary_rename)
            if not args.build_38:
                raise NotImplementedError("Non-build 38 raports are not supported yet.")
            #TODO: add summary files to end of each other, hope that they have the correct columns
            summary_df=pd.concat([summary_df,s_df],axis=0)
                    
    elif args.compare_style=="gwascatalog":
        gwas_df=None
        if os.path.exists("gwas_out_mapping.csv"):
            gwas_df=pd.read_csv("gwas_out_mapping.csv",sep="\t")
        else:
            #create ranges that contain all of the SNPs
            range_df=df.loc[:,["#chrom","pos"]].copy(deep=True)
            range_df.loc[:,"pos2"]=range_df.loc[:,"pos"]
            range_df=range_df.rename(columns={"pos":"pos_rmin","pos2":"pos_rmax"})
            #add 100000 bp to the interval to make the ranges a bit more common
            pad=args.gwascatalog_pad*1000
            range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"]-pad
            range_df.loc[:,"pos_rmax"]=range_df.loc[:,"pos_rmax"]+pad
            range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"].clip(lower=0)
            regions=prune_regions(range_df)
            #use api to get all gwascatalog hits
            print(regions)
            result_lst=[]
            for _,region in regions.iterrows():
                val=gwcatalog_api.get_all_associations(chromosome=region["#chrom"],
                    bp_lower=region["min"],bp_upper=region["max"],p_upper=args.gwascatalog_pval)
                if val:
                    result_lst+=gwcatalog_api.parse_output(val)
            gwas_df=pd.DataFrame(result_lst)
            gwas_df["trait"]=gwas_df["trait"].apply(lambda x:x[0])
            gwas_df.to_csv("gwas_out_mapping.csv",sep="\t")
        #parse hits to a proper form, drop unnecessary information
        #print(gwas_df.columns)
        gwas_cols=["base_pair_location","chromosome","p_value","hm_effect_allele","hm_other_allele",
            "trait","study_accession","hm_code","hm_beta","hm_effect_allele_frequency"]
        gwas_df=gwas_df.loc[:,gwas_cols]
        gwas_rename={"base_pair_location":"pos","chromosome":"#chrom","hm_beta":"beta",
            "p_value":"pval","hm_code":"code","hm_effect_allele":"alt","hm_other_allele":"ref","hm_effect_allele_frequency":"af"}
        gwas_df=gwas_df.rename(columns=gwas_rename)
        
        #gwas_df.loc[:,"code_column"].update(gwas_df.loc[:,"hm_code"])
        tmp_df=gwas_df[["#chrom","pos","ref","alt","pval","code","beta","af","trait"]]
        #assuming code is the same as hm_code, we want to filter out 9, 14, 15, 16, 17, 18
        filter_out_codes=[9, 14, 15, 16, 17, 18]
        tmp_df=tmp_df.loc[~tmp_df.loc[:,"code"].isin(filter_out_codes)]
        tmp_df.loc[:,"#variant"]=create_variant_column(tmp_df,chrom="#chrom",pos="pos",ref="ref",alt="alt")
        #change alleles to a strand using the code from commons
        summary_df=tmp_df
        #create list of unique traits
        unique_efos=list(summary_df["trait"].unique())
        trait_name_map={}
        for key in unique_efos:
            trait_name_map[key]=gwcatalog_api.get_trait_name(key)
        
        summary_df.loc[:,"trait_name"]=summary_df.loc[:,"trait"].apply(lambda x: trait_name_map[x])
        
    else:
        raise NotImplementedError("comparison method '{}' not yet implemented".format(args.compare_style))
    #now we should have df and summary_df
    for _,row in summary_df.iterrows():
        [alt,ref]=map_alleles(row["alt"],row["ref"])
        summary_df.loc[_,"map_ref"]=ref
        summary_df.loc[_,"map_alt"]=alt

    for _,row in df.iterrows():
        [alt,ref]=map_alleles(row["alt"],row["ref"])
        df.loc[_,"map_ref"]=ref
        df.loc[_,"map_alt"]=alt
    
    summary_df.loc[:,"map_variant"]=create_variant_column(summary_df,chrom="#chrom",pos="pos",ref="map_ref",alt="map_alt")
    df.loc[:,"map_variant"]=create_variant_column(df,chrom="#chrom",pos="pos",ref="map_ref",alt="map_alt")
    tmp=pd.merge(df,summary_df,how="inner",on="map_variant")
    #TODO: make options for ld, i.e. if ld_check parameter is supplied, for each of our hits we check if there are any summary stat variants 
    #in ld with them
    if args.ld_check:
        #raise NotImplementedError("LD comparison not yet implemented".format(args.compare_style))
        #preprocess found variants, i.e. divide variants in both our and gwascatalog results to chromosomes
        #Then, we need to build the bim files, I guess, or I could just process some ld panels ready and push them into a bucket
        #Then, for each chromosome set, we want to separate the clear ranges from each other, I guess
        #then, we run ldstore, export correlations, and see if any of the gwascatalog variants are in ld with the gws variants
        #NOTE: ldstore only takes bim file into account if we use range.

        #get variant list, i.e. the list of variants that consists of gws results and summary reuslts
        var_cols=["#variant","pos","#chrom","ref","alt"]
        var_rename={"#variant":"RSID","pos":"position","#chrom":"chromosome","ref":"A_allele","alt":"B_allele"}
        var_lst_df=pd.concat([summary_df.loc[:,var_cols].rename(columns=var_rename),df.loc[:,var_cols].rename(columns=var_rename)  ])
        var_lst_df=var_lst_df.drop_duplicates(subset=["RSID"])
        #get unique chromosomes in the results
        unique_chrom_list=df["#chrom"].unique()
        c1="mkdir temp"
        Popen(shlex.split(c1),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        c2="ln -s ../{}.bed ./temp/temp.bed".format(args.ld_panel_path)
        c3="ln -s ../{}.fam ./temp/temp.fam".format(args.ld_panel_path)
        print(c2)
        print(c3)
        
        Popen(shlex.split(c2),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        Popen(shlex.split(c3),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        ld_df=pd.DataFrame()
        for chromosome in unique_chrom_list:
            pass
            #build bim file, containing both the variants in df and summary_df
            bim_lst=var_lst_df[var_lst_df["chromosome"]==chromosome]
            bim_lst["cm"]=0
            bim_lst=bim_lst[["chromosome","RSID","cm","position","A_allele","B_allele"]]
            #temporary solution: create temp folder, create symbolic links to other files, write bim file to temp folder,
            #do the magic, and delete temp folder afterwards. No modification to existing files, and not a lot of disk space usage
            #still ugly, but at least we don't modify data that is supplied to this script
            range_upper=np.max(bim_lst["position"])
            range_lower=np.min(bim_lst["position"])
            bim_lst.to_csv("temp/temp.bim",sep="\t",index=False,header=False)

            #then, calculate ld
            ldstore_command="ldstore --bplink temp/temp --bcor temp_corr.bcor --ld-thold {}  --incl-range {}-{} --n-threads {}".format(
            args.ld_treshold,
            range_lower,
            range_upper,
            args.ldstore_threads)
            ldstore_merge_command="ldstore --bcor temp_corr.bcor --merge {}".format(args.ldstore_threads)
            ldstore_extract_info="ldstore --bcor temp_corr.bcor --table ld_table.table "
            subprocess.call(shlex.split(ldstore_command), stderr=subprocess.DEVNULL )
            subprocess.call(shlex.split(ldstore_merge_command),stderr=subprocess.DEVNULL )
            subprocess.call(shlex.split(ldstore_extract_info), stderr=subprocess.DEVNULL )
            #read ld_table.table in
            ld_df_=pd.read_csv("ld_table.table",sep="\s+")
            #let's do this policy:
            #   id1 is the one from our results, i.e. keep only those where id1 is our variants
            #   id2 is the catalog entries, i.e. keep only those
            #   Then, announce the results 
            ld_df_=ld_df_[["chromosome", "index1", "RSID1", "position1", "index2", "RSID2", "position2", "correlation"]]
            ld_df_=ld_df_.merge(df["#variant"],how="inner",left_on="RSID1",right_on="#variant")
            ld_df_=ld_df_.merge(summary_df["#variant"],how="inner",left_on="RSID2",right_on="#variant")
            #print(ld_df_.columns)
            ld_df_=ld_df_.drop(columns=["#variant_x","#variant_y"])
            ld_df=pd.concat((ld_df,ld_df_),axis=0)

        c4="rm -r temp/"
        c5="rm temp_corr.bcor* ld_table.table"
        Popen(shlex.split(c4),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        Popen(shlex.split(c5),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        ld_df=ld_df.merge(summary_df[["#variant","pval"]],left_on="RSID2",right_on="#variant",how="inner")
        ld_df=ld_df.drop(columns=["RSID2","index1","index2"])
        ld_df=ld_df.rename({"RSID_1":"gws_variant"})
        ld_df.to_csv("ld_raport_out.csv",sep="\t",index=False)
    #TODO: make a better representation from the data
    tmp.to_csv(args.raport_out,sep="\t",index=False)

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Compare found GWS results to previously found results")
    parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,help="use 'file' or 'gwascatalog'")
    parser.add_argument("--summary-fpath",dest="summary_files",metavar="FILE",nargs="+",help="comparison summary filepaths")
    parser.add_argument("--endpoints",type=str,nargs="+",help="biological endpoint, as many as summaries")
    parser.add_argument("--build-38",dest="build_38",action="store_true",help="Whether is in GRCh38")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",help="The path for the LD panel to determine what samples are in LD with each other")
    parser.add_argument("--raport-out",dest="raport_out",type=str,default="raport_output.csv",help="Raport output path")
    parser.add_argument("--gwascatalog-pval",default=5e-8,help="P-value cutoff for GWASCatalog searches")
    parser.add_argument("--gwascatalog-width-kb",dest="gwascatalog_pad",type=int,default=25,help="gwascatalog range padding")
    parser.add_argument("--ldstore-threads",type=int,default=4,help="Number of threads to use with ldstore")
    parser.add_argument("--ld-treshold",type=float,default=0.4,help="ld treshold")
    args=parser.parse_args()
    compare(args)