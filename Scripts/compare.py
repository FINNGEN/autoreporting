#! /usr/bin/python3

import argparse,shlex,subprocess, glob
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
from autoreporting_utils import *
import gwcatalog_api
import os

def map_alleles(a1,a2):
    """
    Flips alleles to the A strand if necessary and orders them lexicogaphically
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
    """
    Compares our findings to gwascatalog results or supplied summary statistic files
    """
    #load original file
    df=pd.read_csv(args.compare_fname,sep="\t")
    necessary_columns=["#chrom","pos","ref","alt","pval","#variant","locus_id","pos_rmin","pos_rmax"]
    df_cols=df.columns.to_list()
    if not all(nec_col in df_cols for nec_col in necessary_columns):
        Exception("GWS variant file {} did not contain all of the necessary columns:\n{} ".format(args.compare_fname,necessary_columns))
    
    #if using summary file
    if args.compare_style=="file":
        #load summary files
        summary_df=pd.DataFrame()
        for idx in range(0,len(args.summary_files) ):
            s_df=pd.read_csv(args.summary_files[idx],sep="\t")
            #make sure the variant column exists
            summary_cols=s_df.columns
            if "#variant" not in summary_cols:
                s_df.loc[:, "#variant"]=create_variant_column(summary_df)
            
            s_df.loc[:,"trait"] = args.endpoints[idx]
            s_df.loc[:,"trait_name"] = args.endpoints[idx]
            if not args.build_38:
                raise NotImplementedError("Non-build 38 raports are not supported.")
            cols=s_df.columns.to_list()
            necessary_columns=["#chrom","pos","ref","alt","pval","#variant"]
            if not all(nec_col in cols for nec_col in (necessary_columns+["trait","trait_name"])):
                Exception("Summary statistic file {} did not contain all of the necessary columns:\n{} ".format(args.summary_files[idx],necessary_columns))
            #s_df=s_df.loc[:,necessary_columns+["trait","trait_name"]]
            summary_df=pd.concat([summary_df,s_df],axis=0)        
                    
    elif args.compare_style=="gwascatalog":
        gwas_df=None
        if os.path.exists("gwas_out_mapping.csv"):
            gwas_df=pd.read_csv("gwas_out_mapping.csv",sep="\t")
        else:
            #NOTE: change range to be the same that was defined in fetching, either the width or first and last variant in the group if ld
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
            #print(regions)
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
        #tmp_df.to_csv("unfiltered.csv",sep="\t",index=False)
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
    necessary_columns=["#chrom","pos","ref","alt","pval","#variant","trait","trait_name"]
    summary_df=summary_df.loc[:,necessary_columns]
    summary_df.to_csv("summary_df.csv",sep="\t",index=False)

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
    #df.to_csv("df.csv",sep="\t",index=False)
    tmp=pd.merge(df,summary_df.loc[:,["#variant","map_variant","pval","trait","trait_name"]],how="left",on="map_variant")
    tmp=tmp.drop(columns=["map_variant","map_ref","map_alt"])
    tmp=tmp.rename(columns={"#variant_x":"#variant","#variant_y":"#variant_hit","pval_x":"pval","pval_y":"pval_trait"})
    tmp=tmp.sort_values(by=["#chrom","pos","ref","alt","#variant"])
    tmp.to_csv(args.raport_out,sep="\t",index=False)
    
    if args.ld_check:
        #if no groups in base data
        if ("pos_rmin" not in df.columns.to_list()) or ("pos_rmax" not in df.columns.to_list()):
            Exception("ld calculation not supported without grouping. Please supply the flag --group to main.py or gws_fetch.py.") 

        #preprocess found variants, i.e. divide variants in both our and gwascatalog results to chromosomes
        #Then, we need to build the bim files, I guess, or I could just process some ld panels ready and push them into a bucket
        #Then, for each chromosome set, we want to separate the clear ranges from each other, I guess
        #then, we run ldstore, export correlations, and see if any of the gwascatalog variants are in ld with the gws variants
        #NOTE: ldstore only seems to take bim file into account if we use range.

        #get variant list, i.e. the list of variants that consists of gws results and summary results
        var_cols=["#variant","pos","#chrom","ref","alt"]
        var_rename={"#variant":"RSID","pos":"position","#chrom":"chromosome","ref":"A_allele","alt":"B_allele"}
        var_lst_df=pd.concat([summary_df.loc[:,var_cols].rename(columns=var_rename),df.loc[:,var_cols].rename(columns=var_rename)  ])
        var_lst_df=var_lst_df.drop_duplicates(subset=["RSID"])
        #get unique chromosomes in the results
        unique_chrom_list=df["#chrom"].unique()
        unique_locus_list=df["locus_id"].unique()
        c1="mkdir temp"
        Popen(shlex.split(c1),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        c2="ln -s ../{}.bed ./temp/temp.bed".format(args.ld_panel_path)
        c3="ln -s ../{}.fam ./temp/temp.fam".format(args.ld_panel_path)
        
        Popen(shlex.split(c2),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        Popen(shlex.split(c3),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        ld_df=pd.DataFrame()
        for locus in unique_locus_list:
            chromosome=df.loc[df["#variant"]==locus,"#chrom"].unique()[0]
            print("Chromosome {}, group {} ld computation".format(chromosome,locus))
            #build bim file, containing both the variants in df and summary_df
            bim_lst=None
            #get only variants in that group

            bim_lst=var_lst_df.loc[var_lst_df["chromosome"]==chromosome,:].copy()
            #filter those that are in the group
            r_max=df.loc[df["#variant"]==locus,"pos_rmax"].values[0]
            r_min=df.loc[df["#variant"]==locus,"pos_rmin"].values[0]
            bim_lst=bim_lst.loc[(bim_lst["position"]>=r_min)&(bim_lst["position"]<=r_max),:]

            bim_lst.loc[:,"cm"]=0
            bim_lst=bim_lst.loc[:,["chromosome","RSID","cm","position","A_allele","B_allele"]]
            #temporary solution: create temp folder, create symbolic links to other files, write bim file to temp folder,
            #do the magic, and delete temp folder afterwards. No modification to existing files, and not a lot of disk space usage
            #still ugly, but at least we don't modify data that is supplied to this script
            range_upper=np.max(bim_lst["position"])
            range_lower=np.min(bim_lst["position"])
            print("Variant amount: {}".format(bim_lst.shape[0]))
            if (range_upper==range_lower) or bim_lst.shape[0]<=1:
                continue
            bim_lst.to_csv("temp/temp.bim",sep="\t",index=False,header=False)
            #threads set to 1 because it seems our playing around with bim files does not seems to like threads
            #then, calculate ld
            ldstore_command="ldstore --bplink temp/temp --bcor temp_corr.bcor --ld-thold {}  --incl-range {}-{} --n-threads {}".format(
            args.ld_treshold,
            range_lower,
            range_upper,
            1)
            #args.ldstore_threads)
            #ldstore_merge_command="ldstore --bcor temp_corr.bcor --merge {}".format(args.ldstore_threads)
            ldstore_extract_info="ldstore --bcor temp_corr.bcor_1 --table ld_table.table "
            retval=subprocess.call(shlex.split(ldstore_command),stdout=subprocess.DEVNULL )
            if retval!= 0:
                continue
            #subprocess.call(shlex.split(ldstore_merge_command),stdout=subprocess.DEVNULL )
            retval=subprocess.call(shlex.split(ldstore_extract_info),stdout=subprocess.DEVNULL )
            #read ld_table.table in
            ld_df_=pd.read_csv("ld_table.table",sep="\s+")
            ld_df_=ld_df_.loc[:,["RSID1", "RSID2", "correlation"]]
            #constrain id1 to only contain our variants, and id2 to only contain gwascatalog variants
            ld_df_=ld_df_.merge(df["#variant"].drop_duplicates(),how="inner",left_on="RSID1",right_on="#variant")
            ld_df_=ld_df_.merge(summary_df.loc[:,["#variant","pval","trait","trait_name"]],how="inner",left_on="RSID2",right_on="#variant")
            ld_df_=ld_df_.drop(columns=["#variant_x","#variant_y"])
            
            ld_df=pd.concat((ld_df,ld_df_),axis=0)
        c4="rm -r temp/"
        c5="rm ld_table.table "
        corr_files=glob.glob("temp_corr.*")
        Popen(shlex.split(c4),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        Popen(shlex.split(c5)+corr_files,stderr=subprocess.DEVNULL)
        
        ld_out=df.merge(ld_df,how="left",left_on="#variant",right_on="RSID1")
        ld_out=ld_out.drop(columns=["RSID1","map_variant","map_ref","map_alt"]).rename(columns={"pval_x":"pval","pval_y":"pval_trait","RSID2":"#variant_hit"})
        ld_out.to_csv("ld_raport_out.csv",sep="\t",index=False)

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Compare found GWS results to previously found results")
    parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,default="gwascatalog",help="use 'file' or 'gwascatalog'")
    parser.add_argument("--summary-fpath",dest="summary_files",metavar="FILE",nargs="+",help="comparison summary filepaths")
    parser.add_argument("--endpoints",type=str,nargs="+",help="biological endpoint, as many as summaries")
    parser.add_argument("--build-38",dest="build_38",action="store_true",help="Whether supplied comparison summary files are in GRCh38")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",help="The path for the LD panel to determine what samples are in LD with each other")
    parser.add_argument("--raport-out",dest="raport_out",type=str,default="raport_output.csv",help="Raport output path")
    parser.add_argument("--gwascatalog-pval",default=5e-8,help="P-value cutoff for GWASCatalog searches")
    parser.add_argument("--gwascatalog-width-kb",dest="gwascatalog_pad",type=int,default=25,help="gwascatalog range padding")
    parser.add_argument("--ldstore-threads",type=int,default=4,help="Number of threads to use with ldstore")
    parser.add_argument("--ld-treshold",type=float,default=0.4,help="ld treshold")
    args=parser.parse_args()
    compare(args)