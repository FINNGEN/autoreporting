#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd
import numpy as np
from gws_fetch import create_variant_column,prune_regions
import gwcatalog_api


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
        #create ranges that contain all of the SNPs
        range_df=df.loc[:,["#chrom","pos"]].copy(deep=True)
        range_df.loc[:,"pos2"]=range_df.loc[:,"pos"]
        range_df=range_df.rename(columns={"pos":"pos_rmin","pos2":"pos_rmax"})
        #add 100000 bp to the interval to make the ranges a bit more common
        range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"]-25000
        range_df.loc[:,"pos_rmax"]=range_df.loc[:,"pos_rmax"]+25000
        range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"].clip(lower=0)
        regions=prune_regions(range_df)
        #use api to get all gwascatalog hits
        print(regions)
        result_lst=[]
        for _,region in regions.iterrows():
            result_lst+= gwcatalog_api.parse_output(gwcatalog_api.get_all_associations(chromosome=region["#chrom"],
                bp_lower=region["min"],bp_upper=region["max"],p_upper=args.gwascatalog_pval)  )
        gwas_df=pd.DataFrame(result_lst)
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
        #print(tmp_df.columns)
        #tmp_df.to_csv("temp_before_mapping.csv",sep="\t",index=False)
        tmp_df.loc[:,"#variant"]=create_variant_column(tmp_df,chrom="#chrom",pos="pos",ref="ref",alt="alt")
        #change alleles to a strand using the code from commons
        summary_df=tmp_df
        #print(summary_df.loc[0,"trait"][0])
        #create list of unique traits
        summary_df.loc[:,"trait"]=summary_df.loc[:,"trait"].apply(lambda x:",".join(x))
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
        #gather the variant list
        var_cols=["#variant","pos","#chrom","ref","alt"]
        var_rename={"#variant":"RSID","pos":"position","#chrom":"chromosome","ref":"A_allele","alt":"B_allele"}
        var_lst_df=pd.concat([summary_df.loc[:,var_cols].rename(columns=var_rename),df.loc[:,var_cols].rename(columns=var_rename)  ])
        var_lst_df.to_csv("variants_for_ld.csv",sep="\t",index=False)
        raise NotImplementedError("LD comparison not yet implemented".format(args.compare_style))

    #TODO: make a better representation from the data
    tmp.to_csv("temp_joined_mapping.csv",sep="\t",index=False)

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Compare found GWS results to previously found results")
    parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,help="use 'file' or 'gwascatalog'")
    parser.add_argument("--summary-fpath",dest="summary_files",metavar="FILE",nargs="+",help="comparison summary filepaths")
    parser.add_argument("--endpoints",type=str,nargs="+",help="biological endpoint, as many as summaries")
    parser.add_argument("--build-38",dest="build_38",action="store_true",help="Whether is in GRCh38")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--ld-panel-path",dest="ld_panel",help="The path for the LD panel to determine what samples are in LD with each other")
    parser.add_argument("--raport-out",dest="raport_out",type=str,default="raport_output.csv",help="Raport output path")
    parser.add_argument("--gwascatalog-pval",default=5e-8,help="P-value cutoff for GWASCatalog searches")
    
    args=parser.parse_args()
    compare(args)