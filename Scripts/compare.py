#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd
import numpy as np
from gws_fetch import create_variant_column,prune_regions
import gwcatalog_api


def compare(args):
    
    #load original file
    df=pd.read_csv(args.compare_fname,sep="\t")
    #if using summary file
    if args.compare_style=="file":
        #load summary file
        for idx in range(0,len(args.summary_files) ):
            summary_df=pd.read_csv(args.summary_files[idx],sep="\t")
            #make sure the variant column exists
            summary_cols=summary_df.columns
            if "#variant" not in summary_cols:
                summary_df.loc[:, "#variant"]=create_variant_column(summary_df)
            #summary_df.loc[:,"endpoint"] = args.endpoint
            if not args.build_38:
                raise NotImplementedError("Non-build 38 raports are not supported yet.")
            else:
                df.loc[:,args.endpoints[idx]]=np.nan
                df.loc[df["#variant"].isin(summary_df["#variant"]),args.endpoints[idx]]=1
                if args.ld_check:
                    raise NotImplementedError("LD calculations not yet implemented")
        df.to_csv(args.raport_out,sep="\t",index=False)
                    
    elif args.compare_style=="gwascatalog":
        #create ranges that contain all of the SNPs
        range_df=df.loc[:,["#chrom","pos"]].copy(deep=True)
        range_df.loc[:,"pos2"]=range_df.loc[:,"pos"]
        range_df=range_df.rename(columns={"pos":"pos_rmin","pos2":"pos_rmax"})
        print(range_df.columns)
        #add 100000 bp to the interval to make the ranges a bit more common
        range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"]-100000
        range_df.loc[:,"pos_rmax"]=range_df.loc[:,"pos_rmax"]+100000
        range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"].clip(lower=0)
        regions=prune_regions(range_df)
        print(regions)
        #use api to get all gwascatalog hits
        result_lst=[]
        for _,region in regions.iterrows():
            result_lst+= gwcatalog_api.parse_output(gwcatalog_api.get_all_associations(chromosome=region["#chrom"],
                bp_lower=region["min"],bp_upper=region["max"])  )
        gwas_df=pd.DataFrame(result_lst)
        #parse hits to a proper form, drop unnecessary information
        #drop links
        #only keep chrom, pos, ref, alt, p-val, association, study
        #print(gwas_df.columns)
        gwas_cols=["base_pair_location","chromosome","hm_beta","p_value","hm_effect_allele","hm_other_allele",
            "trait","study_accession","effect_allele","other_allele","hm_code"]
        gwas_df=gwas_df.loc[:,gwas_cols]
        gwas_rename={"base_pair_location":"pos","chromosome":"#chrom","hm_beta":"beta",
            "p_value":"pval","hm_effect_allele":"hm_alt","hm_other_allele":"hm_ref","effect_allele":"alt","other_allele":"ref"}
        gwas_df=gwas_df.rename(columns=gwas_rename)
        
        #gwas_df.loc[:,"code_column"].update(gwas_df.loc[:,"hm_code"])
        tmp_df=gwas_df[["#chrom","pos","ref","alt","hm_ref","hm_alt","pval","hm_code"]]
        
        #filter out the ones that are not harmonized

        #assuming code is the same as hm_code, we want to filter out 9, 14, 15, 16, 17, 18
        filter_out_codes=[9, 14, 15, 16, 17, 18]
        tmp_df=tmp_df.loc[~tmp_df.loc[:,"hm_code"].isin(filter_out_codes)]
        print(tmp_df)
        tmp_df.to_csv("temp.csv",sep="\t",index=False)

        #change alleles to a strand using the code from commons

        #report hits in some way
        raise NotImplementedError("comparison method '{}' not yet implemented".format(args.compare_style))
    else:
        raise NotImplementedError("comparison method '{}' not yet implemented".format(args.compare_style))

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Compare found GWS results to previously found results")
    parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,help="use 'file' or 'gwascatalog'")
    parser.add_argument("--summary-fpath",dest="summary_files",metavar="FILE",nargs="+",help="comparison summary filepaths")
    parser.add_argument("--endpoints",type=str,nargs="+",help="biological endpoint, as many as summaries")
    parser.add_argument("--build-38",dest="build_38",action="store_true",help="Whether is in GRCh38")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--raport-out",dest="raport_out",type=str,default="raport_output.csv",help="Raport output path")
    args=parser.parse_args()
    compare(args)