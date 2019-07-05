#! /usr/bin/python3

import argparse,shlex,subprocess, glob
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
from autoreporting_utils import *
import gwcatalog_api
import os
from multiprocessing.dummy import Pool as ThreadPool

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
    
def gwcatalog_call_helper(chrom,bp_lower,bp_upper,p_upper):
    val=gwcatalog_api.get_all_associations(chromosome=chrom,
                    bp_lower=bp_lower,bp_upper=bp_upper,p_upper=p_upper)
    if val:
        return gwcatalog_api.parse_output(val)
    return

def map_dataframe(df,columns):
    for _, row in df.iterrows():
        [alt,ref]=map_alleles(row[ columns["alt"] ],row[ columns["ref"] ])
        df.loc[_,"map_ref"]=ref
        df.loc[_,"map_alt"]=alt
    return df

def create_top_level_report(input_df,input_summary_df,efo_traits,columns):
    """
    Create a top level report from which it is easy to see which loci are novel
    In: df, gwascatalog summary df, traits that appear in matching_pheno_gwas_catalog_hits, column names
    Out: Dataframe with a row for every lead variant in df, with columns locus_id chr, start, end, matching_pheno_gwas_catalog_hits other_gwas_hits  
    """ 
    #copy dfs to make sure that the original dataframes are NOT modified 
    df=input_df.copy()
    summary_df=input_summary_df.copy()

    #create mapping column
    summary_df=map_dataframe(summary_df,columns)
    df=map_dataframe(df,columns)

    summary_df.loc[:,"map_variant"]=create_variant_column(summary_df,chrom=columns["chrom"],pos=columns["pos"],ref="map_ref",alt="map_alt")
    df.loc[:,"map_variant"]=create_variant_column(df,chrom=columns["chrom"],pos=columns["pos"],ref="map_ref",alt="map_alt")
    summary_columns=["map_variant","trait","trait_name"]
    summary_df=summary_df.loc[:,summary_columns]
    merged=df.merge(summary_df,on="map_variant",how="left")
    list_of_loci=list(df["locus_id"].unique())
    #compile new simple top level dataframe
    top_level_columns=["locus_id","chr","start","end","matching_pheno_gwas_catalog_hits","other_gwas_hits"]
    top_level_df=pd.DataFrame(columns=top_level_columns)
    for locus_id in list_of_loci:
        #get variants of this locus
        loc_variants=merged.loc[merged["locus_id"]==locus_id,:]
        #chr,start, end
        chrom=loc_variants[columns["chrom"]].values[0]
        start=np.amin(loc_variants[columns["pos"]])
        end=np.amax(loc_variants[columns["pos"]])
        #find all of the traits that have hits
        all_traits=list(loc_variants["trait"].drop_duplicates().dropna())
        trait_dict={}
        for _,row in loc_variants.loc[:,["trait","trait_name"]].drop_duplicates().dropna().iterrows():
            if row["trait_name"]=="NA":
                trait_dict[row["trait"]]=row["trait"]
            else:
                trait_dict[row["trait"]]=row["trait_name"]
        matching_traits=[str(trait_dict[trait]) for trait in all_traits if trait in efo_traits]
        other_traits=[str(trait_dict[trait]) for trait in all_traits if trait not in efo_traits]
        top_level_df=top_level_df.append({"locus_id":locus_id,"chr":chrom,"start":start,"end":end,"matching_pheno_gwas_catalog_hits":";".join(matching_traits),
        "other_gwas_hits":";".join(other_traits)},ignore_index=True)
    return top_level_df
         
        

def compare(args):
    """
    Compares our findings to gwascatalog results or supplied summary statistic files
    """
    columns={"chrom":args.column_labels[0],"pos":args.column_labels[1],"ref":args.column_labels[2],"alt":args.column_labels[3],"pval":args.column_labels[4]}
    #load original file
    df=pd.read_csv(args.compare_fname,sep="\t")
    necessary_columns=[ columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"],"#variant","locus_id","pos_rmin","pos_rmax"]
    df_cols=df.columns.to_list()
    if not all(nec_col in df_cols for nec_col in necessary_columns):
        Exception("GWS variant file {} did not contain all of the necessary columns:\n{} ".format(args.compare_fname,necessary_columns))
    
    #if using summary file
    summary_df_1=pd.DataFrame()
    summary_df_2=pd.DataFrame()
    if args.compare_style in ["file","both"]:
        #load summary files
        #summary_df=pd.DataFrame()
        for idx in range(0,len(args.summary_files) ):
            s_df=pd.read_csv(args.summary_files[idx],sep="\t")
            #make sure the variant column exists
            summary_cols=s_df.columns
            if "#variant" not in summary_cols:
                s_df.loc[:, "#variant"]=create_variant_column(summary_df)
            s_df.loc[:,"trait"] = args.endpoints[idx]
            s_df.loc[:,"trait_name"] = args.endpoints[idx]
            cols=s_df.columns.to_list()
            necessary_columns=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"],"#variant"]
            if not all(nec_col in cols for nec_col in (necessary_columns+["trait","trait_name"])):
                Exception("Summary statistic file {} did not contain all of the necessary columns:\n{} ".format(args.summary_files[idx],necessary_columns))
            #s_df=s_df.loc[:,necessary_columns+["trait","trait_name"]]
            summary_df_1=pd.concat([summary_df_1,s_df],axis=0,sort=True)        
    if args.compare_style in ["gwascatalog","both"]:
        gwas_df=None
        if os.path.exists("gwas_out_mapping.csv") and args.cache_gwas:
            print("reading gwas results from gwas_out_mapping.csv...")
            gwas_df=pd.read_csv("gwas_out_mapping.csv",sep="\t")
        else:
            rm_gwas_out="rm gwas_out_mapping.csv"
            subprocess.call(shlex.split(rm_gwas_out),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            #NOTE: change range to be the same that was defined in fetching, either the width or first and last variant in the group if ld
            #create ranges that contain all of the SNPs
            range_df=df.loc[:,[columns["chrom"],columns["pos"] ]].copy(deep=True)
            range_df.loc[:,"pos2"]=range_df.loc[:,columns["pos"] ]
            range_df=range_df.rename(columns={columns["pos"]:"pos_rmin","pos2":"pos_rmax"})
            
            pad=args.gwascatalog_pad*1000
            range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"]-pad
            range_df.loc[:,"pos_rmax"]=range_df.loc[:,"pos_rmax"]+pad
            range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"].clip(lower=0)
            regions=prune_regions(range_df,columns=columns)
            #use api to get all gwascatalog hits
            #create call data list
            data_lst=[]
            for _,region in regions.iterrows():
                data_lst.append([region[ columns["chrom"] ], region["min"],region["max"],args.gwascatalog_pval])
            #create worker pool for multithreaded api calls
            r_lst=None
            with ThreadPool(args.gwascatalog_threads) as pool:
                r_lst=pool.starmap(gwcatalog_call_helper,data_lst)
            #remove empties
            r_lst=[r for r in r_lst if r != None]
            result_lst=[]
            for sublst in r_lst:
                result_lst=result_lst+sublst
            gwas_df=pd.DataFrame(result_lst)
            gwas_df["trait"]=gwas_df["trait"].apply(lambda x:x[0])
            gwas_df.to_csv("gwas_out_mapping.csv",sep="\t")
        #parse hits to a proper form, drop unnecessary information
        gwas_cols=["base_pair_location","chromosome","p_value","hm_effect_allele","hm_other_allele",
            "trait","study_accession","hm_code","hm_beta","hm_effect_allele_frequency"]
        gwas_df=gwas_df.loc[:,gwas_cols]
        gwas_rename={"base_pair_location":columns["pos"],"chromosome":columns["chrom"],"hm_beta":"beta",
            "p_value":columns["pval"],"hm_code":"code","hm_effect_allele":columns["alt"],"hm_other_allele":columns["ref"],"hm_effect_allele_frequency":"af"}
        gwas_df=gwas_df.rename(columns=gwas_rename)
        
        tmp_df=gwas_df[[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"],"code","beta","af","trait"]]
        #assuming code is the same as hm_code, we want to filter out 9, 14, 15, 16, 17, 18
        filter_out_codes=[9, 14, 15, 16, 17, 18]
        tmp_df=tmp_df.loc[~tmp_df.loc[:,"code"].isin(filter_out_codes)]
        tmp_df.loc[:,"#variant"]=create_variant_column(tmp_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        #change alleles to a strand using the code from commons
        summary_df_2=tmp_df
        #create list of unique traits
        unique_efos=list(summary_df_2["trait"].unique())
        trait_name_map={}
        for key in unique_efos:
            trait_name_map[key]=gwcatalog_api.get_trait_name(key)
        summary_df_2.loc[:,"trait_name"]=summary_df_2.loc[:,"trait"].apply(lambda x: trait_name_map[x])
        summary_df_2=summary_df_2.drop_duplicates(subset=["#variant","trait"])
    else:
        raise NotImplementedError("comparison method '{}' not yet implemented".format(args.compare_style))
    summary_df=pd.concat([summary_df_1,summary_df_2],sort=True)
    #top level df
    top_df=create_top_level_report(df,summary_df_2,args.efo_traits,columns)
    top_df.to_csv(args.top_report_out,sep="\t",index=False)

    #now we should have df and summary_df
    necessary_columns=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"],"#variant","trait","trait_name"]
    summary_df=summary_df.loc[:,necessary_columns]
    #summary_df.to_csv("summary_df.csv",sep="\t",index=False)

    summary_df=map_dataframe(summary_df,columns)
    df=map_dataframe(df,columns)

    summary_df.loc[:,"map_variant"]=create_variant_column(summary_df,chrom=columns["chrom"],pos=columns["pos"],ref="map_ref",alt="map_alt")
    df.loc[:,"map_variant"]=create_variant_column(df,chrom=columns["chrom"],pos=columns["pos"],ref="map_ref",alt="map_alt")
    #df.to_csv("df.csv",sep="\t",index=False)
    tmp=pd.merge(df,summary_df.loc[:,["#variant","map_variant",columns["pval"],"trait","trait_name"]],how="left",on="map_variant")
    tmp=tmp.drop(columns=["map_variant","map_ref","map_alt"])
    tmp=tmp.rename(columns={"#variant_x":"#variant","#variant_y":"#variant_hit","pval_x":columns["pval"],"pval_y":"pval_trait"})
    tmp=tmp.sort_values(by=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],"#variant"])
    tmp.to_csv(args.raport_out,sep="\t",index=False)
    
    if args.ld_check:
        #if no groups in base data
        if ("pos_rmin" not in df.columns.to_list()) or ("pos_rmax" not in df.columns.to_list()):
            Exception("ld calculation not supported without grouping. Please supply the flag --group to main.py or gws_fetch.py.") 

        #get variant list, i.e. the list of variants that consists of gws results and summary results
        var_cols=["#variant",columns["pos"],columns["chrom"],columns["ref"],columns["alt"]]
        var_rename={"#variant":"RSID",columns["pos"]:"position",columns["chrom"]:"chromosome",columns["ref"]:"A_allele",columns["alt"]:"B_allele"}
        var_lst_df=pd.concat([summary_df.loc[:,var_cols].rename(columns=var_rename),df.loc[:,var_cols].rename(columns=var_rename)  ])
        var_lst_df=var_lst_df.drop_duplicates(subset=["RSID"])
        unique_locus_list=df["locus_id"].unique()
        ld_df=pd.DataFrame()
        for locus in unique_locus_list:
            
            chromosome=df.loc[df["#variant"]==locus,columns["chrom"] ].unique()[0]
            print("Chromosome {}, group {} ld computation".format(chromosome,locus))
            #get group range
            r_max=df.loc[df["#variant"]==locus,"pos_rmax"].values[0]
            r_min=df.loc[df["#variant"]==locus,"pos_rmin"].values[0]
            if r_max == r_min:
                continue
            #calculate ld for that group
            ldstore_command="ldstore --bplink {}_{} --bcor temp_corr.bcor --ld-thold {}  --incl-range {}-{} --n-threads {}".format(
            args.ld_chromosome_panel,
            chromosome,
            args.ld_treshold**0.5,#due to ldstore taking in |r|, not r2
            r_min,
            r_max,
            args.ldstore_threads)
            pr = subprocess.run(shlex.split(ldstore_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='ASCII' )
            if pr.returncode!=0:
                print("LDSTORE FAILURE for locus {}".format(locus)  )
                print(pr.stdout)
                continue 
            #merge the file
            ldstore_merge_command="ldstore --bcor temp_corr.bcor --merge {}".format(args.ldstore_threads)
            subprocess.call(shlex.split(ldstore_merge_command),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL )
            pr = subprocess.run(shlex.split(ldstore_merge_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='ASCII' )
            if pr.returncode!=0:
                print("LDSTORE FAILURE for locus {}".format(locus)  )
                print(pr.stdout)
                continue
            #create list of variants of interest.
            #First, create list of variants that we want to etract the correlations for.
            #This should be all of the variants in the summary df inside the correct range, and the group variants.
            extract_df_1=df.loc[df["locus_id"]==locus,:].copy()
            extract_df_2=summary_df.loc[(summary_df[columns["pos"]] <=r_max) & (summary_df[columns["pos"]] >=r_min)  ,:].copy()
            extract_df=pd.concat([extract_df_1,extract_df_2],sort=True).loc[:,var_cols].rename(columns=var_rename).drop_duplicates().sort_values(by="position")
            extract_df.to_csv("var_lst",sep=" ",index=False)
            #extract variants of interest 
            ldstore_extract_command="ldstore --bcor temp_corr.bcor --table ld_table.table --incl-variants var_lst"
            pr = subprocess.run(shlex.split(ldstore_extract_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='ASCII' )
            if pr.returncode!=0:
                print("LDSTORE FAILURE for locus {}".format(locus)  )
                print(pr.stdout)
                continue
            #read ld_table, make it so rsid1 is for our variants and rsid2 for summary variants
            ld_table=pd.read_csv("ld_table.table",sep="\s+").loc[:,["chromosome","RSID1","RSID2","correlation"]]
            ld_table2=ld_table.copy()
            ld_table2=ld_table2.rename(columns={"RSID1":"temp_rsid2","RSID2":"temp_rsid1"})
            ld_table2=ld_table2.rename(columns={"temp_rsid2":"RSID2","temp_rsid1":"RSID1"})
            ld=pd.concat([ld_table,ld_table2],sort=True)
            #filter
            ld=ld.loc[ld["RSID1"].isin(extract_df_1["#variant"].values),:]
            ld=ld.loc[ld["RSID2"].isin(extract_df_2["#variant"].values),:]
            ld.loc[:,"r2"]=ld["correlation"]*ld["correlation"]
            ld=ld.drop(columns=["correlation"])
            ld_df=pd.concat([ld_df,ld],sort=True)
        #print("ld_df columns:{}".format(ld_df.columns))
        #print("summary_df columns:{}".format(summary_df.columns))
        c5="rm ld_table.table var_lst"
        corr_files=glob.glob("temp_corr.*")
        Popen(shlex.split(c5)+corr_files,stderr=subprocess.DEVNULL)
        ld_df=ld_df.drop_duplicates(subset=["RSID1","RSID2"],keep="first")
        ld_df=ld_df.merge(summary_df.loc[:,["#variant","trait","trait_name"]].rename(columns={"#variant":"RSID2"}),how="inner",on="RSID2")
        ld_df=ld_df.loc[ld_df["RSID1"].isin(df["#variant"].values),:]
        #ld_df.to_csv("ld_out.csv",sep="\t",index=False)
        if not ld_df.empty:
            ld_out=df.merge(ld_df,how="left",left_on="#variant",right_on="RSID1")
            ld_out=ld_out.drop(columns=["RSID1","map_variant","map_ref","map_alt"]).rename(columns={"{}_x".format(columns["pval"]):columns["pval"],"{}_y".format(columns["pval"]):"pval_trait","RSID2":"#variant_hit"})
            ld_out.to_csv(args.ld_raport_out,sep="\t",index=False)
        else:
            print("No variants in ld found, no {} produced.".format(args.ld_raport_out))

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Compare found GWS results to previously found results")
    parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,default="gwascatalog",help="use 'file', 'gwascatalog' or 'both'")
    parser.add_argument("--summary-fpath",dest="summary_files",metavar="FILE",nargs="+",help="comparison summary filepaths")
    parser.add_argument("--endpoints",type=str,nargs="+",help="biological endpoint, as many as summaries")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--ld-chromosome-panel-path",dest="ld_chromosome_panel",help="Path to ld panel, where each chromosome is separated. If path is 'path/panel_#chrom.bed', input 'path/panel' ")
    parser.add_argument("--raport-out",dest="raport_out",type=str,default="raport_out.csv",help="Raport output path")
    parser.add_argument("--ld-raport-out",dest="ld_raport_out",type=str,default="ld_raport_out.csv",help="LD check raport output path")
    parser.add_argument("--gwascatalog-pval",default=5e-8,help="P-value cutoff for GWASCatalog searches")
    parser.add_argument("--gwascatalog-width-kb",dest="gwascatalog_pad",type=int,default=25,help="gwascatalog range padding")
    parser.add_argument("--gwascatalog-threads",dest="gwascatalog_threads",type=int,default=4,help="Number of concurrent queries to GWAScatalog API. Default 4. Increase if the gwascatalog api takes too long.")
    parser.add_argument("--ldstore-threads",type=int,default=4,help="Number of threads to use with ldstore. Default 4")
    parser.add_argument("--ld-treshold",type=float,default=0.4,help="ld treshold")
    parser.add_argument("--cache-gwas",action="store_true",help="save gwascatalog results into gwas_out_mapping.csv and load them from there if it exists. Use only for testing.")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval"],help="Names for data file columns. Default is '#chrom pos ref alt pval'.")
    parser.add_argument("--top-report-out",dest="top_report_out",type=str,default="top_report.csv",help="Top level report filename.")
    parser.add_argument("--efo-codes",dest="efo_traits",type=str,nargs="+",default=[],help="Specific EFO codes to look for in the top level raport")
    args=parser.parse_args()
    compare(args)
