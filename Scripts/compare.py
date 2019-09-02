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
    Author: Pietro, with small changes by Arto
    """
    allele_dict={"T":"A","C":"G","G":"C"}
     # check if the A variant is present
    if 'A' not in a1 + a2 and 'a' not in a1+a2:
        # for both/ref and alt map them to the A strand and order each one lexicographically
        a1 = ''.join([allele_dict[elem.upper()] for elem in a1])
        a2 = ''.join([allele_dict[elem.upper()] for elem in a2])
    # further sorting
    return sorted([a1,a2])

def map_column(df,col_name,columns):
    df_=df.copy()
    if df_.empty:
        cols=list(df_.columns)
        cols.append(col_name)
        return pd.DataFrame(columns=cols)
    df_[["temp_map_ref","temp_map_alt"]]=df_.loc[:,[columns["ref"],columns["alt"] ]].apply(lambda x: map_alleles(*x),axis=1,result_type="expand")
    df_[col_name]=create_variant_column(df_,chrom=columns["chrom"],pos=columns["pos"],ref="temp_map_ref",alt="temp_map_alt")
    df_=df_.drop(columns=["temp_map_ref","temp_map_alt"])
    return df_

def indel_helper(row, chrom,pos,ref,alt):
    row["ref"] = ref
    row["alt"] = alt
    row["chrom"] = chrom
    row["pos"] = pos
    return row

def solve_indels(indel_df,df,columns):
    """ Solve exact matches for indels where they are missing one basepair 
    """
    out_df=pd.DataFrame(columns=indel_df.columns)
    for _, row in indel_df.iterrows():
        #check if our dataframe has an exact match.
        rowpos=int(row["pos"])-1
        rowchrom=str(row["chrom"])
        possible_matches=df.loc[(df[columns["chrom"]].astype(str) == rowchrom)&(df[columns["pos"]] == rowpos )  ,:].copy()
        for __, row2 in possible_matches.iterrows():
            #if alleles are a match s.t. -/TCGA == T/TTCGA, add that to output df and break from the loop
            a1=row["ref"]
            a2=row["alt"]
            b1=row2[columns["ref"]]
            b2=row2[columns["alt"]]
            n_row=row.copy()
            if a1=="-":
                if len(b1)==1 and b2[1:] == a2:
                    n_row=indel_helper(n_row,row2[columns["chrom"]],row2[columns["pos"]],b1,b2)
                    out_df=out_df.append(n_row,sort=True)
                    break
                elif len(b2)==1 and b1[1:] == a2:
                    n_row=indel_helper(n_row,row2[columns["chrom"]],row2[columns["pos"]],b2,b1)
                    out_df=out_df.append(n_row,sort=True)
                    break
            elif a2=="-":
                if len(b1)==1 and b2[1:] == a1:
                    n_row=indel_helper(n_row,row2[columns["chrom"]],row2[columns["pos"]],b2,b1)
                    out_df=out_df.append(n_row,sort=True)
                    break
                elif len(b2)==1 and b1[1:] == a1:
                    n_row=indel_helper(n_row,row2[columns["chrom"]],row2[columns["pos"]],b1,b2)
                    out_df=out_df.append(n_row,sort=True)
                    break
            #else, continue
    return out_df

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
    df=map_column(df,"map_variant",columns)
    if df.empty:
        return pd.DataFrame(columns=["locus_id","chr","start","end","matching_pheno_gwas_catalog_hits","other_gwas_hits"])
    if not summary_df.empty:
        summary_df=map_column(summary_df,"map_variant",columns)
        summary_df=summary_df.loc[:,["map_variant","trait","trait_name"]]
        merged=df.merge(summary_df,on="map_variant",how="left")
    else:
        merged=df.copy()
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
        if not summary_df.empty:
            all_traits=list(loc_variants["trait"].drop_duplicates().dropna())
            trait_dict={}
            for row in loc_variants.loc[:,["trait","trait_name"]].drop_duplicates().dropna().itertuples():
                if row.trait_name =="NA":
                    trait_dict[row.trait]=row.trait
                else:
                    trait_dict[row.trait]=row.trait_name
            matching_traits=[str(trait_dict[trait]) for trait in all_traits if trait in efo_traits]
            other_traits=[str(trait_dict[trait]) for trait in all_traits if trait not in efo_traits]
            top_level_df=top_level_df.append({"locus_id":locus_id,"chr":chrom,"start":start,"end":end,"matching_pheno_gwas_catalog_hits":";".join(matching_traits),
            "other_gwas_hits":";".join(other_traits)},ignore_index=True)
        else:
            top_level_df=top_level_df.append({"locus_id":locus_id,"chr":chrom,"start":start,"end":end,"matching_pheno_gwas_catalog_hits":"",
            "other_gwas_hits":""},ignore_index=True)
    return top_level_df
         
def load_summary_files(summary_fpath,endpoint_fpath,columns):
    necessary_columns=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"],"#variant","trait","trait_name"]
    summary_df_1=pd.DataFrame(columns=necessary_columns)
    with open(summary_fpath,"r") as f:
        s_paths=f.readlines()
        s_paths=[s.strip("\n").strip() for s in s_paths]
    with open(endpoint_fpath,"r") as f:
        endpoints=f.readlines()
        endpoints=[s.strip("\n").strip() for s in endpoints]
    if len(s_paths)!=len(endpoints):
        raise RuntimeError("summary file amount and endpoint amounts are not equal. {} =/= {}".format(len(s_paths),len(endpoints)))
    for idx in range(0,len(s_paths) ):
        s_path=s_paths[idx]
        endpoint=endpoints[idx]
        try:
            s_df=pd.read_csv(s_path,sep="\t")
        except FileNotFoundError as e:
            raise FileNotFoundError("File {} does not exist. Given as line {} in file {} given in argument '--summary-fpath' .".format(s_path,idx+1,summary_fpath))
        #make sure the variant column exists
        summary_cols=s_df.columns
        if "#variant" not in summary_cols:
            s_df.loc[:, "#variant"]=create_variant_column(s_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
        s_df.loc[:,"trait"] = endpoint
        s_df.loc[:,"trait_name"] = endpoint
        cols=s_df.columns.to_list()
        if not all(nec_col in cols for nec_col in necessary_columns):
            Exception("Summary statistic file {} did not contain all of the necessary columns:\n{} ".format(args.summary_files[idx],necessary_columns))
        summary_df_1=pd.concat([summary_df_1,s_df],axis=0,sort=True)
    summary_df_1=summary_df_1.loc[:,necessary_columns].reset_index(drop=True)
    return summary_df_1
    
def load_api_summaries(df, gwascatalog_pad, gwascatalog_pval,gwapi,gwascatalog_threads, columns):
    range_df=df.loc[:,[columns["chrom"],columns["pos"] ]].copy(deep=True)
    range_df.loc[:,"pos2"]=range_df.loc[:,columns["pos"] ]
    range_df=range_df.rename(columns={columns["pos"]:"pos_rmin","pos2":"pos_rmax"})

    pad=gwascatalog_pad*1000
    range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"]-pad
    range_df.loc[:,"pos_rmax"]=range_df.loc[:,"pos_rmax"]+pad
    range_df.loc[:,"pos_rmin"]=range_df.loc[:,"pos_rmin"].clip(lower=0)
    regions=prune_regions(range_df,columns=columns)
    #use api to get all gwascatalog hits
    data_lst=[]
    for _,region in regions.iterrows():
        data_lst.append([region[ columns["chrom"] ], region["min"],region["max"],gwascatalog_pval])
    r_lst=None
    threads=gwascatalog_threads
    with ThreadPool(threads) as pool:
        r_lst=pool.starmap(gwapi.get_associations,data_lst)
    #remove empties, flatten
    r_lst=[r for r in r_lst if r != None]
    result_lst=[i for sublist in r_lst for i in sublist]
    gwas_df=pd.DataFrame(result_lst)
    #resolve indels
    if not gwas_df.empty:
        indel_idx=(gwas_df["ref"]=="-")|(gwas_df["alt"]=="-")
        indels=solve_indels(gwas_df.loc[indel_idx,:],df,columns)
        gwas_df=gwas_df.loc[~indel_idx,:]
        gwas_df=pd.concat([gwas_df,indels]).reset_index(drop=True)
    return gwas_df

def extract_ld_variants(df,summary_df,locus,args,columns):
    if df.loc[df["locus_id"]==locus,"pos_rmax"].shape[0]<=1:
        return
    chromosome=df.loc[df["#variant"]==locus,columns["chrom"] ].unique()[0]
    print("Chromosome {}, group {} ld computation, variant amount {}".format(chromosome,locus,df.loc[df["locus_id"]==locus,"pos_rmax"].shape[0]))
    #get group range
    r_max=df.loc[df["locus_id"]==locus,"pos_rmax"].values[0]
    r_min=df.loc[df["locus_id"]==locus,"pos_rmin"].values[0]
    if r_max == r_min:
        return
    #calculate ld for that group
    threads=min(args.ldstore_threads,df.loc[df["locus_id"]==locus,"pos_rmax"].shape[0]) # ldstore silently errors out when #variants < #threads, do this to alleviate that 
    ldstore_command="ldstore --bplink {}temp_chrom --bcor {}temp_corr.bcor --ld-thold {}  --incl-range {}-{} --n-threads {}".format(
        args.prefix, args.prefix, args.ld_treshold**0.5, r_min, r_max, threads)
    pr = subprocess.run(shlex.split(ldstore_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='ASCII' )
    if pr.returncode!=0:
        print("LDSTORE FAILURE for locus {}".format(locus)  )
        print(pr.stdout)
        return
    #merge the file
    ldstore_merge_command="ldstore --bcor {}temp_corr.bcor --merge {}".format(args.prefix, threads)
    pr = subprocess.run(shlex.split(ldstore_merge_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='ASCII' )
    if pr.returncode!=0:
        print("LDSTORE FAILURE for locus {}".format(locus)  )
        print(pr.stdout)
        return
    #check if ldstore merged file exists, if not, throw an error
    if not os.path.exists("{}temp_corr.bcor".format(args.prefix)):
        raise FileNotFoundError("The LD correlation file {} does not exist. Check that the chromosome index is correct, e.g. 23 in both summary statistic and LD panel instead of 23 and X.".format(args.prefix+"temp_corr.bcor"))
    #create list of variants of interest.
    var_cols=["#variant",columns["pos"],columns["chrom"],columns["ref"],columns["alt"]]
    var_rename={"#variant":"RSID",columns["pos"]:"position",columns["chrom"]:"chromosome",columns["ref"]:"A_allele",columns["alt"]:"B_allele"}
    extract_df_1=df.loc[df["locus_id"]==locus,:].copy()
    extract_df_2=summary_df.loc[(summary_df[columns["pos"]] <=r_max) & (summary_df[columns["pos"]] >=r_min)  ,:].copy()
    extract_df=pd.concat([extract_df_1,extract_df_2],sort=True).loc[:,var_cols].rename(columns=var_rename).drop_duplicates().sort_values(by="position")
    ### different datatypes caused drop duplicates to not recognize duplicates and ldstore failed in the next step
    extract_df = extract_df.astype(str)
    extract_df = extract_df.drop_duplicates()
    var_lst_name="{}var_lst".format(args.prefix)
    extract_df.to_csv(var_lst_name,sep=" ",index=False)
    #extract variants of interest 
    ldstore_extract_command="ldstore --bcor {}temp_corr.bcor --table {}ld_table.table --incl-variants {}".format(args.prefix, args.prefix, var_lst_name)
    pr = subprocess.run(shlex.split(ldstore_extract_command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='ASCII' )
    if pr.returncode!=0:
        print("LDSTORE FAILURE for locus {}".format(locus)  )
        print(pr.stdout)
        return
    #read ld_table, make it so rsid1 is for our variants and rsid2 for summary variants
    if not os.path.exists("{}ld_table.table".format(args.prefix)):
        raise FileNotFoundError("The LD association table {} does not exist. Check that the chromosome index is correct, e.g. 23 in both summary statistic and LD panel instead of 23 and X.".format(args.prefix+"ld_table.table"))
    ld_table=pd.read_csv("{}ld_table.table".format(args.prefix),sep="\s+").loc[:,["chromosome","RSID1","RSID2","correlation"]]
    if ld_table.empty:
        return
    ld_table2=ld_table.copy()
    ld_table2=ld_table2.rename(columns={"RSID1":"RSID2","RSID2":"RSID1"})
    ld=pd.concat([ld_table,ld_table2],sort=True).reset_index(drop=True)
    ld[columns["chrom"]]=ld["RSID2"].apply(lambda x:x.strip("chr").split("_")[0] )
    ld[columns["pos"]]=ld["RSID2"].apply(lambda x:x.strip("chr").split("_")[1] )
    ld[columns["ref"]]=ld["RSID2"].apply(lambda x:x.strip("chr").split("_")[2] )
    ld[columns["alt"]]=ld["RSID2"].apply(lambda x:x.strip("chr").split("_")[3] )
    #filter
    ld=ld.merge(extract_df_1.loc[:,["#variant"] ].rename(columns={"#variant":"RSID1"}),how="inner",on="RSID1")
    ld=map_column(ld,"RSID2_map",columns)
    #is #variant mapped to A strand? maybe, but this should be checked.
    ld=ld.merge(extract_df_2.loc[:,["#variant"] ].rename(columns={"#variant":"RSID2_map"}),how="inner",on="RSID2_map")
    ld.loc[:,"r2"]=ld["correlation"]*ld["correlation"]
    ld=ld.drop(columns=["correlation",columns["chrom"],columns["pos"],columns["ref"],columns["alt"] ])
    #remove temporary files
    corr_files=glob.glob("{}temp_corr.*".format(args.prefix))
    rmcmd="rm {}ld_table.table {}var_lst".format(args.prefix, args.prefix, args.prefix)
    Popen(shlex.split(rmcmd)+corr_files,stderr=subprocess.DEVNULL,stdout=subprocess.DEVNULL)
    return ld

def compare(args):
    """
    Compares our findings to gwascatalog results and/or supplied summary statistic files
    """
    columns={"chrom":args.column_labels[0],"pos":args.column_labels[1],"ref":args.column_labels[2],"alt":args.column_labels[3],"pval":args.column_labels[4]}
    #load original file
    df=pd.read_csv(args.compare_fname,sep="\t")
    if df.empty:
        print("No variants, {}, {} and {} will not be produced".format(args.top_report_out, args.raport_out,args.ld_raport_out))
        return
    necessary_columns=[ columns["chrom"],columns["pos"],columns["ref"],columns["alt"],columns["pval"],"#variant","locus_id","pos_rmin","pos_rmax"]
    df_cols=df.columns.to_list()
    if not all(nec_col in df_cols for nec_col in necessary_columns):
        Exception("GWS variant file {} did not contain all of the necessary columns:\n{} ".format(args.compare_fname,necessary_columns))
    summary_df_1=pd.DataFrame()
    summary_df_2=pd.DataFrame()
    #building summaries, external files and/or gwascatalog
    if args.compare_style in ["file","both"]:
        summary_df_1=load_summary_files(args.summary_fpath,args.endpoints,columns)
    if args.compare_style in ["gwascatalog","both"]:
        gwas_df=None
        gwapi=None
        if os.path.exists("{}gwas_out_mapping.csv".format(args.prefix)) and args.cache_gwas:
            print("reading gwas results from gwas_out_mapping.csv...")
            gwas_df=pd.read_csv("{}gwas_out_mapping.csv".format(args.prefix),sep="\t")
            gwapi=gwcatalog_api.GwasApi()
        else:
            rm_gwas_out="rm {}gwas_out_mapping.csv".format(args.prefix)
            subprocess.call(shlex.split(rm_gwas_out),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            threads=args.gwascatalog_threads
            if args.database_choice=="local":
                gwapi=gwcatalog_api.LocalDB(args.localdb_path)
                threads=1
            elif args.database_choice=="summary_stats":
                gwapi=gwcatalog_api.SummaryApi()
            else:
                gwapi=gwcatalog_api.GwasApi()
            gwas_df=load_api_summaries(df,args.gwascatalog_pad,args.gwascatalog_pval,gwapi,threads,columns)
            if args.cache_gwas:
                gwas_df.to_csv("{}gwas_out_mapping.csv".format(args.prefix),sep="\t",index=False)
        gwas_rename={"chrom":columns["chrom"],"pos":columns["pos"],"ref":columns["ref"],"alt":columns["alt"],"pval":columns["pval"]}
        gwas_df=gwas_df.rename(columns=gwas_rename)
        if gwas_df.empty:
            summary_df_2=gwas_df
        else:    
            #filter out invalid values
            filter_out_codes=[9, 14, 15, 16, 17, 18]
            gwas_df=gwas_df.loc[~gwas_df.loc[:,"code"].isin(filter_out_codes)]
            gwas_df.loc[:,"#variant"]=create_variant_column(gwas_df,chrom=columns["chrom"],pos=columns["pos"],ref=columns["ref"],alt=columns["alt"])
            summary_df_2=gwas_df
            if args.database_choice != "local":
                unique_efos=list(summary_df_2["trait"].unique())
                trait_name_map={}
                for key in unique_efos:
                    trait_name_map[key]=gwapi.get_trait(key)
                summary_df_2.loc[:,"trait_name"]=summary_df_2.loc[:,"trait"].apply(lambda x: trait_name_map[x])
            summary_df_2=summary_df_2.drop_duplicates(subset=["#variant","trait"])
    if args.compare_style not in ["file","gwascatalog","both"]:
        raise NotImplementedError("comparison method '{}' not yet implemented".format(args.compare_style))
    summary_df=pd.concat([summary_df_1,summary_df_2],sort=True)
    #top level df
    if summary_df_2.empty:
        print("no hits from gwascatalog, therefore top report will NOT have any traits shown.")
    top_df=create_top_level_report(df,summary_df_2,args.efo_traits,columns)
    top_df.to_csv(args.top_report_out,sep="\t",index=False)
    if summary_df.empty:
        #just abort, output the top report but no merging summary df cause it doesn't exist
        print("No summary variants, raport will be incomplete")
        raport_out_df=df.copy()
        raport_out_df["variant_hit"]="NA"
        raport_out_df["pval_trait"]="NA"
        raport_out_df["trait"]="NA"
        raport_out_df["trait_name"]="NA"
    else:
        summary_df.to_csv("{}summary_df.csv".format(args.prefix),sep="\t",index=False)
        summary_df=map_column(summary_df,"map_variant",columns)
        df=map_column(df,"map_variant",columns)
        necessary_columns=[columns["pval"],"#variant","map_variant","trait","trait_name"]
        raport_out_df=pd.merge(df,summary_df.loc[:,necessary_columns],how="left",on="map_variant")
        raport_out_df=raport_out_df.drop(columns=["map_variant"])
        raport_out_df=raport_out_df.rename(columns={"#variant_x":"#variant","#variant_y":"#variant_hit","pval_x":columns["pval"],"pval_y":"pval_trait"})
        raport_out_df=raport_out_df.sort_values(by=[columns["chrom"],columns["pos"],columns["ref"],columns["alt"],"#variant"])
    raport_out_df.to_csv(args.raport_out,sep="\t",index=False)
    #Calculate ld between our variants and external variants
    if args.ld_check and (not summary_df.empty):
        #if no groups in base data
        if (("pos_rmin" not in df.columns.to_list()) or ("pos_rmax" not in df.columns.to_list())):
            Exception("ld calculation not supported without grouping. Please supply the flag --group to main.py or gws_fetch.py.") 
        unique_locus_list=df["locus_id"].unique()
        ld_df=pd.DataFrame()
        df.to_csv("{}df.csv".format(args.prefix),index=False,sep="\t")
        #create chromosome list and group loci based on those
        chrom_lst=  sorted([*{*[s.split("_")[0].strip("chr") for s in unique_locus_list]}])
        for chrom in chrom_lst:
            print("------------LD for groups in chromosome {}------------".format(chrom))
            groups=df[df[columns["chrom"]].astype(str) == chrom ].loc[:,"locus_id"].unique()
            plink_cmd="plink --bfile {} --output-chr M --chr {} --make-bed --out {}temp_chrom --memory {}".format( args.ld_panel_path, chrom ,args.prefix,args.plink_mem)
            pr=subprocess.run(shlex.split(plink_cmd),stdout=PIPE,stderr=subprocess.STDOUT)
            if pr.returncode!=0:
                print("PLINK FAILURE for chromosome {}. Error code {}".format(chrom,pr.returncode)  )
                print(pr.stdout)
                continue
            for gr in groups:
                ld=extract_ld_variants(df,summary_df,gr,args,columns)
                if type(ld)==type(None):
                    continue
                ld_df=pd.concat([ld_df,ld],sort=True)
        c5="rm "
        plink_files=glob.glob("{}temp_chrom.*".format(args.prefix))
        Popen(shlex.split(c5)+plink_files,stderr=subprocess.DEVNULL)
        if not ld_df.empty:
            ld_df=ld_df.drop_duplicates(subset=["RSID1","RSID2"],keep="first")
            ld_df=ld_df.merge(summary_df.loc[:,["map_variant","trait","trait_name"]].rename(columns={"map_variant":"RSID2_map"}),how="inner",on="RSID2_map")
            ld_out=df.merge(ld_df,how="inner",left_on="#variant",right_on="RSID1")
            ld_out=ld_out.drop(columns=["RSID1","map_variant","RSID2_map"]).rename(columns={"{}_x".format(columns["pval"]):columns["pval"],"{}_y".format(columns["pval"]):"pval_trait","RSID2":"#variant_hit"})
            ld_out.to_csv(args.ld_raport_out,sep="\t",index=False)
        else:
            print("No variants in ld found, no {} produced.".format(args.ld_raport_out))

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Compare found GWS results to previously found results")
    parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,default="gwascatalog",help="use 'file', 'gwascatalog' or 'both'")
    parser.add_argument("--summary-fpath",dest="summary_fpath",type=str,help="Summary listing file path.")
    parser.add_argument("--endpoint-fpath",dest="endpoints",type=str,help="Endpoint listing file path.")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--plink-memory", dest="plink_mem", type=int, default=12000, help="plink memory for ld clumping, in MB")
    #parser.add_argument("--ld-chromosome-panel-path",dest="ld_chromosome_panel",help="Path to ld panel, where each chromosome is separated. If path is 'path/panel_#chrom.bed', input 'path/panel' ")
    parser.add_argument("--ld-panel-path",dest="ld_panel_path",type=str,help="Filename to the genotype data for ld calculation, without suffix")
    parser.add_argument("--prefix",dest="prefix",type=str,default="",help="output and temporary file prefix. Default value is the base name (no path and no file extensions) of input file. ")
    parser.add_argument("--report-out",dest="raport_out",type=str,default="report_out.csv",help="Report output path")
    parser.add_argument("--ld-report-out",dest="ld_raport_out",type=str,default="ld_report_out.csv",help="LD check report output path")
    parser.add_argument("--gwascatalog-pval",type=str,default="5e-8",help="P-value cutoff for GWASCatalog searches")
    parser.add_argument("--gwascatalog-width-kb",dest="gwascatalog_pad",type=int,default=25,help="gwascatalog range padding")
    parser.add_argument("--gwascatalog-threads",dest="gwascatalog_threads",type=int,default=4,help="Number of concurrent queries to GWAScatalog API. Default 4. Increase if the gwascatalog api takes too long.")
    parser.add_argument("--ldstore-threads",type=int,default=4,help="Number of threads to use with ldstore. Default 4")
    parser.add_argument("--ld-treshold",type=float,default=0.9,help="ld treshold for including ld associations in ld report")
    parser.add_argument("--cache-gwas",action="store_true",help="save gwascatalog results into gwas_out_mapping.csv and load them from there if it exists. Use only for testing.")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval"],help="Names for data file columns. Default is '#chrom pos ref alt pval'.")
    parser.add_argument("--top-report-out",dest="top_report_out",type=str,default="top_report.csv",help="Top level report filename.")
    parser.add_argument("--efo-codes",dest="efo_traits",type=str,nargs="+",default=[],help="Specific EFO codes to look for in the top level report")
    parser.add_argument("--local-gwascatalog",dest='localdb_path',type=str,help="Path to local GWAS Catalog DB.")
    parser.add_argument("--db",dest="database_choice",type=str,default="gwas",help="Database to use for comparison. use 'local','gwas' or 'summary_stats'.")
    args=parser.parse_args()
    if args.prefix!="":
        args.prefix=args.prefix+"."
    args.raport_out = "{}{}".format(args.prefix,args.raport_out)
    args.top_report_out = "{}{}".format(args.prefix,args.top_report_out)
    args.ld_raport_out = "{}{}".format(args.prefix,args.ld_raport_out)
    compare(args)
