import argparse,shlex,subprocess, glob, re
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
from autoreporting_utils import *
import os
from data_access import datafactory
from typing import Dict, List

def top_report_lead_cols(top_df: pd.DataFrame, report_df: pd.DataFrame, columns: Dict[str,str], variant_col: str, lead_cols: List[str]) -> pd.DataFrame:
    """Returns lead variant columns
    Args:
        top_df: Dataframe
        report_df: Full dataframe
        columns: column names
        variant_col (str): Variant col name
        lead_cols (List[str]) : Columns that are added to the top report
        TODO
    Returns:
        (pd.DataFrame): dataframe with lead variant columns
    """
    cols=[variant_col]+lead_cols
    annotation_df = report_df[cols]
    out = pd.merge(top_df[variant_col], annotation_df, how="left",on=variant_col).drop_duplicates(subset=[variant_col],keep="first")
    return out


def create_top_level_report(report_df,efo_traits,columns,grouping_method,significance_threshold,strict_ld_threshold, extra_cols):
    """
    Create a top level report from which it is easy to see which loci are novel
    In: report_out df, traits that appear in matching_pheno_gwas_catalog_hits, column names
    Out: Dataframe with a row for every lead variant in df, with columns locus_id chr, start, end, matching_pheno_gwas_catalog_hits other_gwas_hits  
    """ 

    pheno_cols = ["phenotype", "longname", "n_cases", "n_controls"]
    pheno_col_rename = {
                        "longname":"phenotype",
                        "phenotype":"phenotype_abbreviation",
                        "n_cases":"Cases",
                        "n_controls":"Controls"
    }

    group_cols = ["locus_id", columns["chrom"], columns["pos"], columns["ref"], columns["alt"], columns["pval"]]
    group_cols_rename = {columns["chrom"]: "chrom",
                         columns["pos"]: "pos",
                         columns["ref"]: "ref",
                         columns["alt"]: "alt",
                         columns["pval"]: "pval"}

    lead_var_cols = ["beta_previous_release","pval_previous_release",
                     "most_severe_consequence","most_severe_gene",
                     "GENOME_FI_enrichment_nfe_est"]+extra_cols
    lead_col_d: Dict = {k:"lead_{}".format(k) for k in lead_var_cols}
    lead_col_d["GENOME_FI_enrichment_nfe_est"]="lead_enrichment"

    gnomad_add_cols = ["functional_category","enrichment_nfsee","fin.AF","fin.AN","fin.AC", "fin.homozygote_count", "fet_nfsee.odds_ratio", "fet_nfsee.p_value", "nfsee.AC", "nfsee.AN", "nfsee.AF", "nfsee.homozygote_count"]
    gnomad_add_cols_rename = {k:"gnomAD_{}".format(k) for k in gnomad_add_cols}
    cs_cols = ["cs_id", "cs_size", "cs_log10bf", "cs_number", "cs_region","good_cs"]
    cs_cols_rename = {"cs_log10bf":"cs_log_bayes_factor" }

    aggregated_cols = ["credible_set_min_r2_value","start", "end", "found_associations_strict", "found_associations_relaxed",
                       "credible_set_variants", "functional_variants_strict",
                       "functional_variants_relaxed", "specific_efo_trait_associations_strict",
                       "specific_efo_trait_associations_relaxed","n_ld_partners_0_8"]

    lead_cols=list(lead_col_d.values())
    gnomad_cols = list(gnomad_add_cols_rename.values())

    top_level_cols=["phenotype",
                    "phenotype_abbreviation",
                    "locus_id",
                    "Cases",
                    "Controls",
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                    "pval"]+\
                    lead_cols+\
                    gnomad_cols+\
                    ["cs_id",
                     "cs_size",
                     "cs_log_bayes_factor",
                     "cs_number",
                     "cs_region",
                     "good_cs"]+\
                    aggregated_cols

    df=report_df.copy()

    if df.empty:
        return pd.DataFrame(columns=top_level_cols)

    list_of_loci=list(df["locus_id"].unique())

    lead_var_df = df.loc[df["locus_id"]==df["#variant"],["#variant"]].drop_duplicates()
    lead_var_idx = lead_var_df.index

    #get the pheno_cols dataframe
    try:
        pheno_info_df = df.loc[lead_var_idx,pheno_cols+["#variant"]].rename(columns=pheno_col_rename)
    except:
        pheno_info_df = pd.DataFrame(columns=["phenotype","phenotype_abbreviation","Cases","Controls","#variant"])
    #get the group cols dataframe
    try:
        group_info_df = df.loc[lead_var_idx,group_cols+["#variant"]]
        group_info_df=group_info_df.rename(columns=group_cols_rename)
    except:
        group_info_df = pd.DataFrame(columns = ["#variant","locus_id","chrom","pos","ref","alt"])
    
    #get cs cols dataframe
    try:
        cs_info_df = df.loc[lead_var_idx,cs_cols+["#variant"]]
        cs_info_df=cs_info_df.rename(columns=cs_cols_rename)
    except:
        cs_info_df= pd.DataFrame(columns=cs_cols+["#variant"]).rename(columns=cs_cols_rename)
    #get lead cols dataframe
    try:
        lead_info_df = df.loc[lead_var_idx,lead_var_cols+["#variant"]]
        lead_info_df=lead_info_df.rename(columns= lead_col_d)
    except:
        lead_info_df = pd.DataFrame(columns = lead_cols+["#variant"])
    #get gnomad cols
    try:
        gnomad_info_df = df.loc[lead_var_idx,gnomad_add_cols+["#variant"]]
        gnomad_info_df = gnomad_info_df.rename(columns=gnomad_add_cols_rename)
    except:
        gnomad_info_df=pd.DataFrame(columns=gnomad_cols+["#variant"])

    top_level_df=pd.DataFrame(columns=["#variant"]+aggregated_cols)

    
    for locus_id in list_of_loci:
        # The row is a dict which will contain the aggregated values for a single group
        row = {}
        # Separate group variants from all variants
        loc_variants=df.loc[df["locus_id"]==locus_id,:]

        # Create strict group. The definition changes based on the grouping method.
        try:
            locus_cs_id = loc_variants.loc[loc_variants["#variant"] == locus_id,"cs_id"].values[0]
        except:
            raise Exception("lead variant not in group, results incorrect.")
        strict_group=None
        if grouping_method == "cred":
            #strict group = variants in the credible set of this locus
            strict_group = loc_variants[loc_variants["cs_id"]==locus_cs_id].copy()
        elif grouping_method == "ld":
            #strict group = variants that have low enough pvalue and high enough correlation with lead variant
            strict_group = loc_variants[(loc_variants[columns["pval"]]<=significance_threshold) & (loc_variants["r2_to_lead"]>=strict_ld_threshold )].copy()
        else:
            #variants that have low enough pvalue
            strict_group = loc_variants[loc_variants[columns["pval"]]<=significance_threshold].copy()
        
        row['#variant']=locus_id
        row["start"]=int(np.amin(loc_variants[columns["pos"]]))
        row["end"]=int(np.amax(loc_variants[columns["pos"]]))
        credset_vars = strict_group.loc[strict_group["cs_id"]==locus_cs_id,["#variant","cs_prob","r2_to_lead"]].drop_duplicates()
        cred_set=";".join( "{}|{:.3g}|{:.3g}".format(t._1,t.cs_prob,t.r2_to_lead) for t in  credset_vars.itertuples() )
        # Get credible set variants in relazed & strict group, as well as functional variants. 
        # Try because it is possible that functional data was skipped.
        try:
            func_s = loc_variants.loc[~loc_variants["functional_category"].isna(),["#variant","functional_category","r2_to_lead"] ].drop_duplicates()
            func_set=";".join("{}|{}|{:.3g}".format(t._1,t.functional_category,t.r2_to_lead) for t in  func_s.itertuples())
            func_s_strict = strict_group.loc[~strict_group["functional_category"].isna(),["#variant","functional_category","r2_to_lead"] ].drop_duplicates()
            func_set_strict=";".join("{}|{}|{:.3g}".format(t._1,t.functional_category,t.r2_to_lead) for t in  func_s_strict.itertuples())
        except:
            func_set=""
            func_set_strict=""
        row["credible_set_variants"]=cred_set
        row["functional_variants_strict"]=func_set_strict
        row["functional_variants_relaxed"]=func_set
        # Solve traits for relaxed and strict groups.
        all_traits=(
                    loc_variants[["trait","trait_name","r2_to_lead"]]
                    .sort_values(by="r2_to_lead",ascending=False)
                    .drop_duplicates(subset=["trait","trait_name"], keep="first")
                    .dropna()
                    )
        strict_traits=(
                    strict_group[["trait","trait_name","r2_to_lead"]]
                    .sort_values(by="r2_to_lead",ascending=False)
                    .drop_duplicates(subset=["trait","trait_name"],keep="first")
                    .dropna()
                    )
        
        other_traits_relaxed = all_traits[ ~all_traits["trait"].isin(efo_traits) ].copy()
        other_traits_strict = strict_traits[ ~strict_traits["trait"].isin(efo_traits) ].copy()

        row["found_associations_relaxed"]=";".join( "{}|{:.3g}".format(t.trait_name,t.r2_to_lead) for t in other_traits_relaxed.itertuples() )
        row["found_associations_strict"]=";".join( "{}|{:.3g}".format(t.trait_name,t.r2_to_lead) for t in other_traits_strict.itertuples() )
        
        matching_traits_relaxed= all_traits[ all_traits["trait"].isin(efo_traits) ].copy()
        matching_traits_strict = strict_traits[ strict_traits["trait"].isin(efo_traits) ].copy()
        row["specific_efo_trait_associations_relaxed"]=";".join( "{}|{:.3g}".format(t.trait_name,t.r2_to_lead) for t in matching_traits_relaxed.itertuples() )
        row["specific_efo_trait_associations_strict"]=";".join( "{}|{:.3g}".format(t.trait_name,t.r2_to_lead) for t in matching_traits_strict.itertuples() )
        try:
            row["credible_set_min_r2_value"] = np.nanmin(credset_vars.loc[:, "r2_to_lead" ].values)
        except:
            row["credible_set_min_r2_value"] = np.nan

        #N ld partners with LD >0.8
        r2_thresh_0_8 = 0.8
        n_ld_gt_0_8 = loc_variants[((loc_variants["r2_to_lead"]>r2_thresh_0_8) & (loc_variants["cs_id"]!=locus_cs_id)),"r2_to_lead"].count()
        row["n_ld_partners_0_8"] = n_ld_gt_0_8
        top_level_df=top_level_df.append(row,ignore_index=True)

    #merge the different dataframes to top_level_df
    top_level_df=top_level_df.merge(pheno_info_df,on="#variant",how="left")
    top_level_df=top_level_df.merge(group_info_df,on="#variant",how="left")
    top_level_df=top_level_df.merge(cs_info_df,on="#variant",how="left")
    top_level_df=top_level_df.merge(lead_info_df,on="#variant",how="left")
    top_level_df = top_level_df.merge(gnomad_info_df,on="#variant",how="left")
    return top_level_df[top_level_cols]


if __name__ == "__main__":
    #read parameters
    parser=argparse.ArgumentParser(description="Create top report")
    parser.add_argument("report_fname",type=str,help="autoreporting report file")
    parser.add_argument("--sign-treshold",dest="sig_treshold",type=float,help="Signifigance treshold",default=5e-8)
    parser.add_argument("--grouping-method",dest="grouping_method",type=str,default="simple",help="Decide grouping method, simple, ld or cred, default simple")
    parser.add_argument("--column-labels",dest="column_labels",metavar=("CHROM","POS","REF","ALT","PVAL"),nargs=5,default=["#chrom","pos","ref","alt","pval"],help="Names for data file columns. Default is '#chrom pos ref alt pval'.")
    parser.add_argument("--extra-cols",dest="extra_cols",nargs="*",default=[],help="extra columns in the summary statistic you want to add to the results")
    parser.add_argument("--efo-codes",dest="efo_traits",type=str,nargs="+",default=[],help="Specific EFO codes to look for in the top level report")
    parser.add_argument("--strict-group-r2",dest="strict_group_r2",type=float,default=0.5,help="R^2 threshold for including variants in strict groups in top report")
    parser.add_argument("--top-report-out",dest="top_report_out",type=str,default="top_report.tsv",help="Top level report filename.")
    args=parser.parse_args()
    #read data
    report_df = pd.read_csv(args.report_fname,sep="\t")
    print("read file")
    #create top report
    columns = columns_from_arguments(args.column_labels)
    top_df=create_top_level_report(report_df,efo_traits=args.efo_traits,columns=columns,grouping_method= args.grouping_method,significance_threshold=args.sig_treshold,strict_ld_threshold=args.strict_group_r2, extra_cols=args.extra_cols)
    print("made top report, shape:",top_df.shape)
    #write to file
    top_df.fillna("NA").replace("","NA").to_csv(args.top_report_out,sep="\t",index=False,float_format="%.3g")