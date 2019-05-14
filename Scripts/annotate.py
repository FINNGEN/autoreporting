#! /usr/bin/python3

import argparse,shlex,subprocess
from subprocess import Popen, PIPE
import pandas as pd 
import numpy as np
import tabix
#TODO: make a system for making sure we can calculate all necessary fields,
#e.g by checking that the columns exist
def tabix_command(df,path,prefix=""):
    tcall="tabix "+path+" -h " 
    for t in df.itertuples():
        tcall=tcall+" "+str(prefix)+str(t._1)+":"+str(t.pos)+"-"+str(t.pos)+" "
    return tcall

def pytabix(tb,chrom,start,end):
    try:
        retval=tb.querys("{}:{}-{}".format(chrom,start,end))
        return retval
    except tabix.TabixError:
        return []

def load_tb_df(df,tb,fpath,chrom_prefix=""):
    tbxlst=[]
    for _,row in df.iterrows():
        tbxlst=tbxlst+list(pytabix(tb,"{}{}".format(chrom_prefix,row["#chrom"]),int(row["pos"]),int(row["pos"]) ) )
    header=get_gzip_header(fpath)
    out_df=pd.DataFrame(tbxlst,columns=header )
    out_df[out_df.columns]=out_df[out_df.columns].apply(pd.to_numeric,errors="ignore")
    return out_df

def get_gzip_header(fname):
    """"Returns header for gzipped tsvs, as that is not currently possible using pytabix
    In: file path of gzipped tsv
    Out: header of tsv as a list of column names"""
    #create gzip call
    gzip_call=shlex.split("gzip -cd {}".format(fname))
    head_call=shlex.split("head -n 1")
    out=[]
    with Popen(gzip_call,stdout=PIPE) as gz:
        with Popen(head_call,stdin=gz.stdout,stdout=PIPE) as hd:
            for line in hd.stdout:
                out.append(line.decode().strip().split("\t"))
    return out[0]

def create_variant_column(df,chrom="#chrom",pos="pos",ref="ref",alt="alt"):
    return df.apply( lambda x: "chr{}_{}_{}_{}".format(x[chrom],x[pos],x[ref],x[alt]) ,axis=1)

def calculate_enrichment(gnomad_df,fi_af_col,count_nfe_lst,number_nfe_lst):
    """Calculate enrichment for finns vs the other group, 
    which is defined in the column name lists
    In: gnomad dataframe, finnish allele frequency column,
    other group allele count and number columns
    Out: enrichment column"""
    nfe_counts=gnomad_df.loc[:,count_nfe_lst].sum(axis=1)
    nfe_numbers=gnomad_df.loc[:,number_nfe_lst].sum(axis=1)
    finn_freq=gnomad_df.loc[:,fi_af_col]
    enrichment=nfe_numbers*finn_freq/nfe_counts
    # clip enrichment
    enrichment=enrichment.clip(lower=0.0,upper=1e6)
    return enrichment

def create_rename_dict(list_of_names, prefix):
    d={}
    for value in list_of_names:
        d[value]="{}{}".format(prefix,value)
    return d

def annotate(args):
    #TODO
    """Annotation function info
    What it does
    In:
    Out:"""
    
    #columns that we want to take from gnomad and finngen annotations
    gnomad_gen_cols=["AF_fin","AF_nfe","AF_nfe_est","AF_nfe_nwe","AF_nfe_onf","AF_nfe_seu","FI_enrichment_nfe","FI_enrichment_nfe_est"]
    gnomad_exo_cols=["AF_nfe_bgr","AF_fin","AF_nfe","AF_nfe_est","AF_nfe_swe","AF_nfe_nwe","AF_nfe_onf",\
        "AF_nfe_seu","FI_enrichment_nfe","FI_enrichment_nfe_est","FI_enrichment_nfe_swe","FI_enrichment_nfe_est_swe"]
    finngen_cols=["most_severe_gene","most_severe_consequence"]

    #load main file
    df=pd.read_csv(args.annotate_fpath,sep="\t")
    original_cols=df.columns.values.tolist()
    #load gnomad_genomes
    tb_g=tabix.open(args.gnomad_genome_path)
    gnomad_genomes=load_tb_df(df,tb_g,args.gnomad_genome_path)

    #load gnomad_exomes
    tb_e=tabix.open(args.gnomad_exome_path)
    gnomad_exomes=load_tb_df(df,tb_e,args.gnomad_exome_path)

    #load finngen annotations
    tb_f=tabix.open(args.finngen_path)
    fg_df=load_tb_df(df,tb_f,args.finngen_path,chrom_prefix="chr")
    fg_df=fg_df.drop_duplicates(subset=["#variant"])

    if not gnomad_genomes.empty:
        gnomad_genomes=gnomad_genomes.drop_duplicates(subset=["#CHROM","POS","REF","ALT"]).rename(columns={"#CHROM":"#chrom","POS":"pos","REF":"ref","ALT":"alt"})
        gnomad_genomes.loc[:,"#variant"]=create_variant_column(gnomad_genomes)
        #calculate enrichment for gnomad genomes, nfe, nfe without est
        gn_gen_nfe_counts=["AC_nfe","AC_nfe_est","AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
        gn_gen_nfe_nums=["AN_nfe","AN_nfe_est","AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
        gn_gen_nfe_est_counts=["AC_nfe","AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
        gn_gen_nfe_est_nums=["AN_nfe","AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
        #calculate enrichment
        
        gnomad_genomes.loc[:,"FI_enrichment_nfe"]=calculate_enrichment(gnomad_genomes,"AF_fin",gn_gen_nfe_counts,gn_gen_nfe_nums)
        gnomad_genomes.loc[:,"FI_enrichment_nfe_est"]=calculate_enrichment(gnomad_genomes,"AF_fin",gn_gen_nfe_est_counts,gn_gen_nfe_est_nums)
    else:
        for val in ["FI_enrichment_nfe","FI_enrichment_nfe_est","#variant"]:
            gnomad_genomes[val]=None
    
    if not gnomad_exomes.empty:
        gnomad_exomes=gnomad_exomes.drop_duplicates(subset=["#CHROM","POS","REF","ALT"]).rename(columns={"#CHROM":"#chrom","POS":"pos","REF":"ref","ALT":"alt"})
        gnomad_exomes.loc[:,"#variant"]=create_variant_column(gnomad_exomes)
        #calculate enrichment for gnomax exomes, nfe, nfe without est, nfe without swe, nfe without est, swe?
        gn_exo_nfe_counts=["AC_nfe","AC_nfe_bgr","AC_nfe_est","AC_nfe_onf","AC_nfe_seu","AC_nfe_swe"]
        gn_exo_nfe_nums=["AN_nfe","AN_nfe_bgr","AN_nfe_est","AN_nfe_onf","AN_nfe_seu","AN_nfe_swe"]
        gn_exo_nfe_est_counts=["AC_nfe","AC_nfe_bgr","AC_nfe_onf","AC_nfe_seu","AC_nfe_swe"]
        gn_exo_nfe_est_nums=["AN_nfe","AN_nfe_bgr","AN_nfe_onf","AN_nfe_seu","AN_nfe_swe"]
        gn_exo_nfe_swe_counts=["AC_nfe","AC_nfe_bgr","AC_nfe_est","AC_nfe_onf","AC_nfe_seu"]
        gn_exo_nfe_swe_nums=["AN_nfe","AN_nfe_bgr","AN_nfe_est","AN_nfe_onf","AN_nfe_seu"]
        gn_exo_nfe_est_swe_counts=["AC_nfe","AC_nfe_bgr","AC_nfe_onf","AC_nfe_seu"]
        gn_exo_nfe_est_swe_nums=["AN_nfe","AN_nfe_bgr","AN_nfe_onf","AN_nfe_seu"]
    
        gnomad_exomes.loc[:,"FI_enrichment_nfe"]=calculate_enrichment(gnomad_exomes,"AF_fin",gn_exo_nfe_counts,gn_exo_nfe_nums)
        gnomad_exomes.loc[:,"FI_enrichment_nfe_est"]=calculate_enrichment(gnomad_exomes,"AF_fin",gn_exo_nfe_est_counts,gn_exo_nfe_est_nums)
        gnomad_exomes.loc[:,"FI_enrichment_nfe_swe"]=calculate_enrichment(gnomad_exomes,"AF_fin",gn_exo_nfe_swe_counts,gn_exo_nfe_swe_nums)
        gnomad_exomes.loc[:,"FI_enrichment_nfe_est_swe"]=calculate_enrichment(gnomad_exomes,"AF_fin",gn_exo_nfe_est_swe_counts,gn_exo_nfe_est_swe_nums)
    else:
        for val in ["FI_enrichment_nfe","FI_enrichment_nfe_est","FI_enrichment_nfe_swe","FI_enrichment_nfe_est_swe","#variant"]:
            gnomad_exomes[val]=None
    
        
        
    
    #rename gnomad_exomes and gnomad_genomes
    gnomad_genomes=gnomad_genomes.loc[:,["#variant"]+gnomad_gen_cols]
    gn_gen_rename_d=create_rename_dict(gnomad_gen_cols,"GENOME_")
    gnomad_genomes=gnomad_genomes.rename(columns=gn_gen_rename_d)
    gnomad_exomes=gnomad_exomes.loc[:,["#variant"]+gnomad_exo_cols]

    gn_exo_rename_d=create_rename_dict(gnomad_exo_cols,"EXOME_")
    gnomad_exomes=gnomad_exomes.rename(columns=gn_exo_rename_d)
    
    fg_df=fg_df.rename(columns={"gene":"most_severe_gene","most_severe":"most_severe_consequence"})
    fg_cols=fg_df.columns.values.tolist()
    info=list(filter(lambda s: "INFO_" in s,fg_cols))
    imp=list(filter(lambda s: "IMP_" in s,fg_cols))
    af=list(filter(lambda s: "AF_" in s,fg_cols))

    fg_batch_lst=info+imp+af
    fg_batch_rename=create_rename_dict(fg_batch_lst,"FG_")
    fg_df=fg_df.rename(columns=fg_batch_rename)
    fg_batch_lst=["FG_{}".format(x) for x in fg_batch_lst]
    if args.batch_freq:
        finngen_cols=finngen_cols+fg_batch_lst
    fg_df=fg_df.loc[:,["#variant"]+finngen_cols]

    #merge the wanted columns into df
    df=df.merge(gnomad_genomes,how="left",on="#variant")
    df=df.merge(gnomad_exomes,how="left",on="#variant")
    
    df=df.merge(fg_df,how="left",on="#variant")

    df.to_csv(path_or_buf=args.out_fname,sep="\t",na_rep="NA",index=False)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Annotate results using gnoMAD and additional annotations")
    parser.add_argument("annotate_fpath",type=str,help="Filepath of the results to be annotated")
    parser.add_argument("--gnomad-genome-path",dest="gnomad_genome_path",type=str,help="Gnomad genome annotation file filepath")
    parser.add_argument("--gnomad-exome-path",dest="gnomad_exome_path",type=str,help="Gnomad exome annotation file filepath")
    parser.add_argument("--include-batch-freq",dest="batch_freq",action="store_true",help="Include batch frequencies from finngen annotations")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default=None,help="Finngen annotation file filepath")
    parser.add_argument("-o","--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    args=parser.parse_args()
    annotate(args)