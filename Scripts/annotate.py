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
    return tb.querys("{}:{}-{}".format(chrom,start,end))

def get_gzip_header(fname):
    """"Returns header for gzipped tsvs, as that is not currently possible using pytabix
    In: file path of gzipped tsv
    Out: header of tsv as a list of column names"""
    #create gzip call
    gzip_call=shlex.split("gzip -cd {}".format(fname))
    head_call=shlex.split("head -n 1")
    gzip_process=Popen(gzip_call,stdout=PIPE)
    head_proc=Popen(head_call,stdin=gzip_process.stdout,stdout=PIPE)
    out=[]
    for line in head_proc.stdout:
        out.append(line.decode().strip().split("\t"))
    return out[0]

def create_variant_column(df,chrom="#chrom",pos="pos",ref="ref",alt="alt"):
    return df.apply( lambda x: "chr{}_{}_{}_{}".format(x[chrom],x[pos],x[ref],x[alt]) ,axis=1)

def annotate(args):
    #load main file
    df=pd.read_csv(args.annotate_fpath,sep="\t")
    original_cols=df.columns.values.tolist()
    #create tabix command
    tbxlst=[]
    tb=tabix.open(args.gnomad_path)
    for _,row in df.iterrows():
                tbxlst=tbxlst+list(pytabix(tb,row["#chrom"],int(row["pos"]),int(row["pos"]) ) )
    tbxheader=get_gzip_header(args.gnomad_path)
    annotation_tbx=pd.DataFrame(tbxlst,columns=tbxheader )
    annotation_tbx[annotation_tbx.columns]=annotation_tbx[annotation_tbx.columns].apply(pd.to_numeric,errors="ignore")
    annotation_tbx=annotation_tbx.drop_duplicates(subset=["#CHROM","POS","REF","ALT"]).rename(columns={"#CHROM":"#chrom","POS":"pos","REF":"ref","ALT":"alt"})
    #add gnomad annotations, by fetching the correct rows using tabix and calculating the necessary features
    allele_frequency_groups=["AF_fin","AF_nfe","AF_nfe_est","AF_nfe_nwe","AF_nfe_onf","AF_nfe_seu"]
    #join based on chrom, pos, ref, and alt
    df=df.merge(annotation_tbx,on=["#chrom","pos","ref","alt"],how="left")
    all_counts_nfe=["AC_nfe","AC_nfe_est","AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
    all_number_nfe=["AN_nfe","AN_nfe_est","AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
    df.loc[:,"AC_sum_nfe"]=df[all_counts_nfe].sum(axis=1)
    df.loc[:,"AN_sum_nfe"]=df[all_number_nfe].sum(axis=1)
    df.loc[:,"FI_enrichment_nfe"]=df.loc[:,"AN_sum_nfe"]*df.loc[:,"AF_fin"]/df.loc[:,"AC_sum_nfe"]
    df.loc[:,"FI_enrichment_nfe"]=df.loc[:,"FI_enrichment_nfe"].clip(lower=0.0,upper=1e6)#replace inf with 1 000 000
    #same for non-estonian europeans 
    all_counts_nfe=["AC_nfe","AC_nfe_nwe","AC_nfe_onf","AC_nfe_seu"]
    all_number_nfe=["AN_nfe","AN_nfe_nwe","AN_nfe_onf","AN_nfe_seu"]
    df.loc[:,"AC_sum_nfe_est"]=df[all_counts_nfe].sum(axis=1)
    df.loc[:,"AN_sum_nfe_est"]=df[all_number_nfe].sum(axis=1)
    df.loc[:,"FI_enrichment_nfe_no_est"]=df.loc[:,"AN_sum_nfe_est"]*df.loc[:,"AF_fin"]/df.loc[:,"AC_sum_nfe_est"]
    df.loc[:,"FI_enrichment_nfe_no_est"]=df.loc[:,"FI_enrichment_nfe_no_est"].clip(0.0,1e6)#replace inf with 1 000 000
    #TODO: add enrichment for non-finnish,estonian,swedist europeans
    #NOTE: the swedish counts o not exist in gnomad.genomes.r2.1.sites.liftover.b38.finngen.r2pos.af.ac.an.tsv.gz 
    #add finngen annotations, using tabix
    tbxlst=[]
    tb=tabix.open(args.finngen_path)
    for _,row in df.iterrows():
        tbxlst=tbxlst+list(pytabix(tb,"chr{}".format(row["#chrom"]),int(row["pos"]),int(row["pos"]) ) )
    tbxheader=get_gzip_header(args.finngen_path)
    fg_df=pd.DataFrame(tbxlst,columns=tbxheader )
    fg_df[fg_df.columns]=fg_df[fg_df.columns].apply(pd.to_numeric,errors="ignore")
    df.loc[:,"#variant"]=create_variant_column(df)
    fg_df=fg_df.drop_duplicates(subset=["#variant"])
    #add #variant column in df
    #join based on variant
    fg_cols=fg_df.columns.values.tolist()
    colist=list(set(fg_cols)-set(["chrom","pos","ref","alt"]))
    fg_df=fg_df.loc[:,]
    df=df.merge(fg_df.loc[:,colist],on="#variant",how="left")
    fg_out=["most_severe","INFO","gene"]
    
    #get batch columns, for INFO, IMP, AF
    info=list(filter(lambda s: "INFO_" in s,fg_cols))
    imp=list(filter(lambda s: "IMP_" in s,fg_cols))
    af=list(filter(lambda s: "AF_" in s,fg_cols))
    fg_batches=af+info+imp
    
    out_columns=original_cols+allele_frequency_groups+["FI_enrichment_nfe","FI_enrichment_nfe_no_est"]+fg_out
    if args.batch_freq:
        out_columns+=fg_batches

    df=df.loc[:,out_columns]
    df=df.rename(index=str,columns={"gene":"most_severe_gene"})
    df.to_csv(path_or_buf=args.out_fname,sep="\t",index=False)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Annotate results using gnoMAD and additional annotations")
    parser.add_argument("annotate_fpath",type=str,help="Filepath of the results to be annotated")
    parser.add_argument("-g","--gnomad-path",dest="gnomad_path",type=str,help="Gnomad annotation file filepath")
    parser.add_argument("--include-batch-freq",dest="batch_freq",action="store_true",help="Include batch frequencies from finngen annotations")
    parser.add_argument("--finngen-path",dest="finngen_path",type=str,default=None,help="Finngen annotation file filepath")
    parser.add_argument("-o","--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    args=parser.parse_args()
    annotate(args)