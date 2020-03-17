#!/usr/bin/env python3

import pandas as pd, numpy as np
import json, argparse

def process_phenos(input_array_fname, num_per_worker):
    #read in the data matrix. No header, separator \t
    phenoname="pheno"
    ssname="ss"
    credsetname="credset"
    data=pd.read_csv(input_array_fname,header=None, sep="\t",names=[phenoname,ssname,credsetname])
    data=data.fillna("")
    #create boolean is credset column
    data["is_credset"]=~(data[credsetname]=="")
    #summstat arr, summstat tbi arr, credset arr
    cr_phenos = [getattr(t,phenoname) for t in data.itertuples() if getattr(t,credsetname) != "" ]
    ld_phenos = [getattr(t,phenoname) for t in data.itertuples() if getattr(t,credsetname) == "" ]
    phenolines=[]
    summstatlines=[]
    summstattbilines=[]
    credsetarrlines=[]
    cr_data = data[data[phenoname].isin(cr_phenos)].reset_index(drop=True)
    ld_data = data[data[phenoname].isin(ld_phenos)].reset_index(drop=True)
    for i in range(cr_data.shape[0]//num_per_worker+1):
        idx_start = i*num_per_worker
        idx_end = (i+1)*num_per_worker-1 # -1 because of pandas indexing being start and end inclusive
        tmp_dat = cr_data.loc[idx_start:idx_end,:]
        additional_tabs = len([a for a in list(tmp_dat[phenoname]) if a == ""]) 
        phenos = "\t".join( list(tmp_dat[phenoname]) ) +"".join( ["\t"] * additional_tabs ) + "\n"
        summ_stat = "\t".join( list(tmp_dat[ssname]) ) +"".join( ["\t"] * additional_tabs ) + "\n"
        summ_stat_tb = "\t".join([a+".tbi" for a in tmp_dat[ssname] ]) +"".join( ["\t"] * additional_tabs ) + "\n"
        credset = "\t".join(list(tmp_dat[credsetname])) +"".join( ["\t"] * additional_tabs ) + "\n"
        phenolines.append(phenos)
        summstatlines.append(summ_stat)
        summstattbilines.append(summ_stat_tb)
        credsetarrlines.append(credset)

    for i in range(ld_data.shape[0]//num_per_worker+1):
        idx_start = i*num_per_worker
        idx_end = (i+1)*num_per_worker-1 # -1 because of pandas indexing being start and end inclusive
        tmp_dat = ld_data.loc[idx_start:idx_end,:]
        additional_tabs = len([a for a in list(tmp_dat[phenoname]) if a == ""]) 
        phenos = "\t".join( list(tmp_dat[phenoname]) ) +"".join( ["\t"] * additional_tabs ) + "\n"
        summ_stat = "\t".join( list(tmp_dat[ssname]) ) +"".join( ["\t"] * additional_tabs ) + "\n"
        summ_stat_tb = "\t".join([a+".tbi" for a in tmp_dat[ssname] ]) +"".join( ["\t"] * additional_tabs ) + "\n"
        credset = "\t\n" #use at least one \t so wdl does not remove this cell
        phenolines.append(phenos)
        summstatlines.append(summ_stat)
        summstattbilines.append(summ_stat_tb)
        credsetarrlines.append(credset)
    #write them into files
    with open("pheno_array","w") as f:
        f.writelines(phenolines)
    with open("summ_array","w") as f:
        f.writelines(summstatlines)
    with open("summ_tb_array","w") as f:
        f.writelines(summstattbilines)
    with open("credset_array","w") as f:
        f.writelines(credsetarrlines)

    
if __name__ == "__main__":
    ap = argparse.ArgumentParser("preprocess an input array (pheno,ss,credset) into multiple Array[Array[String]]s")
    ap.add_argument("file",type=str,help="array filename")
    ap.add_argument("--n",type=int,required=True,help="num_workers")
    args = ap.parse_args()
    process_phenos(args.file, args.n)