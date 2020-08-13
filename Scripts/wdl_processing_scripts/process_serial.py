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

    phenolines=[]
    summstatlines=[]
    summstattbilines=[]
    credsetarrlines=[]
    #num of lines is data.shape[0]//num_per_worker plus one if they don't divide evenly
    line_amount = data.shape[0]//num_per_worker
    if data.shape[0]%num_per_worker:
        line_amount = line_amount+1

    for i in range(line_amount):
        idx_start = i*num_per_worker
        idx_end = (i+1)*num_per_worker-1 # -1 because of pandas indexing being start and end inclusive
        tmp_dat = data.loc[idx_start:idx_end,:]
        
        phenos = "\t".join( list(tmp_dat[phenoname]) ) + "\n"
        summ_stat = "\t".join( list(tmp_dat[ssname]) ) + "\n"
        summ_stat_tb = "\t".join([a+".tbi" for a in tmp_dat[ssname] ]) + "\n"
        credset = "\t".join(list(tmp_dat[credsetname])) + "\n"
        phenolines.append(phenos)
        summstatlines.append(summ_stat)
        summstattbilines.append(summ_stat_tb)
        credsetarrlines.append(credset)

    #strip empty lines
    phenolines=[a for a in phenolines if a not in ["","\n"]]
    summstatlines=[a  for a in summstatlines if a not in ["","\n"]]
    summstattbilines=[a for a in summstattbilines if a not in ["","\n"]]
    credsetarrlines=[a for a in credsetarrlines if a not in ["","\n"]]
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