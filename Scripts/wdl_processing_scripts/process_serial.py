import pandas as pd, numpy as np
import json, argparse

def process_phenos(input_array_fname,num_per_worker):
    #read in the data matrix. No header, separator \t
    data=pd.read_csv(input_array_fname,header=None, sep="\t",names=["pheno","ss","credset"])
    data=data.fillna("")
    #create boolean is credset column
    data["is_credset"]=~(data["credset"]=="")
    #boolmap
    bool_map = {row["pheno"]:row["is_credset"] for (idx,row) in data.iterrows()}
    #summstat arr, summstat tbi arr, credset arr
    phenolines=[]
    summstatlines=[]
    summstattbilines=[]
    credsetarrlines=[]
    for i in range(data.shape[0]//num_per_worker+1):
        idx_start = i*num_per_worker
        idx_end = (i+1)*num_per_worker-1 # -1 because of pandas indexing being start and end inclusive
        tmp_dat = data.loc[idx_start:idx_end,:]
        phenos = "\t".join( list(tmp_dat["pheno"]) ) + "\n"
        summ_stat = "\t".join( list(tmp_dat["ss"]) ) + "\n"
        summ_stat_tb = "\t".join([a+".tbi" for a in tmp_dat["ss"] ]) + "\n"
        credset_grp = [a for a in list(tmp_dat["credset"]) if a != "" ] 
        credset = "\t".join(credset_grp) + "\n"
        phenolines.append(phenos)
        summstatlines.append(summ_stat)
        summstattbilines.append(summ_stat_tb)
        credsetarrlines.append(credset)
    #write them into files
    with open("boolmap","w") as f:
        f.write(json.dumps(bool_map))
    with open("pheno_array","w") as f:
        f.writelines(phenolines)
    with open("summ_array","w") as f:
        f.writelines(summstatlines)
    with open("summ_tb_array","w") as f:
        f.writelines(summstattbilines)
    with open("credset_array","w") as f:
        f.writelines(credsetarrlines)

    
if __name__ == "__main__":
    ap = argparse.ArgumentParser("preprocess an input array (pheno,ss,credset) into multiple Array[Array[String]]s and a Map[String,Boolean]")
    ap.add_argument("file",type=str,help="array filename")
    ap.add_argument("--n",type=int,required=True,help="num_workers")
    args = ap.parse_args()
    process_phenos(args.file, args.n)