#! /usr/bin/env python3

import pandas as pd, numpy as np
import json, argparse, re, datetime

def write_readme(filelist_fname,output_name,release_number,input_array_fname,readme_template_fname):
    #read in file list
    with open(filelist_fname,"r") as f:
        flist=f.readlines()
    in_arr = pd.read_csv(input_array_fname,sep="\t",names=["pheno","ss","cr"])
    num_cr=np.sum(in_arr["cr"].notna())
    num_phenos=in_arr.shape[0]
    reg_group = "\.top\.out"
    reg_report = "\.report\.out"
    group_files = [a for a in flist if re.search(reg_group,a)]
    report_files = [a for a in flist if re.search(reg_report,a)]
    #calculate/get things to replace in readme
    #number of group reports
    group_num=len(group_files)
    report_num=len(report_files)
    today=str(datetime.date.today())
    replace_dict={}
    replace_dict["{RELEASE_NUMBER}"] = str(release_number)
    replace_dict["{DATE}"] = today
    replace_dict["{RESULT_PHENO_NUM}"] = str(report_num)
    replace_dict["{CR_PHENO_NUM}"] = str(num_cr)
    replace_dict["{PHENO_NUM}"] = str(num_phenos)
    #replace keywords in readme
    with open(readme_template_fname,"r") as f:
        readme_lines = f.readlines()
    readme_output = []
    for readme_line in readme_lines:
        scratchline = readme_line
        for key,value in replace_dict.items():
            scratchline = re.sub(key,value,scratchline)
        readme_output.append(scratchline)
    #output readme
    with open(output_name,"w") as f:
        f.writelines(readme_output)

if __name__ == "__main__":
    ap = argparse.ArgumentParser("Process Autoreporting WDL results into a folder with correct README and stuff")
    ap.add_argument("--release",type=str, required=True,help="Release number to include in README")
    ap.add_argument("input_files",type=str,help="Filelist")
    ap.add_argument("--in-array",required=True,type=str,help="input array")
    ap.add_argument("--out",required=True,type=str,help="Output filename")
    ap.add_argument("--readme-template",required=True,type=str,help="README template filename")
    args = ap.parse_args()
    write_readme(args.input_files, args.out,args.release,args.in_array,args.readme_template)