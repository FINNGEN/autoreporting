#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd
import numpy as np
from gws_fetch import create_variant_column


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
                    

    else:
        raise NotImplementedError("comparison method '{}' not yet implemented".format(args.compare_style))

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Compare found GWS results to previously found results")
    parser.add_argument("compare_fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,help="use database or file")
    parser.add_argument("--summary-fpath",dest="summary_files",metavar="FILE",nargs="+",help="comparison summary filepaths")
    parser.add_argument("--endpoints",type=str,nargs="+",help="biological endpoint, as many as summaries")
    parser.add_argument("--build-38",dest="build_38",action="store_true",help="Whether is in GRCh38")
    parser.add_argument("--check-for-ld",dest="ld_check",action="store_true",help="Whether to check for ld between the summary statistics and GWS results")
    parser.add_argument("--raport-out",dest="raport_out",type=str,default="raport_output.csv",help="Raport output path")
    args=parser.parse_args()
    compare(args)