#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd 
import numpy as np

def main(args):
    pass

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Annotate results using gnoMAD and additional annotations")
    parser.add_argument("fpath",type=str,help="Filepath of the results to be annotated")
    parser.add_argument("-o","--out-fname",dest="out_fname",type=str,default="out.csv",help="Output filename, default is out.csv")
    args=parser.parse_args()
    fetch_gws(args)