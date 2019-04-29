#! /usr/bin/python3

import argparse,shlex,subprocess
import pandas as pd
import numpy as np


def main(args):
    #if using summary file
    if args.compare_style=="file":
        #load original file

        #load summary file
    else:
        raise NotImplementedError("comparison method '{}' not yet implemented".format(args.compare_style))

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Compare found GWS results to previously found results")
    parser.add_argument("fname",type=str,help="GWS result file")
    parser.add_argument("--compare-style",type=str,help="use database or file")
    parser.add_argument("--summary-fpath",type=str,help="comparison summary filepath")
    
    args=parser.parse_args()
    main(args)