"""Script for concatenating autoreporting reports into one file."""

import os
import glob
import csv
import argparse
import re
import pandas as pd


def concat_files(folderpath: str, output_fname: str, release: str, header: bool, file_separator: str, remove_annotations: bool) -> str:
    """Concatenate a list of autoreporting reports into one file for mySQL import
    The files are appended into one file. It is important that all of the files have the same columns.
    Args:
        folderpath (str): Input folder path
        output_fname (str): output file name
        release (str): release name. Added as one of the columns in the resulting file.
        header (bool): To include header in output file or not.
        file_separator (str). File separator for output file.
    Returns:
        (str): Comma-separated columns. can be used for the --columns arg for gcloud sql import
    """
    inputdata = glob.glob("{}/*".format(folderpath))
    length = len(inputdata)
    for idx, fname in enumerate(inputdata):
        #print progress
        print("\rProgress: {:>6.3g}%, file {:>6d} of {:>6d}".format(100*idx/length,idx,length),end="")
        
        #read in report
        temp_df = pd.read_csv(fname, sep="\t")

        #add phenotype column as well as the release data
        temp_df["phenotype"] = fname.split("/")[-1].split(".")[0]
        temp_df["rel"] = release

        #remove invalid characters ('.', '#') from column names
        rename_dict = {}
        for c in temp_df.columns:
            if ("#" in c) or ("." in c):
                rename_dict[c] = c.replace("#", "" ).replace(".","")
        temp_df = temp_df.rename(columns=rename_dict)

        if remove_annotations:
            # annotation columns have been prefixed with FG_ and GENOME_ and EXOME_ prefixes. At least in the beginning, remove FG_. 
            # GENOME_ and EXOME_ do not change between releases, at least if the script is not changed.
            columns_with_fg = [a for a in list(temp_df.columns) if bool(re.match(r'^FG_', a))]
            temp_df = temp_df.drop(columns=columns_with_fg)
        
        #write the report df to a single file. In case it is not the first file, it is appended after the existing file.
        if idx == 0:
            temp_df.to_csv(output_fname,
                           sep=file_separator,
                           mode="w",
                           na_rep="NA",
                           header=header,
                           index=False,
                           quoting=csv.QUOTE_NONNUMERIC)
        else:
            temp_df.to_csv(output_fname,
                           sep=file_separator,
                           mode="a",
                           na_rep="NA",
                           header=False,
                           index=False,
                           quoting=csv.QUOTE_NONNUMERIC)

    print()
    return ",".join([a for a in temp_df.columns])


if __name__=="__main__":
    parser=argparse.ArgumentParser("Concatenate variant reports into a single file")
    parser.add_argument("--report-folder",type=str,required=True,help="variant folder")
    parser.add_argument("--output-fname",type=str,required=True,help="output filename")
    parser.add_argument("--release",type=str,required=True,help="Release column value to give to variants")
    parser.add_argument("--header",action="store_true",default=False,help="Include header in file with this flag. By default output file has no header.")
    parser.add_argument("--file-separator",type=str,default=",",help="File separator, by default ',' for compatibility with mysql csv import.")
    parser.add_argument("--remove-annotations",action="store_true",help="Remove annotation columns, as those can change based on release. Use if importing to multiple-release table.")
    args = parser.parse_args()
    cols = concat_files(args.report_folder, 
                        args.output_fname,
                        args.release,
                        args.header,
                        args.file_separator,
                        args.remove_annotations)
    print("Results written to {}".format(args.output_fname))
    print("Column names:")
    print(cols)