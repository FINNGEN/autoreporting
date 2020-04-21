import pandas as pd
import os, glob, csv
import argparse


def concat_files(folderpath: str, output_fname: str, release: str, header: bool, file_separator: str) -> str:
    """Concatenate a list of files to 
    """
    pass
    inputdata = glob.glob("{}/*".format(folderpath))
    i=0
    length=len(inputdata)
    for fname in inputdata:
        print("\rProgress: {:>6.3g}%, file {:>6d} of {:>6d}".format(100*i/length,i,length),end="")
        temp_df=pd.read_csv(fname,sep="\t")
        temp_df["phenotype"]=fname.split("/")[-1].split(".")[0]
        temp_df["rel"]=release
        rename_dict={}
        for c in temp_df.columns:
            if ("#" in c) or ("." in c):
                rename_dict[c] = c.replace("#", "" ).replace(".","")
        temp_df=temp_df.rename(columns=rename_dict)
        if i == 0:
            temp_df.to_csv(output_fname,sep=file_separator,mode="w",na_rep="NA", header=header,index=False,quoting=csv.QUOTE_NONNUMERIC)
        else:
            temp_df.to_csv(output_fname,sep=file_separator,mode="a",na_rep="NA", header=False,index=False,quoting=csv.QUOTE_NONNUMERIC)
        i=i+1
    print()
    return ",".join([a for a in temp_df.columns])


if __name__=="__main__":
    parser=argparse.ArgumentParser("Concatenate variant reports into a single file")
    parser.add_argument("--report-folder",type=str,required=True,help="variant folder")
    parser.add_argument("--output-fname",type=str,required=True,help="output filename")
    parser.add_argument("--release",type=str,required=True,help="Release column value to give to variants")
    parser.add_argument("--header",action="store_true",default=False,help="Include header in file with this flag. By default output file has no header.")
    parser.add_argument("--file-separator",type=str,default=",",help="File separator, by default ',' for compatibility with mysql csv import.")
    args=parser.parse_args()
    cols=concat_files(args.report_folder, 
        args.output_fname,
        args.release,
        args.header,
        args.file_separator)
    print("Results written to {}".format(args.output_fname))
    print("Column names:")
    print(cols)