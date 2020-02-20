import argparse
import pandas as pd
from typing import List

def main(filelist: List[str], columns: List[str], input_sep: str) -> pd.DataFrame :
    column_dict = {
                columns[0]: "chrom",
                columns[1]: "pos",
                columns[2]: "ref",
                columns[3]: "alt",
                columns[4]: "pval",
                columns[5]: "beta",
                columns[6]: "se",
                columns[7]: "study_doi",
                columns[8]: "trait"}
    data= pd.DataFrame(columns=list( column_dict.values() )  )
    
    for f in filelist:
        datafile = pd.read_csv(f,sep=input_sep)
        #check that all of the columns exist 
        if not all(col in datafile.columns for col in list( column_dict.keys() ) ):
            raise KeyError("Datafile {} does not contain all necessary columns! missing columns: {}".format(f, [a for a in column_dict.keys() if a not in datafile.columns]  ) )
        #rename columns
        datafile = datafile.rename(columns=column_dict)
        data=pd.concat([data,datafile],sort=False, ignore_index=True,axis="index")
    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compile phenotype<->trait association files into one single datafile")
    parser.add_argument("files",nargs="+",help="input files")
    parser.add_argument("--columns",type=str, nargs=9, required=True,help="Required columns, i.e. chrom, pos, ref, alt, pval, beta, se, study_doi, trait. Other columns will be included, but won't get renamed nor checked for existence." )
    parser.add_argument("--input-sep",type=str, default="\t", help="Input file separator")
    parser.add_argument("--out",type=str,required=True,help="Output file name")
    args=parser.parse_args()
    output = main(args.files, args.columns, args.input_sep)
    output.to_csv(args.out, sep="\t",index=False,na_rep="NA")
    print("Wrote output file to {}".format(args.out) )