import argparse
import pandas as pd, numpy as np

def main(args):
    df = pd.read_csv(args.file,sep=args.sep,header=None,names=["fg","efo"])
    df=df.replace("^\s$",np.nan,regex=True)
    df.to_csv(args.out, sep="\t",index=False, na_rep="NA",header=False)


if __name__ == "__main__":
    parser=argparse.ArgumentParser("DEPRECATED Create a NA-filled fg-to-efo mapping file from fg-efo map file.")
    parser.add_argument("file",type=str,help="Mapping tsv containing the FG phenotype on first column and EFO codes in second. No header.")
    parser.add_argument("--fg-column",type=str,required=True,help="FinnGen phenotype name column" )
    parser.add_argument("--efo-column",type=str,required=True,help="EFO code column" )
    parser.add_argument("--out",type=str,required=True,help="Output filename" )
    parser.add_argument("--sep",type=str, default="\t",help="Input file separator")
    args = parser.parse_args()
    main(args)