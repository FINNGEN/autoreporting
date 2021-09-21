#!/usr/bin/env python3
"""Filter top report loci
Filter out two types of loci:
1) loci that are not part of FinnGen: They have passed meta-analysis, autoreporting etc because they are very significant in UKBB, Estonia studies, but they aren't part of the FinnGen 
2) rare loci that probably are part of another, stronger locus.
"""
import pandas as pd
import numpy as np
import scipy.stats as stats 
import argparse
from typing import NamedTuple, List, Dict, Optional

class Locus(NamedTuple):
    locid: str
    c: str
    p: int
    r: str
    a: str
    pval: float
    af: float

class NotInFinnGen(NamedTuple):
    locid:str
    pval: str

class FilteredOut(NamedTuple):
    locid:str
    pval:str
    fg_af: str
    filtering_locus: str
    filtering_ld: str
    filtering_dist: str


def part_of_stronger_hit(locus:Locus,
    data:pd.DataFrame,
    r2_threshold:float,
    region_width:int,
    af_threshold: float,
    )->Optional[FilteredOut]:
    """
    Check if a locus is near a stronger signal, and if it is, return a FilteredOut object with information about the (to be removed) locus and its remover
    """
    #get values from line
    pos_min = locus.p-region_width
    pos_max = locus.p+region_width

    if locus.af < af_threshold:
        #larger hits 
        #   - are in same chromosome, within the region
        #   - have a smaller r2 than the r2 threshold
        larger_hits_in_region= ( (data["locus_id"]!=locus.locid)&
            (data["chrom"].astype(str)==locus.c) &
            (data["pos"] <= pos_max) &
            (data["pos"] >= pos_min) &
            (data["lead_r2_threshold"].astype(float) <= r2_threshold) )

        if larger_hits_in_region.any():
            possible_strong_hits = data.loc[larger_hits_in_region,:]
            most_significant_hit = possible_strong_hits.loc[possible_strong_hits["pval"].idxmin(),:]
            sig_locus_id = most_significant_hit["locus_id"]
            sig_ld = most_significant_hit["lead_r2_threshold"]
            sig_distance = abs(most_significant_hit["pos"]-locus.p)
            return FilteredOut(
                locus.locid,
                str(locus.pval),
                str(locus.af),
                str(sig_locus_id),
                str(sig_ld),
                str(sig_distance)
            )
    return None

def main(summstat:str,output:str,width:int,r2_threshold:float,af_col:str,af_threshold:float,fg_specific_column:str):
    """
    Filter out bad rows that are in presence of a row that most likely is responsible for those rows
    """
    #if pos +-2MB, r2_thresh < 0.02, self AF < 0.01 then filter out
    data=pd.read_csv(summstat,sep="\t")
    #we require the following columns: chrom, pos, ref, alt, pval, locus_id, af_col,fg_specific_column
    required_cols = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "pval",
        "locus_id",
        af_col,
        fg_specific_column
    ]
    if not all([a in data.columns for a in required_cols]):
        raise Exception(f"Not all required columns in dataframe: columns {[a for a in required_cols if a not in data.columns]} missing from top report")
    
    #fill in lead_r2_threshold
    if "lead_r2_threshold" not in data.columns:
        data["lead_r2_threshold"] = 5.0/stats.chi2.isf(data["pval"],df=1)
    
    not_in_fg=[]
    filtered_out = []
    
    with open(summstat) as infile:
        with open(output,"w") as outfile:

            headerline = infile.readline()
            header_ord = {a:i for i,a in enumerate(headerline.strip("\n").split("\t"))}

            outfile.write(headerline)

            for line in infile:

                line_cols = line.strip("\n").split("\t")

                #if finngen-specific values are not present for this line, then the locus is not in finngen data and we don't care about it -> filter out
                if line_cols[header_ord[fg_specific_column]] == "NA":
                    not_in_fg.append(NotInFinnGen(
                        line_cols[header_ord["locus_id"]],
                        line_cols[header_ord["pval"]]
                    ))
                    continue

                #build locus

                locus = Locus(
                    str(line_cols[header_ord["locus_id"]]),
                    str(line_cols[header_ord["chrom"]]),
                    int(line_cols[header_ord["pos"]]),
                    str(line_cols[header_ord["ref"]]),
                    str(line_cols[header_ord["alt"]]),
                    float(line_cols[header_ord["pval"]]),
                    float(line_cols[header_ord[af_col]])
                )

                stronger_hit = part_of_stronger_hit(locus,data,r2_threshold,width,af_threshold)
                if stronger_hit:
                    filtered_out.append(stronger_hit)
                else:
                    outfile.write(line)
                
    with open(output+".not_in_finngen","w") as nofg:
        not_in_finngen_header = "\t".join(["locus","pval"])
        nofg.write(f"{not_in_finngen_header}\n")
        for loc in not_in_fg:
            join = "\t".join([loc.locid,loc.pval])
            nofg.write(f"{join}\n")

    with open(output+".filtered","w") as weak_hits:
        weak_hits_header = "\t".join([
            "locus",
            "pval",
            "af",
            "filterer",
            "filterer_ld_threshold",
            "filterer_dist"
        ])
        weak_hits.write(f"{weak_hits_header}\n")

        for loc in filtered_out:
            join = "\t".join([
                loc.locid,
                loc.pval,
                loc.fg_af,
                loc.filtering_locus,
                loc.filtering_ld,
                loc.filtering_dist
            ])
            weak_hits.write(f"{join}\n")

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("top_report")
    parser.add_argument("--output")
    parser.add_argument("--width-kb",type=int,default=2000,help="region width in kb")
    parser.add_argument("--r2-threshold",type=float,default=0.02,help="r2 threshold for stronger hits")
    parser.add_argument("--af-column",required=True,help="column to determine whether a locus is rare")
    parser.add_argument("--af-threshold",type=float,default=0.01,help="af threshlod for weak hits")
    parser.add_argument("--fg-specific-column",required=True,help="This column is used to determine whether the variant is in finngen or not")
    args=parser.parse_args()
    width_bp = args.width_kb*1000
    main(args.top_report,args.output,width_bp,args.r2_threshold,args.af_column,args.af_threshold,args.fg_specific_column)
