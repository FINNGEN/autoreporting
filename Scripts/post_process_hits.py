#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse
from typing import NamedTuple, Optional, List
from data_access.db import Variant, LDData, LDAccess
import data_access.linkage as linkage
import data_access.db as db

def max_r2_correlation(variant:Variant,data_variants:List[Variant],ld_data:List[LDData])->Optional[LDData]:
    """Get strongest r2 correlation between a locus and other loci in its region
    """
    #filter to variants on the top report
    ld_data = [a for a in ld_data if ((a.variant1 in data_variants) and (a.variant2 in data_variants))]
    #if for some reason the variant in question is in variant2, add flipped LDData so it's not filtered away
    ld_data.extend([LDData(a.variant2,a.variant2,a.r2) for a in ld_data])
    #filter the LD so that 1) the first variant in pairs is the variant, and 2) exclude self-correlation(which is by definition strongest)
    ld_data = [a for a in ld_data if ( (a.variant1 == variant) and (a.variant2 != variant) )]
    if ld_data:
        ld_max = max(ld_data,key=lambda x: x.r2)
        return ld_max
    return None


def is_strongest_association(var:Variant,loc_id: str, pval_threshold:float,region_width: int,data:pd.DataFrame)->bool:
    """Check if the locus is the strongest association in the region.  
    """
    chrom = str(var.chrom)
    pos_min = int(var.pos)-region_width
    pos_max = int(var.pos)+region_width
    stronger_hits_idx = ((data["pval"]<pval_threshold)&
            (data["chrom"].astype(str)==chrom)&
            (data["pos"].astype(int)>pos_min)&
            (data["pos"].astype(int)<pos_max)&
            (data["locus_id"].astype(str) != loc_id))
    if stronger_hits_idx.any():
        return False
    return True


def main(data:pd.DataFrame,region_width:int,ld_api:LDAccess)->pd.DataFrame:
    """Calculate r2 between rows, and check if the locus it the strongest hit in region
    """
    #go through the data row by row
    #for each row, get LD data
    #filter that by the variants in current chromosome
    #pick the largest value or leave empty
    #assign r2_max column and r2_max_var columns to row, or create new data, whatever
    #check for strongest_hit by checking if there are any variants within the same chromosome with distance less than region_width and pval<my_pval
    #assign that to the row
    out_data = data.copy(deep=True)
    out_data["max_r2"]=0.0
    out_data["max_r2_hit"]=""
    out_data["strongest_hit"]=False
    #create list of variants in our data
    list_of_datavars = [Variant(str(d.chrom),int(d.pos),str(d.ref),str(d.alt)) for d in data.itertuples()]

    for idx,dtuple in enumerate(data.itertuples()):
        print(idx)
        #create easier to work with Variant from the row named tuple
        variant = Variant(str(dtuple.chrom),int(dtuple.pos),dtuple.ref,dtuple.alt)


        r2=0.0
        max_var = ""

        # Part 1: calculate r2
        ld_data = ld_api.get_range(variant,region_width)

        ld_max = max_r2_correlation(variant,list_of_datavars,ld_data)

        if ld_max:
            r2 = ld_max.r2
            max_var = f"chr{ld_max.variant2.chrom}_{ld_max.variant2.pos}_{ld_max.variant2.ref}_{ld_max.variant2.alt}"

        # Part 2: calculate whether there is any variant that is a stronger association than our variant "in the region"
        strongest_hit = is_strongest_association(variant,dtuple.locus_id,dtuple.pval,region_width,data)

        # write results for this locus to data
        out_data.loc[idx,"max_r2"] = r2
        out_data.loc[idx,"max_r2_hit"] = max_var
        out_data.loc[idx,"strongest_hit"] = strongest_hit
    return out_data

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("input_data",help="input top report")
    parser.add_argument("--output",required=True)
    parser.add_argument("--region-width-kb",type=int,default=2000,help="region width in kb")
    parser.add_argument("--ld-panel-path",required=True)
    parser.add_argument("--plink-memory",type=int,default=17000)
    args=parser.parse_args()
    #load prerequisites
    ld_api = linkage.PlinkLD(args.ld_panel_path,args.plink_memory)
    region_width_bp = args.region_width_kb*1000
    data=pd.read_csv(args.input_data,sep="\t")
    #calc
    out_data=main(data, region_width_bp, ld_api)
    #write output
    out_data.fillna("NA").replace("","NA").to_csv(args.output,sep="\t",index=False,float_format="%.3g")