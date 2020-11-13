import abc
import argparse,shlex,subprocess, glob, time
from subprocess import Popen, PIPE
from typing import List, Text, Dict,Any, Optional
import pandas as pd, numpy as np
from data_access.gwcatalog_api import try_request, ResourceNotFound, ResponseFailure
from data_access.db import LDAccess, LDData, Variant


class OnlineLD(LDAccess):
    def __init__(self,url):
        self.url=url
    
    def __parse_variant(self,variant):
        parsed_variant="chr{}".format( variant.replace(":","_"))
        return parsed_variant

    def __get_range(self, chrom, pos, ref, alt, window, ld_threshold=None) -> List[LDData]:
        """Get LD data for a range around one variant
        """
        window = max(min(window, 5000000), 100000)#range in api.finngen.fi is [100 000, 5 000 000]
        variant="{}:{}:{}:{}".format(chrom, pos, ref, alt)
        params={"variant":variant,"panel":"sisu3","variant":variant,"window":window}
        if ld_threshold:
            params["r2_thresh"]=ld_threshold
        try:
            data=try_request("GET",url=self.url,params=params)
        except ResourceNotFound as e:
            print("LD data not found (status code {}) with url {} and params {}.".format(e.parameters["status_code"],self.url,params) )
            return None
        except ResponseFailure as e:
            print("Error with request.")
            print(e)
            return None
        #parse data
        ld_data=data.json()["ld"]
        ld_out=[]
        for d in ld_data:
            v1 = "chr"+d["variation1"].replace(":","_")
            c1 = d["variation1"].split(":")[0]
            p1 = int(d["variation1"].split(":")[1])
            ref1 = d["variation1"].split(":")[2]
            alt1 = d["variation1"].split(":")[3]

            v2 = "chr"+d["variation2"].replace(":","_")
            c2 = d["variation2"].split(":")[0]
            p2 = int(d["variation2"].split(":")[1])
            ref2 = d["variation2"].split(":")[2]
            alt2 = d["variation2"].split(":")[3]

            r2 = float(d["r2"])

            ld_out.append(
                LDData(
                    Variant(v1, c1, p1, ref1, alt1),
                    Variant(v2, c2, p2, ref2, alt2),
                    r2)
                )
        return ld_out

    def get_ranges(self, variants: List[Variant], window: int, ld_threshold: Optional[float]=None) -> List[LDData]:
        ld_data=[]
        for v in variants:
            variant_ld = self.__get_range(v.chrom,v.pos,v.ref,v.alt,window*2,ld_threshold)
            ld_data=ld_data + variant_ld
        return ld_data


class PlinkLD(LDAccess):
    def __init__(self,path,memory):
        self.path=path
        self.memory=memory

    def get_ranges(self, variants: List[Variant], window: int, ld_threshold: Optional[float]=None) -> List[LDData]:
        if not ld_threshold:
            ld_threshold = 0.0
        variants = [Variant(
            a.variant.replace("chr23","chrX"),
            a.chrom.replace("23","X"),
            a.pos,
            a.ref,
            a.alt) for a in variants] 
        #remove duplicates of variants, because PLINK errors with that
        nodups = sorted(set(variants))
        #get chromosomes, order variants under them
        chroms = sorted(set([a.chrom for a in variants]))
        chromdict={}
        for c in chroms:
            chromdict[c] = [var for var in nodups if var.chrom == c]
        ld_data=pd.DataFrame()
        r_kb=window//1000
        for chromosome, variants in chromdict.items():
            var_df = pd.DataFrame(variants)
            var_df=var_df["variant"]
            plink_prefix = "plink{}{}".format(chromosome,len(variants))
            plink_name="{}_variants".format(plink_prefix)
            var_df.to_csv(plink_name,sep="\t",index=False, header=False)
            plink_cmd="plink --allow-extra-chr --chr {} --bfile {} --r2 --ld-snp-list {} --ld-window-r2 {} --ld-window-kb {} --ld-window 100000 --out {} --memory {}".format(
                chromosome,
                self.path,
                plink_name,
                ld_threshold,
                r_kb,
                plink_prefix,
                self.memory
            )
            pr = subprocess.Popen(shlex.split(plink_cmd),stdout=PIPE,stderr=subprocess.STDOUT,encoding='ASCII')
            pr.wait()
            plink_log = pr.stdout.readlines()
            if pr.returncode != 0:
                if any([ True for a in plink_log if "Error: No valid variants specified by --ld-snp/--ld-snps/--ld-snp-list." in a]):
                    print("PLINK FAILURE. Variants not found in LD panel for chromosome {}".format(chromosome))
                print("PLINK FAILURE. Error code {}".format(pr.returncode)  )
                print(*plink_log, sep="\n")
            else:
                ld_data_=pd.read_csv("{}.ld".format(plink_prefix),sep="\s+")
                ld_data=pd.concat([ld_data,ld_data_],axis="index",sort=False,ignore_index=True)
            cleanup_cmd= "rm {}".format(plink_name)
            plink_files = glob.glob( "{}.*".format(plink_prefix) )
            subprocess.call(shlex.split(cleanup_cmd)+plink_files, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if ld_data.empty:
            return []
        ld_data = ld_data.astype({"BP_A":int,"BP_B":int,"R2":float})
        #convert to List[LDData]
        ld_data=ld_data.to_dict('records')
        ld_data = [LDData(
            Variant(d['SNP_A'].replace("chrX","chr23"),
                    str(d['CHR_A']).replace("X","23"),
                    int(d['BP_A']),
                    d['SNP_A'].split("_")[2],
                    d['SNP_A'].split("_")[3]),
            Variant(d['SNP_B'].replace("chrX","chr23"),
                    str(d['CHR_B']).replace("X","23"),
                    int(d['BP_B']),
                    d['SNP_B'].split("_")[2],
                    d['SNP_B'].split("_")[3]),
            d['R2'])
                   for d in ld_data]
        return ld_data
