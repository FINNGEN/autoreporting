import abc
import argparse,shlex,subprocess, glob, time
from subprocess import Popen, PIPE
from typing import List, Text, Dict,Any, Optional
import pandas as pd, numpy as np # type: ignore
from Scripts.data_access.gwcatalog_api import try_request, ResourceNotFound, ResponseFailure
from Scripts.data_access.db import LDAccess, LDData, Variant
from Scripts.autoreporting_utils import create_variant_column


class OnlineLD(LDAccess):
    def __init__(self,url):
        self.url=url

    def get_range(self, variant: Variant, bp_range: int, ld_threshold: Optional[float]=None) -> List[LDData]:
        """Get LD data for a range around one variant
        """
        window = 2* bp_range 
        window = max(min(window, 5000000), 100000)#range in api.finngen.fi is [100 000, 5 000 000]
        variant_str="{}:{}:{}:{}".format(variant.chrom, variant.pos, variant.ref, variant.alt)
        params={"variant":variant_str,"panel":"sisu3","variant":variant_str,"window":window}
        if ld_threshold:
            params["r2_thresh"]=ld_threshold
        try:
            data=try_request("GET",url=self.url,params=params)
        except ResourceNotFound as e:
            print("LD data not found (status code {}) with url {} and params {}.".format(e.parameters["status_code"],self.url,params) )
            return [LDData(
                variant,
                variant,
                1.0
            )]
        except ResponseFailure as e:
            print("Error with request.")
            print(e)
            raise e
        #parse data
        ld_data=data.json()["ld"]
        ld_out=[]
        for d in ld_data:
            c1 = d["variation1"].split(":")[0]
            p1 = int(d["variation1"].split(":")[1])
            ref1 = d["variation1"].split(":")[2]
            alt1 = d["variation1"].split(":")[3]

            c2 = d["variation2"].split(":")[0]
            p2 = int(d["variation2"].split(":")[1])
            ref2 = d["variation2"].split(":")[2]
            alt2 = d["variation2"].split(":")[3]

            r2 = float(d["r2"])

            ld_out.append(
                LDData(
                    Variant(c1, p1, ref1, alt1),
                    Variant(c2, p2, ref2, alt2),
                    r2)
                )
        return ld_out



class PlinkLD(LDAccess):
    def __init__(self,path,memory):
        self.path=path
        self.memory=memory

    def get_range(self, variant: Variant, bp_range: int, ld_threshold_:Optional[float]=None)->List[LDData]:
        if not ld_threshold_:
            ld_threshold = 0.0
        else:
            ld_threshold = ld_threshold_
        chromosome = variant.chrom
        snp = "chr{}_{}_{}_{}".format(
            variant.chrom.replace("23","X"),
            variant.pos,
            variant.ref,
            variant.alt
        )
        kb_range = int(bp_range/1000)
        plink_name = f"plink_{snp}_ld"
        plink_cmd = f"plink --allow-extra-chr --bfile {self.path} --chr {chromosome} --r2 gz --ld-snp {snp} --ld-window-r2 {ld_threshold} --ld-window-kb {kb_range} --ld-window 100000 --out {plink_name} --memory {self.memory}"
        pr = subprocess.Popen(shlex.split(plink_cmd),stdout=PIPE,stderr=subprocess.STDOUT,encoding='ASCII')
        pr.wait()
        plink_log = pr.stdout.readlines() if pr.stdout is not None else []
        if pr.returncode != 0:
            if any([ True for a in plink_log if "Error: No valid variants specified by --ld-snp/--ld-snps/--ld-snp-list." in a]):
                print("PLINK FAILURE. Variants not found in LD panel for chromosome {}".format(chromosome))
            print("PLINK FAILURE. Error code {}".format(pr.returncode)  )
            print(*plink_log, sep="\n")
            return [LDData(variant,variant,1.0)]
        else:
            #columns are: CHR_A, BP_A, SNP_A, CHR_B. BP_B, SNP_B, R2
            ld_df=pd.read_csv("{}.ld.gz".format(plink_name),delim_whitespace=True,compression="gzip")
        #clean up the files
        cleanup_cmd= "rm {}".format(plink_name)
        plink_files = glob.glob( "{}.*".format(plink_name) )
        subprocess.call(shlex.split(cleanup_cmd)+plink_files, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        ld_data = [
            LDData(
                Variant(str(v["CHR_A"]).replace("X","23"), int(v["BP_A"]), v["SNP_A"].split("_")[2], v["SNP_A"].split("_")[3]),
                Variant(str(v["CHR_B"]).replace("X","23"), int(v["BP_B"]), v["SNP_B"].split("_")[2], v["SNP_B"].split("_")[3]),
                float(v["R2"])
            )
            for v in ld_df.to_dict('records')
        ]
        return ld_data
