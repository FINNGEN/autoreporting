import abc
import argparse,shlex,subprocess, glob, time
from subprocess import Popen, PIPE
from typing import List, Text, Dict,Any
import pandas as pd, numpy as np
from data_access.gwcatalog_api import try_request, ResourceNotFound, ResponseFailure
from data_access.db import LDAccess


class OnlineLD(LDAccess):
    def __init__(self,url):
        self.url=url
    
    def __parse_variant(self,variant):
        parsed_variant="chr{}".format( variant.replace(":","_"))
        return parsed_variant

    def __get_range(self, chrom,pos,ref,alt,window,ld_threshold = None):
        window = max(min(window, 5000000), 100000)#range in api.finngen.fi is [100 000, 5 000 000]
        variant="{}:{}:{}:{}".format(chrom, pos, ref, alt)
        params={"variant":variant,"panel":"sisu3","variant":variant,"window":window}
        if ld_threshold:
            params["r2_thresh"]=ld_threshold
        try:
            data=try_request("GET",url=self.url,params=params)
        except ResourceNotFound as e:
            print("LD data not found (status code {}) with url {} and params {}.".format(e.parameters["status_code"],self.url,params) )
            return pd.DataFrame(columns=["variation1", "variation2", "r2"])
        except ResponseFailure as e:
            print("Error with request.")
            print(e)
            return pd.DataFrame(columns=["variation1", "variation2", "r2"])
        #parse data
        ld_data=data.json()["ld"]
        ld_data=pd.DataFrame(ld_data)
        if ld_data.empty:
            return pd.DataFrame(columns=["variation1", "variation2", "r2"])
        ld_data=ld_data[ ["variation1", "variation2", "r2"] ]
        ld_data=ld_data.append({"variation1":variant,"variation2":variant,"r2":1.00},ignore_index=True)
        return ld_data

    def get_ranges(self, variants, window,ld_threshold = None):
        data=variants[ [ "chr", "pos", "ref", "alt" ] ]
        ld_data=pd.DataFrame()
        for idx, row in variants.iterrows():
            variant_ld = self.__get_range(row["chr"],row["pos"],row["ref"],row["alt"],window,ld_threshold)
            ld_data=pd.concat([ld_data,variant_ld],ignore_index=True,sort=False)
        ld_data["variant_1"] = ld_data["variation1"].apply(self.__parse_variant)
        ld_data["variant_2"] = ld_data["variation2"].apply(self.__parse_variant)
        ld_data["chrom_1"] = ld_data["variation1"].apply(lambda x: x.split(":")[0])
        ld_data["pos_1"] = ld_data["variation1"].apply(lambda x: x.split(":")[1])
        ld_data["chrom_2"] = ld_data["variation2"].apply(lambda x: x.split(":")[0])
        ld_data["pos_2"] = ld_data["variation2"].apply(lambda x: x.split(":")[1])
        ld_data=ld_data.astype(dtype={"pos_1":int,"pos_2":int,"r2":float})
        returncolumns=["chrom_1","pos_1","variant_1","chrom_2","pos_2","variant_2","r2"]
        ld_data=ld_data[returncolumns]
        return ld_data


class PlinkLD(LDAccess):
    def __init__(self,path,memory):
        self.path=path
        self.memory=memory

    def get_ranges(self, variants, window, ld_threshold = None):
        if not ld_threshold:
            ld_threshold = 0.0
        #assume columns are: [chrom, pos, ref, alt,#variant], and window is int
        ld_data=pd.DataFrame()
        chromosomes = variants["chr"].unique()
        r=window//1000
        for chrom in chromosomes:
            chrom_variants = variants.loc[variants["chr"] == chrom,"#variant" ].copy()
            plink_prefix = "plink{}{}".format(chrom,chrom_variants.shape[0])
            plink_name="{}_variants".format(plink_prefix)
            chrom_variants.to_csv(plink_name,sep="\t",index=False, header=False)
            plink_cmd="plink --allow-extra-chr --chr {} --bfile {} --r2 --ld-snp-list {} --ld-window-r2 {} --ld-window-kb {} --ld-window 100000 --out {} --memory {}".format(
                chrom,
                self.path,
                plink_name,
                ld_threshold,
                r,
                plink_prefix,
                self.memory
            )
            pr = subprocess.Popen(shlex.split(plink_cmd),stdout=PIPE,stderr=subprocess.STDOUT,encoding='ASCII')
            pr.wait()
            plink_log = pr.stdout.readlines()
            if pr.returncode != 0:
                print("PLINK FAILURE. Error code {}".format(pr.returncode)  )
                [print(l) for l in plink_log]
                raise ValueError("Plink r2 calculation returned code {}".format(pr.returncode))
            ld_data_=pd.read_csv("{}.ld".format(plink_prefix),sep="\s+")
            ld_data=pd.concat([ld_data,ld_data_],axis="index",sort=False,ignore_index=True)
            cleanup_cmd= "rm {}".format(plink_name)
            plink_files = glob.glob( "{}.*".format(plink_prefix) )
            subprocess.call(shlex.split(cleanup_cmd)+plink_files, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if ld_data.empty:
            return pd.DataFrame(columns=["variant_1","pos_1","chrom_1","variant_2","pos_2","chrom_2","r2"])
        ld_data=ld_data.rename(columns={"SNP_A":"variant_1","BP_A":"pos_1","CHR_A":"chrom_1","SNP_B":"variant_2","BP_B":"pos_2","CHR_B":"chrom_2","R2":"r2"})
        return ld_data.astype({"pos_1":int,"pos_2":int,"r2":float})
