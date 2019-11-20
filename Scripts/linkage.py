import abc
import argparse,shlex,subprocess, glob
from subprocess import Popen, PIPE
from typing import List, Text, Dict,Any
import pandas as pd, numpy as np
from gwcatalog_api import try_request


class LDAccess(object):
    """
    Abstract object for getting LD
    """

    @abc.abstractmethod
    def get_ranges(self, variants: pd.DataFrame, window: int, ld_threshold: float) -> pd.DataFrame:#List[ Dict[str, Any ]]:
        """
        Return LD for multiple variant ranges
        In: variant data, i.e. a dataframe with columns [chr, pos, ref, alt, #variant], a window,ld threshold
        Out: Dataframe with LD information
        """
        return


def parse_variant(variant):
    parsed_variant="chr{}".format( "_".join(variant.split(":")) )
    return parsed_variant


class OnlineLD(LDAccess):
    def __init__(self):
        self.url="http://api.finngen.fi/api/ld"
    
    def __get_range(self, chrom,pos,ref,alt,window,ld_threshold):
        variant="{}:{}:{}:{}".format(chrom, pos, ref, alt)
        params={"variant":variant,"panel":"sisu3","variant":variant,"window":window,"r2_thresh":ld_threshold}
        data=try_request(self.url,params=params)
        if type(data) == type(None):
            print("LD data not found.") #TODO: handle this better. Perhaps return more information from try_request
            return None
        #parse data
        ld_data=data.json()["ld"]
        ld_data=pd.DataFrame(ld_data)
        ld_data=ld_data[ ["variation1", "variation2", "r2"] ]
        return ld_data

    def get_ranges(self, variants, window,ld_threshold):
        #The api works with requesting chr:pos:ref:alt keys
        data=variants[ [ "chr", "pos", "ref", "alt" ] ]
        ld_data=pd.DataFrame()
        for idx, row in variants.iterrows():
            variant_ld = self.__get_range(row["chr"],row["pos"],row["ref"],row["alt"],window,ld_threshold)
            ld_data=pd.concat([ld_data,variant_ld],ignore_index=True,sort=False)
        #parse the LD data. parsing single columns at a time because result_type="expand" adds horrible (8x current running time) overhead.
        ld_data["variant_1"] = ld_data["variation1"].apply(parse_variant)
        ld_data["variant_2"] = ld_data["variation2"].apply(parse_variant)
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

    def get_ranges(self, variants, window, ld_threshold):
        #assume columns are: [chrom, pos, ref, alt,#variant], and window is int
        #calculate per-chromosome to save up on memory. Performance should be similar, if not faster.
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
        ld_data=ld_data.rename(columns={"SNP_A":"variant_1","BP_A":"pos_1","CHR_A":"chrom_1","SNP_B":"variant_2","BP_B":"pos_2","CHR_B":"chrom_2","R2":"r2"})
        print(ld_data.columns)
        return ld_data.astype({"pos_1":int,"pos_2":int,"r2":float})
        #raise NotImplementedError("Not implemented yet")

# class LDStoreLD(LDAccess):
#     """
#     Might be a bit tricky, as the ld panel has to be separated into chromosomes.
#     """
#     def __get_range(self,chrom,pos,ref,alt,window):
#         raise NotImplementedError("Not implemented yet. Use get_ranges.")

#     def get_ranges(self, variants, window):
#         raise NotImplementedError("Not implemented yet")


#c=OnlineLD()
#c=PlinkLD("../../../imputation_panel/wgs_all",17000)
#data=pd.DataFrame({ "chr":["X"], "pos":[130355577], "ref": ["A"], "alt":["C"],"#variant":"chrX_130355577_A_C"  })
#print(c.get_ranges(data,1500000,0.5))