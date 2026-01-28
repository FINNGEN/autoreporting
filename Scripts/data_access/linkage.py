import shlex,subprocess, glob,os, time,sys
from subprocess import PIPE
from typing import List,  Optional
import pandas as pd, numpy as np # type: ignore
from data_access.gwcatalog_api import try_request, ResourceNotFound, ResponseFailure
from data_access.db import LDAccess, LDData, Variant
import pysam


class OnlineLD(LDAccess):
    def __init__(self,url,panel="sisu42"):
        self.url=url
        self.panel=panel

    def get_range(self, variant: Variant, bp_range: int, ld_threshold: Optional[float]=None) -> List[LDData]:
        """Get LD data for a range around one variant
        """
        window = 2* bp_range 
        window = max(min(window, 5000000), 100000)#range in api.finngen.fi is [100 000, 5 000 000]
        variant_str="{}:{}:{}:{}".format(variant.chrom, variant.pos, variant.ref, variant.alt)
        params={"variant":variant_str,"panel":self.panel,"variant":variant_str,"window":window}
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

class TabixLD(LDAccess):
    def __init__(self,path_template:str):
        self.token_refresh_time = time.time()
        self.path_template = path_template
        ## NOTE: we assume that data is in chromosomes 1..22,X, and that the files are named so
        chroms = [str(a).replace("23","X") for a in range(1,24)]
        paths={str(a):path_template.replace("{CHROM}",f"{a}") for a in chroms}
        exists = [(os.path.exists(a),a) for a in paths.values()]
        # check that files are available
        self.paths = {a:b for a,b in paths.items()}
        if not all([a[0] for a in exists]) or not all([a[1].startswith("gs://") for a in exists]):
            print(f"Warning: When loading TabixLD, some chromosome files were not found for template {path_template}: {[a[1] for a in exists if not a[0]]}")
        self.chrompos = [
            "#chrom",
            "pos"
        ]
        self.sequences = [a for a in paths.keys()]
        self.fileobjects = {a:pysam.TabixFile(b,encoding="utf-8") for a,b in self.paths.items()}
        self.header = self.fileobjects[self.sequences[0]].header[0].split("\t")
        self.hdi= {a:i for i,a in enumerate(self.header)}

        
    
    def get_range(self, variant: Variant, bp_range: int, ld_threshold_:Optional[float]=None)->List[LDData]:
        if self.path_template.startswith("gs://"):
            # refresh tokens every 20 minutes
            current_time = time.time()
            if (current_time - self.token_refresh_time) > 1200.0:
                print("Refreshing GCS_AUTH_TOKEN",file=sys.stderr)
                self.refresh_token()
                self.token_refresh_time = current_time
        variant_id=f"chr{variant.chrom.replace('chr','')}_{variant.pos}_{variant.ref}_{variant.alt}"
        start = max(0,variant.pos-bp_range)
        end = variant.pos+bp_range
        sequence = variant.chrom.replace("23","X")
        if sequence not in self.sequences:
            raise Exception(f"Error in fetching LD for variant {variant}: File for chromosome {sequence} not in available files.")
        if not ld_threshold_:
            ld_threshold = 0.0
        else:
            ld_threshold = ld_threshold_
        data = []
        tries = 0
        while True:
            try:
                data = []
                iter = self.fileobjects[sequence].fetch(sequence, start,end)
                for l in iter:
                    cols = l.split("\t")
                    if cols[self.hdi["variant1"]].replace("chrX","chr23") == variant_id:
                        data.append(cols)
                break
            except:
                print(f"Error loading data from region {sequence}:{start}-{end}")
                #assume error is in accessing over gcp, so we retry
                if tries > 5:
                    raise Exception(f"Accessing LD region {sequence}:{start}-{end} from file {self.paths[sequence]} failed after {tries} tries!")
                
                print(f"Error when accessing data from LD tabix file, assuming error is with access over GCP. Waiting {2**tries} seconds, recycling fileobject.",file=sys.stderr)
                time.sleep(2**tries)
                self.restart_fileobject(sequence)
                tries +=1

        lddata = []
        variant1 = variant
        for d in data:
            v2_str = d[self.hdi["variant2"]].split("_")
            variant2 = Variant(v2_str[0].replace("chr","").replace("X","23"),int(v2_str[1]),v2_str[2],v2_str[3])
            r2 = float(d[self.hdi["r2"]])
            if variant2.pos < end and variant2.pos >= start and r2 > ld_threshold:
                lddata.append(LDData(variant1,variant2,r2))
        lddata.append(
            LDData(variant1,variant1,1.0)
        )
        return lddata
        
    def refresh_token(self):
        tmp = subprocess.run(shlex.split("gcloud auth print-access-token"),capture_output=True,encoding="utf-8")
        os.environ["GCS_OAUTH_TOKEN"] = tmp.stdout.strip()

    def restart_fileobject(self,sequence:str):
        #close current fobj
        self.fileobjects[sequence].close()
        self.refresh_token()
        #try to reopen the fileobject, with exponential backoff.
        tries = 0
        while True:
            try:
                fobj = pysam.TabixFile(self.paths[sequence],encoding="utf-8")
                self.fileobjects[sequence] = fobj
                break
            except:
                if tries > 5:
                    raise Exception(f"Could not reopen tabix fileobject {self.paths[sequence]} with multiple tries.")
                print(f"Error refreshing fileobject for tabixfile {self.paths[sequence]}, waiting {2**tries} seconds.")
                time.sleep(2**tries)
                tries += 1


    def close(self):
        #NOTE: always call close after the resource is no longer needed
        for _,b in self.fileobjects.items():
            b.close()
        