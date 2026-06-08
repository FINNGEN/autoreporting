import shlex,subprocess, glob,os, time,sys
import multiprocessing as mp
from subprocess import PIPE
from typing import Dict, List,  Optional
import pandas as pd, numpy as np # type: ignore
from data_access.gwcatalog_api import try_request, ResourceNotFound, ResponseFailure
from data_access.db import LDAccess, LDData, Variant
import pysam

MAX_RETRIES=7

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
    def __init__(self,path_template:str, assume_variant1_indexed:bool=True):
        self.token_refresh_time = time.time()
        self.path_template = path_template
        # when True, fetch only a 1bp window at the lead position instead of the whole
        # bp_range window. Valid only if the LD tabix is indexed by variant1's position
        # (all of a lead's partner rows then sit at exactly that position) — true for the
        # finngen LD panels, where it is the dominant speedup (avoids scanning the whole
        # bp_range window). Defaults on; set False for a panel not indexed by variant1 pos.
        self.assume_variant1_indexed = assume_variant1_indexed
        ## NOTE: we assume that data is in chromosomes 1..22,X, and that the files are named so
        chroms = [str(a).replace("23","X") for a in range(1,24)]
        paths={str(a):path_template.replace("{CHROM}",f"{a}") for a in chroms}
        exists = [(os.path.exists(a),a) for a in paths.values()]
        # check that files are available
        self.paths = {a:b for a,b in paths.items()}
        # if paths don't exist and they aren't gs paths, raise warning
        if (not all([a[0] for a in exists])) and (not all([a[1].startswith("gs://") for a in exists])):
            print(f"Warning: When loading TabixLD, some chromosome files were not found for template {path_template}: {[a[1] for a in exists if not a[0]]}",file=sys.stderr)
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
        # chromosome spelling as stored in the file (X is "X", never "23")
        chrom_clean = variant.chrom.replace("chr","")
        file_chrom = "X" if chrom_clean in ("23","X") else chrom_clean
        # variant1 id exactly as it appears in the file; lets us reject non-matching lines
        # with a cheap substring test before paying for str.split, and avoids a per-row replace.
        file_variant_id = f"chr{file_chrom}_{variant.pos}_{variant.ref}_{variant.alt}"
        start = max(0,variant.pos-bp_range)
        end = variant.pos+bp_range
        sequence = file_chrom
        if sequence not in self.sequences:
            raise Exception(f"Error in fetching LD for variant {variant}: File for chromosome {sequence} not in available files.")
        if not ld_threshold_:
            ld_threshold = 0.0
        else:
            ld_threshold = ld_threshold_
        # when the file is indexed by variant1 pos, all rows for this lead sit at lead.pos,
        # so a 1bp fetch suffices and avoids scanning every variant1 in the window.
        if self.assume_variant1_indexed:
            fetch_start = max(0,variant.pos-1)
            fetch_end = variant.pos
        else:
            fetch_start = start
            fetch_end = end
        i_v1 = self.hdi["variant1"]
        data = []
        tries = 0
        while True:
            try:
                data = []
                iter = self.fileobjects[sequence].fetch(sequence, fetch_start,fetch_end)
                for l in iter:
                    # cheap reject: most rows in the window belong to other lead variants
                    if file_variant_id not in l:
                        continue
                    cols = l.split("\t")
                    if cols[i_v1] == file_variant_id:
                        data.append(cols)
                break
            except:
                print(f"Error loading data from region {sequence}:{fetch_start}-{fetch_end}",file=sys.stderr)
                #assume error is in accessing over gcp, so we retry
                if tries > MAX_RETRIES:
                    raise Exception(f"Accessing LD region {sequence}:{fetch_start}-{fetch_end} from file {self.paths[sequence]} failed after {tries} tries!")

                print(f"Error when accessing data from LD tabix file, assuming error is with access over GCP. Waiting {2**tries} seconds, recycling fileobject.",file=sys.stderr)
                time.sleep(2**tries)
                self.restart_fileobject(sequence)
                tries +=1

        # parse only the matched rows (outside the retry block, so parse errors don't trigger retries)
        i_v2 = self.hdi["variant2"]
        i_r2 = self.hdi["r2"]
        lddata = []
        variant1 = variant
        for d in data:
            r2 = float(d[i_r2])
            if r2 <= ld_threshold:
                continue
            v2_str = d[i_v2].split("_")
            pos2 = int(v2_str[1])
            if pos2 < start or pos2 >= end:
                continue
            variant2 = Variant(v2_str[0].replace("chr","").replace("X","23"),pos2,v2_str[2],v2_str[3])
            lddata.append(LDData(variant1,variant2,r2))
        lddata.append(
            LDData(variant1,variant1,1.0)
        )
        return lddata

    def get_ranges(self, variants: List[Variant], bp_range: int,
                   thresholds: Optional[List[float]] = None,
                   bp_ranges: Optional[List[int]] = None,
                   workers: int = 1) -> Dict[Variant, List[LDData]]:
        """Fetch LD for many leads, optionally across worker processes.

        The LD result for a lead depends only on (lead, range, threshold), never on any
        grouping state, so all fetches are independent and safe to parallelize. Pysam
        TabixFile objects aren't picklable, so each worker re-opens its own files via the
        Pool initializer.
        """
        tasks = []
        for i, v in enumerate(variants):
            r = bp_ranges[i] if bp_ranges is not None else bp_range
            t = thresholds[i] if thresholds is not None else None
            tasks.append((v, r, t))
        total = len(tasks)
        if workers is None or workers <= 1 or total <= 1:
            return {v: self.get_range(v, r, t) for v, r, t in tasks}
        # for gs:// panels make sure a fresh token exists to hand to the workers
        if self.path_template.startswith("gs://"):
            self.refresh_token()
        gcs_token = os.environ.get("GCS_OAUTH_TOKEN")
        n = min(workers, total)
        chunksize = max(1, total // (n * 4))
        out: Dict[Variant, List[LDData]] = {}
        done = 0
        with mp.Pool(n, initializer=_ld_worker_init,
                     initargs=(self.path_template, self.assume_variant1_indexed, gcs_token)) as pool:
            for variant, ld in pool.imap_unordered(_ld_worker_fetch, tasks, chunksize=chunksize):
                out[variant] = ld
                done += 1
                if done % 500 == 0 or done == total:
                    print(f"LD prefetch: {done}/{total} variants fetched", flush=True)
        return out
        
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
                if tries > MAX_RETRIES:
                    raise Exception(f"Could not reopen tabix fileobject {self.paths[sequence]} with multiple tries.")
                print(f"Error refreshing fileobject for tabixfile {self.paths[sequence]}, waiting {2**tries} seconds.",file=sys.stderr)
                time.sleep(2**tries)
                tries += 1


    def close(self):
        #NOTE: always call close after the resource is no longer needed
        for _,b in self.fileobjects.items():
            b.close()


# per-worker TabixLD, opened once in each Pool worker (pysam file objects aren't picklable)
_WORKER_LD: Optional[TabixLD] = None

def _ld_worker_init(path_template, assume_variant1_indexed, gcs_token):
    global _WORKER_LD
    if gcs_token:
        os.environ["GCS_OAUTH_TOKEN"] = gcs_token
    _WORKER_LD = TabixLD(path_template, assume_variant1_indexed=assume_variant1_indexed)

def _ld_worker_fetch(task):
    variant, bp_range, threshold = task
    return (variant, _WORKER_LD.get_range(variant, bp_range, threshold))
