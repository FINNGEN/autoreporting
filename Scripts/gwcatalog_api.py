import json, requests
import time
import abc
from typing import List, Text, Dict,Any

def parse_output(dumplst):
    rows=[]
    for d in dumplst:
        for k in d.keys():
            rows.append(d[k])
    return rows

def get_trait_name(trait):
    base_url="https://www.ebi.ac.uk/gwas/rest/api/efoTraits/"
    trait_=trait.upper()
    r=requests.get(url=base_url+trait_)
    if r.status_code == 404:
        print("Trait {} not found in GWASCatalog".format(trait))
        return "NA"
    elif r.status_code != 200:
        print("Request for trait {} returned status code {}".format(trait,r.status_code))
        return "NA"
    else:
        return r.json()["trait"]

def try_request(url, params=None,timeout=5):
    r=requests.get(url, params=params)
    tries=0
    while r.status_code not in (200,400,404):
        time.sleep(1)
        print("Request returned code {} with attempt {}".format(r.status_code,tries+1))
        r=requests.get(url, params=params)
        if tries >= timeout:
            break
        tries+=1
    return r

class ExtDB(object):
    @abc.abstractmethod
    def get_associations(self,chromosome: str,start: int,end: int,pval: float,size: int)-> List[Dict[str,Any]]:
        """ Return associations of range chr:start-end that have pval smaller than pval. Get results in at most size sized chunks.
            Args: chromosome start end pval size
            Returns: List of Dictionaries with elements "chrom":chromosome "pos":position "ref":ref_allele "alt":alt_allele "phenotype":phenotype_code "harmonized":harmonization_code
        """
        return

    def get_trait(self, trait_code : str)-> Text:
        """ Return trait given trait code
            Args: trait_code
            Returns: Trait name
        """
        return



class SummaryApi(ExtDB):
    """ 
    Gwas Catalog summary statistic API
    """
    def get_associations(self,chromosome,start,end,pval=5e-8,size=1000):
        p_lower=1e-324
        base_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/"
        association="/associations"
        payload={"p_upper":pval,"p_lower":p_lower,
            "reveal":"all","bp_lower":start,"bp_upper":end,"start":0,"size":size}
        r=try_request(url,params=payload)
        if r.status_code not in {200,400,404}:
            print("Request returned status code {}.\n url:{}".format(r.status_code,r.url))
            return
        if r.status_code in [400,404]:
            print("No variants found in {}:{}-{}".format(chromosome,bp_lower,bp_upper))
            return
        dump=r.json()
        dumplst=[]
        if "_links" not in dump.keys():
            return
        dumplst.append(dump["_embedded"]["associations"])
        i=1
        while "next" in dump["_links"].keys():
            time.sleep(1)
            r=try_request(dump["_links"]["next"]["href"],params={"reveal":"all"})
            #print(r.url)
            if r.status_code != 200:
                print("Request {} with params {} returned status code {}".format(url,payload,r.status_code))
                break
            dump=r.json()
            if "_links" not in dump.keys():
                break
            dumplst.append(dump["_embedded"]["associations"])
            i+=1
        retval=parse_output(dumplst)
        return retval

    def get_trait(self, trait_code):
        base_url="https://www.ebi.ac.uk/gwas/rest/api/efoTraits/"
        trait_=trait_code.upper()
        r=requests.get(url=base_url+trait_)
        if r.status_code == 404:
            print("Trait {} not found in GWASCatalog".format(trait))
            return trait_code
        elif r.status_code != 200:
            print("Request for trait {} returned status code {}".format(trait,r.status_code))
            return trait_code
        else:
            return r.json()["trait"]

class GwasApi(ExtDB):
    """ 
    Gwas Catalog + ensembl api, returning values identical to the gwas catalog website.
    """

    def in_chunks(self,lst, chunk_size):
        return (lst[pos:pos + chunk_size] for pos in range(0, len(lst), chunk_size))

    def parse_ensembl(self,json_data):
        out=[]
        for key in json_data.keys():
            minor_allele=json_data[key]["minor_allele"]
            other_allele=json_data[key]["ancestral_allele"]
            rsid=key
            out.append({"rsid":rsid,"ref":minor_allele,"alt":other_allele})
        return out

    def get_associations(self,chromosome,start,end,pval,size):
        url="https://www.ebi.ac.uk/gwas/api/search/downloads?q=chromosomeName: {} AND chromosomePosition:[ {} TO {}"\
            "]&pvalfilter={}&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true".format(
            chr,bp_lower,bp_upper,p_limit)
        gwcat_response=req.get(url)
        s_io=StringIO(gwcat_response.text)
        df=pd.read_csv(s,sep="\t")
        rsids=list(df["SNPS"])
        rsids=[a for a in rsids if ' x ' not in a] 
        ensembl_url="https://rest.ensembl.org/variation/human"
        out=[]
        for rsid_chunk in in_chunks(rsids,200):
            data={"ids":rsid_chunk}
            ensembl_response=req.post(url=ensembl_url,data=data)
            out=out+ parse_ensembl(ensembl_response.json()) 
        rsid_df=pd.DataFrame(out)
        df_out=df.merge(rsid_df,how="inner",left_on="SNPS",right_on="rsid")
        cols=["SNPS","CHR_ID","CHR_POS","ref","alt","P-VALUE","MAPPED_TRAIT","MAPPED_TRAIT_URI"]
        retval=df_out.loc[:,cols].copy()
        retval.loc[:,"trait"]=retval.loc[:,"MAPPED_TRAIT_URI"].apply(lambda x: x.split("/").pop())
        retval.loc[:,"harmonized"]=20
        rename={"CHR_ID":"chrom","CHR_POS":"pos","P-VALUE":"pval"}
        retval=retval.rename(columns=rename)
        retcols=["chrom","pos","ref","alt","pval","trait","harmonized"]
        return

    def get_trait(self, trait_code):
        return



def get_all_associations(chromosome=1,bp_lower=1000000,bp_upper=2000000,p_upper=5e-8,p_lower=2.4704e-324,size=1000):
    base_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/"
    association="/associations"
    payload={"p_upper":p_upper,"p_lower":p_lower,
        "reveal":"all","bp_lower":bp_lower,"bp_upper":bp_upper,"start":0,"size":size}
    url="{}{}{}".format(base_url,chromosome,association)
    r=try_request(url,params=payload)
    print(r.url)
    if r.status_code not in {200,400,404}:
        print("Request returned status code {}.\n url:{}".format(r.status_code,r.url))
        return
    if r.status_code in [400,404]:
        print("No variants found in {}:{}-{}".format(chromosome,bp_lower,bp_upper))
        return
    dump=r.json()
    dumplst=[]
    if "_links" not in dump.keys():
        return None
    dumplst.append(dump["_embedded"]["associations"])
    i=1
    while "next" in dump["_links"].keys():
        time.sleep(5)
        #payload["start"]+=20
        print(i)
        r=try_request(dump["_links"]["next"]["href"],params={"reveal":"all"})
        #r=requests.get(url,params=payload)
        print(r.url)
        if r.status_code != 200:
            print("Request {} with params {} returned status code {}".format(url,payload,r.status_code))
            break
        dump=r.json()
        if "_links" not in dump.keys():
            break
        dumplst.append(dump["_embedded"]["associations"])
        i+=1
    return dumplst
    
