import json, requests
import time
import abc
from typing import List, Text, Dict,Any
from io import StringIO
import pandas as pd, numpy as np

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
        print("Request {} with parameters {} returned code {} with attempt {}".format(url, params, r.status_code,tries+1))
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

    @abc.abstractmethod
    def get_trait(self, trait_code : str)-> str:
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
        p_lower=1e-323
        base_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/"
        association="/associations"
        url="{}{}{}".format(base_url,chromosome,association)
        payload={"p_upper":pval,"p_lower":p_lower,
            "reveal":"all","bp_lower":start,"bp_upper":end,"start":0,"size":size}
        r=try_request(url,params=payload)
        if r.status_code not in {200,400,404}:
            print("Request returned status code {}.\n url:{}".format(r.status_code,r.url))
            return
        if r.status_code in [400,404]:
            print("No variants found in {}:{}-{}".format(chromosome,start,end))
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
        retval_=parse_output(dumplst)
        retval=[]
        for record in retval_:
            retval.append({"chrom":record["chromosome"],"pos":record["base_pair_location"],"ref":record["hm_effect_allele"],
            "alt":record["hm_other_allele"],"pval":record["p_value"],"trait":record["trait"][0],"code":record["hm_code"]})
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

def parse_efo(code):
    if type(code) != type("asd"):
        print("INVALID TYPE OF EFO CODE with code {}, type {}".format(code,type(code)))
        return "NAN"
    else: 
        return code.split("/").pop()

class GwasApi(ExtDB):
    """ 
    Gwas Catalog + ensembl api, returning values identical to the gwas catalog website.
    """

    def __in_chunks(self,lst, chunk_size):
        return (lst[pos:pos + chunk_size] for pos in range(0, len(lst), chunk_size))

    def __parse_ensembl(self,json_data):
        out=[]
        for key in json_data.keys():
            #minor_allele=json_data[key]["minor_allele"]
            #other_allele=json_data[key]["ancestral_allele"]
            alleles=json_data[key]["mappings"][0]["allele_string"].split("/")
            other_allele=alleles[0]
            minor_allele=alleles[1]
            rsid=key
            out.append({"rsid":rsid,"ref":minor_allele,"alt":other_allele})
        return out

    def get_associations(self,chromosome,start,end,pval=5e-8,size=1000):
        url="https://www.ebi.ac.uk/gwas/api/search/downloads?q=chromosomeName: {} AND chromosomePosition:[ {} TO {}"\
            "]&pvalfilter={}&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true".format(
            chromosome,start,end,pval)
        gwcat_response=requests.get(url)
        s_io=StringIO(gwcat_response.text)
        df=pd.read_csv(s_io,sep="\t")
        if df.empty:
            return
        rsids=list(df["SNPS"])
        rsids=[a for a in rsids if ' x ' not in a] 
        ensembl_url="https://rest.ensembl.org/variation/human"
        out=[]
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        for rsid_chunk in self.__in_chunks(rsids,200):
            #list_str="', '".join(rsid_chunk)
            list_str='["{}"]'.format('", "'.join(rsid_chunk))
            data='{{ "ids":{} }}'.format(list_str) #{"ids":rsid_chunk}
            ensembl_response=requests.post(url=ensembl_url,headers=headers,data=data)
            if ensembl_response.status_code != 200:
                print(ensembl_response.text)
            out=out+ self.__parse_ensembl(ensembl_response.json()) 
        rsid_df=pd.DataFrame(out)
        df_out=df.merge(rsid_df,how="inner",left_on="SNPS",right_on="rsid")
        cols=["SNPS","CHR_ID","CHR_POS","ref","alt","P-VALUE","MAPPED_TRAIT","MAPPED_TRAIT_URI"]
        tmpdf=df_out.loc[:,cols].copy()
        #deal with multiple efo codes in retval trait uri column
        retval=pd.DataFrame()
        for _,row in tmpdf.iterrows():
            if type(row["MAPPED_TRAIT_URI"]) == type("string"):
                if "," in row["MAPPED_TRAIT_URI"]:
                    efos=row["MAPPED_TRAIT_URI"].split(",")
                    efos=[e.strip() for e in efos]
                    for efo in efos:
                        new_row=row.copy()
                        new_row["MAPPED_TRAIT_URI"]=efo
                        retval=retval.append(new_row)
                else:
                    retval=retval.append(row)
            else:
                retval=retval.append(row)
        retval=retval.reset_index()
        retval.loc[:,"trait"]=retval.loc[:,"MAPPED_TRAIT_URI"].apply(lambda x: parse_efo(x))
        retval.loc[:,"code"]=20
        rename={"CHR_ID":"chrom","CHR_POS":"pos","P-VALUE":"pval"}
        retval=retval.rename(columns=rename)
        retval=retval.astype(dtype={"chrom":int})
        retval=retval.astype(dtype={"chrom":str,"pos":int,"ref":str,"alt":str,"pval":float,"trait":str,"code":int})
        retcols=["chrom","pos","ref","alt","pval","trait","code"]
        return retval.loc[:,retcols].to_dict("records")

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
    
