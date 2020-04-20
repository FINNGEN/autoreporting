import json, requests
import time
import abc
from typing import List, Text, Dict,Any
from io import StringIO
import pandas as pd, numpy as np
from data_access.db import ExtDB
from multiprocessing.dummy import Pool as ThreadPool

class RequestError(Exception):
    """
    Unrecoverable error baseclass for try_request errors
    """
    pass

class ResourceNotFound(RequestError):
    """
    Value does not exist in database or it was not found for some other reason
    """
    def __init__(self,parameters=None):
        super(RequestError,self).__init__()
        self.parameters=parameters
    pass

class ResponseFailure(RequestError):
    """
    The request aborted due to some error
    """
    def __init__(self,parameters=None):
        super(RequestError,self).__init__()
        self.parameters=parameters

class SummaryApi(ExtDB):
    """ 
    Gwas Catalog summary statistic API
    """
    def __init__(self, pval_threshold, padding, threads):
        self.__get_trait=get_trait_name
        self.pval_threshold = float(pval_threshold)
        self.pad = int(padding)
        self.threads = threads
        self.result_size=100

    def __get_associations(self,chromosome,start,end):
        pval=self.pval_threshold
        start=max(0,int(start)-self.pad)
        end=int(end)+self.pad
        size=self.result_size
        p_lower=1e-323
        base_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/"
        association="/associations"
        url="{}{}{}".format(base_url,chromosome,association)
        payload={"p_upper":pval,"p_lower":p_lower,
            "reveal":"all","bp_lower":start,"bp_upper":end,"start":0,"size":size}
        try:
            r=try_request("GET",url,params=payload)
        except ResourceNotFound:
            print("No variants found in {}:{}-{}".format(chromosome,start,end))
            return []
        except ResponseFailure:
            print("Request failure for {}:{}-{}".format(chromosome,start,end))
            return []
        dump=r.json()
        dumplst=[]
        if "_links" not in dump.keys():
            return []
        dumplst.append(dump["_embedded"]["associations"])
        i=1
        while "next" in dump["_links"].keys():
            try:
                r=try_request("GET",url,params=payload)
            except ResourceNotFound:
                print("No variants found in {}:{}-{}".format(chromosome,start,end))
                return []
            except ResponseFailure:
                print("Request failure for {}:{}-{}".format(chromosome,start,end))
                return []
            if type(r) == type(None):
                break
            dump=r.json()
            dumplst.append(dump["_embedded"]["associations"])
            if "_links" not in dump.keys():
                break
            i+=1
        retval_=parse_output(dumplst)
        retval=[]
        for record in retval_:
            if record["hm_code"] not in [9, 14, 15, 16, 17, 18]:
                retval.append({"chrom":record["chromosome"],"pos":record["base_pair_location"],"ref":record["hm_effect_allele"],
                "alt":record["hm_other_allele"],"pval":record["p_value"],"trait":record["trait"][0]})
        for record in retval:
            record["trait_name"] = self.__get_trait(record["trait"])
        return retval

    def associations_for_regions(self, regions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Return associations for a list of regions of type {"chrom": str, "min": int, "max": int }
        Args:
            regions (List[Dict[str, Any]]): The list of regions for which associations are queried
        """
        data_lst=[]
        for region in regions:
            data_lst.append([region["chrom"],region["min"],region["max"]])
        results = None
        with ThreadPool(self.threads) as pool:
            results=pool.starmap(self.__get_associations,data_lst)
        results=[r for r in results if r != None]
        results=[i for sublist in results for i in sublist]
        return results

class LocalDB(ExtDB):

    def __init__(self, db_path: str, pval_threshold: float, padding: int):
        self.__get_trait=get_trait_name
        self.pad = int(padding)
        self.pval_threshold = float(pval_threshold)
        try:
            self.df=pd.read_csv(db_path,sep="\t",low_memory=False)
        except FileNotFoundError as err:
            raise FileNotFoundError("argument {} to flag '--local-gwascatalog' not found: Does the file exist?".format(db_path))
        self.df=self.df.astype(str)
        self.df=self.df.loc[ ~ self.df["CHR_POS"].str.contains(";") ,:]
        self.df=self.df.loc[ ~ self.df["CHR_POS"].str.contains("x") ,:]
        self.df["CHR_POS"]=pd.to_numeric(self.df["CHR_POS"],errors="coerce")
        self.df["P-VALUE"]=pd.to_numeric(self.df["P-VALUE"],errors="coerce")
        self.df["PVALUE_MLOG"]=pd.to_numeric(self.df["PVALUE_MLOG"],errors="coerce")
        self.df=self.df.dropna(axis="index",subset=["CHR_POS","CHR_ID","P-VALUE","PVALUE_MLOG"])
        self.df=self.df.astype({"CHR_POS":int,"P-VALUE":float})
        self.df=self.df.loc[self.df["P-VALUE"]<=self.pval_threshold ,:] #filter the df now by pval
    
    def __get_associations(self, chromosome: str, start: int, end: int):
        start=max(0,int(start)-self.pad)
        end=int(end)+self.pad
        #filter based on chromosome, start, end, pval
        df=self.df.loc[self.df["CHR_ID"]==str(chromosome),:].copy()
        df=df.loc[ (df["CHR_POS"] >=start)& (df["CHR_POS"] <=end ) ,:]
        if df.empty:
            return []
        cols=["SNPS","CHR_ID","CHR_POS","P-VALUE","PVALUE_MLOG","MAPPED_TRAIT","MAPPED_TRAIT_URI","LINK","STUDY"]
        tmpdf=df.loc[:,cols].copy()
        #deal with multiple efo codes in retval trait uri column
        retval = split_traits(tmpdf)
        if retval.empty:
            return []
        retval=retval.reset_index()
        retval.loc[:,"trait"]=retval.loc[:,"MAPPED_TRAIT_URI"].apply(lambda x: parse_efo(x))
        rename={"CHR_ID":"chrom","CHR_POS":"pos","P-VALUE":"pval","PVALUE_MLOG":"pval_mlog","MAPPED_TRAIT":"trait_name","STUDY":"study","LINK":"study_link"}
        retval=retval.rename(columns=rename)
        retval=retval.astype(dtype={"chrom":str,"pos":int,"pval":float,"trait":str})
        retcols=["chrom","pos","SNPS","pval","pval_mlog","trait","trait_name","study","study_link"]
        return retval.loc[:,retcols].to_dict("records")

    def associations_for_regions(self, regions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Get associations for a list of regions
        Args:
            regions(List[Dict[str, Any]]): A list of dictionaries with keys 'chrom', 'min', 'max'
        Returns:
        List[Dict[str, Any]]: A List of associations for variants in the regions
        """
        out= []
        for region in regions:
            try:
                out.extend(self.__get_associations(region["chrom"],region["min"],region["max"]))
            except:
                pass
        #TODO: get rid of this ugly format switching
        result_df = pd.DataFrame(out)
        #TODO: add ensembl data
        #filter rsids
        rsids = [a["SNPS"] for a in out if '  x  ' not in a["SNPS"]]
        rsids=[a.split(";")[0].strip() for a in rsids]
        rsid_out = get_rsid_alleles_ensembl(rsids)
        rsid_df =  pd.DataFrame(rsid_out,columns=["rsid","ref","alt"])
        df_out = result_df.merge(rsid_df,how="inner",left_on="SNPS",right_on="rsid").drop(column=["rsid"])
        return df_out.to_dict("records")

class GwasApi(ExtDB):
    """ 
    Gwas Catalog + ensembl api, returning values identical to the gwas catalog website.
    """

    def __init__(self, pval_threshold, padding, threads):
        self.__get_trait=get_trait_name
        self.pval_threshold = float(pval_threshold)
        self.pad = int(padding)
        self.threads = threads

    def __get_associations(self,chromosome,start,end):
        pval=self.pval_threshold
        start=max(0,int(start)-self.pad)
        end=int(end)+self.pad
        url="https://www.ebi.ac.uk/gwas/api/search/downloads?q=chromosomeName: {} AND chromosomePosition:[ {} TO {}"\
            "]&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true".format(
            chromosome,start,end)
        #gwcat_response=requests.get(url)
        try:
            gwcat_response=try_request("GET",url=url)
        except ResourceNotFound:
            return []
        except ResponseFailure as e:
            print(e,e.parameters)
            return []
        s_io=StringIO(gwcat_response.text)
        df=pd.read_csv(s_io,sep="\t")
        if df.empty:
            return []
        cols=["SNPS","CHR_ID","CHR_POS","P-VALUE","PVALUE_MLOG","MAPPED_TRAIT","MAPPED_TRAIT_URI","LINK","STUDY"]
        tmpdf=df.loc[:,cols].copy()
        #deal with multiple efo codes in retval trait uri column
        retval = split_traits(tmpdf)
        if retval.empty:
            return []
        retval=retval.reset_index()
        retval.loc[:,"trait"]=retval.loc[:,"MAPPED_TRAIT_URI"].apply(lambda x: parse_efo(x))
        rename={"CHR_ID":"chrom","CHR_POS":"pos","P-VALUE":"pval","PVALUE_MLOG":"pval_mlog","MAPPED_TRAIT":"trait_name","STUDY":"study","LINK":"study_link"}
        retval=retval.rename(columns=rename)
        retval=retval.astype(dtype={"chrom":str,"pos":int,"pval":float,"trait":str})
        retcols=["chrom","pos","SNPS","pval","pval_mlog","trait","trait_name","study","study_link"]
        return retval.loc[:,retcols].to_dict("records")

    def associations_for_regions(self, regions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        data_lst=[]
        for region in regions:
            data_lst.append([region["chrom"],region["min"],region["max"]])
        results = None
        with ThreadPool(self.threads) as pool:
            results=pool.starmap(self.__get_associations,data_lst)
        results=[r for r in results if r != None]
        results=[i for sublist in results for i in sublist]
        #TODO: get rid of this ugly format switching
        result_df = pd.DataFrame(results)
        #TODO: add ensembl data
        #filter rsids
        rsids = [a["SNPS"] for a in results if '  x  ' not in a["SNPS"]]
        rsids=[a.split(";")[0].strip() for a in rsids]
        rsid_out = get_rsid_alleles_ensembl(rsids)
        rsid_df =  pd.DataFrame(rsid_out,columns=["rsid","ref","alt"])
        df_out = result_df.merge(rsid_df,how="inner",left_on="SNPS",right_on="rsid").drop(column=["rsid"])
        return df_out.to_dict("records")

def parse_output(dumplst):
    rows=[]
    for d in dumplst:
        for k in d.keys():
            rows.append(d[k])
    return rows

def parse_float(number):
    """
    Parse number so that it only consists of a whole number part and an exponent, i.e. 5.1e-8 -> 51e-9
    In: any floating point number (though for the purposes of this project positive values do not make sense. Should still work.)
    Out: a string representation of that float with no dots.
    Range of operation: 0 < value < 1
    Values smaller than 0 will be '0', and values larger than 1 will be '1'
    """
    if number <= 0.0:
        return '0'
    elif number >= 1.0:
        return '1'
    numpart="".join("{:e}".format(number).split("e")[0].strip("0").split("."))
    exponent=int("{:e}".format(number).split("e")[1])-len(numpart)+1
    retval="{}e{}".format(numpart,exponent)
    return retval

def get_trait_name(trait):
    base_url="https://www.ebi.ac.uk/gwas/rest/api/efoTraits/"
    try:
        r=try_request("GET",url="{}{}".format(base_url,trait) )
    except ResourceNotFound:
        print("Trait {} not found in GWASCatalog".format(trait))
        return "NA"
    except ResponseFailure:
        print("Request for trait {} failed".format(trait))
        return "NA"
    except:#unhandled exceptions
        raise
    else:
        return r.json()["trait"]

def try_request(method,url,headers="",data="", params={},retry_count=5):
    """
    Make a GET/POST request and handle possible errors
    In: method, url, headers,data,parameters, number fo times to retry the request
    Out: Request or None
    """
    try:
        r=requests.request(method=method, url=url, headers=headers,params=params,data=data)
        tries=2
        while r.status_code not in (200,400,404):
            time.sleep(0.5)
            r=requests.request(method=method, url=url, headers=headers,params=params,data=data)
            if tries >= retry_count:
                print("{} Request {} with parameters {}; headers {}; data {}; returned code {} with attempt {}".format(method,
                    url,
                    params,
                    headers,
                    data,
                    r.status_code
                    ,tries+1))
                break
            tries+=1
    except Exception as e:
        print("Request caused an exception:{}".format(e))
        raise ResponseFailure(parameters=e)
    if r.status_code in (400,404):
        raise ResourceNotFound(parameters={"url":r.url,"params":params,"headers":headers,"data":data,"status_code":r.status_code})
    if r.status_code not in (200,):
        raise ResponseFailure(parameters={"url":r.url,"params":params,"headers":headers,"data":data,"status_code":r.status_code})
    return r

def parse_efo(code):
    if type(code) != type("string"):
        print("Invalid EFO code:{}".format(code) )
        return "NAN"
    else: 
        return code.split("/").pop()

def in_chunks(lst, chunk_size):
    return (lst[pos:pos + chunk_size] for pos in range(0, len(lst), chunk_size))

def parse_ensembl(json_data):
    out=[]
    for key in json_data.keys():
        alleles=json_data[key]["mappings"][0]["allele_string"].split("/")
        other_allele=alleles[0]
        minor_allele=alleles[1]
        rsid=key
        out.append({"rsid":rsid,"ref":minor_allele,"alt":other_allele})
    return out

def get_rsid_alleles_ensembl(rsids):
    """
    For a list of rsids, return list of variant information (alleles). Uses the ensembl human genomic variation API.
    In: list of RSIDs
    Out: list of dicts, with fields ['rsid','ref','alt']
    """
    ensembl_url="https://rest.ensembl.org/variation/human"
    out=[]
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    for rsid_chunk in in_chunks(rsids,200):
        list_str='["{}"]'.format('", "'.join(rsid_chunk))
        data='{{ "ids":{} }}'.format(list_str)
        try:
            r = try_request("POST",url=ensembl_url,headers=headers,data=data)
        except ResourceNotFound as e:
            print("{} Request to {} with params {} returned status code {} ".format( "POST",ensembl_url, e.parameters, e.parameters["status_code"] ))
        except ResponseFailure:
            print("No valid response from Ensembl API. Matches regarding the following RSIDS may not be available:{}".format(", ".join(rsid_chunk)))
        else:
            try:
                out=out + parse_ensembl(r.json())
            except ValueError:
                pass
    return out

def split_traits(df,traitname="MAPPED_TRAIT",traituriname="MAPPED_TRAIT_URI"):
    """
    Split gwascatalog trait column to single traits. 
    Does this by duplicating rows with multiple traits to multiple rows, each containing only one trait.
    In: Dataframe containing information and the trait columns, 'MAPPED_TRAIT' and 'MAPPED_TRAIT_URI'
    Out: Same dataframe with the information and columns 'MAPPED_TRAIT' and 'MAPPED_TRAIT_URI'
    """
    retval=[]
    if (traitname not in df.columns) or (traituriname not in df.columns) :
        return None
    for row in df.itertuples(index=False):
        if type(getattr(row,traituriname)) == type("string"):
            if "," in getattr(row,traituriname):
                efos=getattr(row,traituriname).split(",")
                traits=getattr(row,traitname).split(",")
                efos=[e.strip() for e in efos]
                traits=[t.strip() for t in traits]
                for idx, (efo,trait) in enumerate(zip(efos,traits)):
                    new_row=row._replace(**{traitname:trait, traituriname:efo})
                    retval.append(new_row)
            else:
                retval.append(row)
        else:
            retval.append(row)
    retval=pd.DataFrame.from_records(retval,columns=df.columns)
    return retval
