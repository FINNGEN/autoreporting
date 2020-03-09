import json, requests
import time
import abc
from typing import List, Text, Dict,Any
from io import StringIO
import pandas as pd, numpy as np

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


class ExtDB(object):

    @abc.abstractmethod
    def get_associations(self,chromosome: str,start: int,end: int,pval: float,size: int)-> List[Dict[str,Any]]:
        """ Return associations of range chr:start-end that have pval smaller than pval. Get results in at most size sized chunks.
            Args: chromosome start end pval size
            Returns: List of Dictionaries with elements "chrom":chromosome "pos":position "ref":ref_allele "alt":alt_allele "pval":p-value "trait":phenotype_code
        """
        return

    @abc.abstractmethod
    def get_trait(self, trait_code : str)-> str:
        """ Return trait given trait code
            Args: trait_code
            Returns: Trait name
        """
        return

def parse_output(dumplst):
    rows=[]
    for d in dumplst:
        for k in d.keys():
            rows.append(d[k])
    return rows

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

class SummaryApi(ExtDB):
    """ 
    Gwas Catalog summary statistic API
    """
    def get_associations(self,chromosome,start,end,pval=5e-8,size=1000):
        pval=float(pval)
        start=int(start)
        end=int(end)
        size=int(size)
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
            return
        except ResponseFailure:
            print("Request failure for {}:{}-{}".format(chromosome,start,end))
            return
        dump=r.json()
        dumplst=[]
        if "_links" not in dump.keys():
            return
        dumplst.append(dump["_embedded"]["associations"])
        i=1
        while "next" in dump["_links"].keys():
            try:
                r=try_request("GET",url,params=payload)
            except ResourceNotFound:
                print("No variants found in {}:{}-{}".format(chromosome,start,end))
                return
            except ResponseFailure:
                print("Request failure for {}:{}-{}".format(chromosome,start,end))
                return
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
            if record["code"] not in [9, 14, 15, 16, 17, 18]:
                retval.append({"chrom":record["chromosome"],"pos":record["base_pair_location"],"ref":record["hm_effect_allele"],
                "alt":record["hm_other_allele"],"pval":record["p_value"],"trait":record["trait"][0]})
        return retval

    def get_trait(self, trait_code):
        return get_trait_name(trait_code)

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

class LocalDB(ExtDB):

    def __init__(self,db_path):
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
        #self.df=self.df.astype(dtype={"CHR_POS":np.int64,"P-VALUE":np.float})
        self.df=self.df.dropna(axis="index",subset=["CHR_POS","CHR_ID","P-VALUE","PVALUE_MLOG"])
    
    def get_associations(self,chromosome,start,end,pval=5e-8,size=1000):
        pval=float(pval)
        start=int(start)
        end=int(end)
        size=int(size)
        #filter based on chromosome, start, end, pval
        df=self.df.loc[self.df["CHR_ID"]==str(chromosome),:].copy()
        df=df.astype({"CHR_POS":int,"P-VALUE":float})
        df=df.loc[ (df["CHR_POS"] >=int(start))& (df["CHR_POS"] <=int(end) ) ,:]
        df=df.loc[df["P-VALUE"]<=float(pval) ,:]
        rsids=list(df["SNPS"])
        rsids=[a for a in rsids if ' x ' not in a] #filter out variant*variant-interaction associations
        rsids=[a.split(";")[0].strip() for a in rsids]#some have multiple rsid codes.
        out = get_rsid_alleles_ensembl(rsids)
        rsid_df=pd.DataFrame(out,columns=["rsid","ref","alt"])
        if rsid_df.empty:
            return
        df_out=df.merge(rsid_df,how="inner",left_on="SNPS",right_on="rsid")
        cols=["SNPS","CHR_ID","CHR_POS","ref","alt","P-VALUE","PVALUE_MLOG","MAPPED_TRAIT","MAPPED_TRAIT_URI","LINK","STUDY"]
        tmpdf=df_out.loc[:,cols].copy()
        #deal with multiple efo codes in retval trait uri column
        retval = split_traits(tmpdf)
        if retval.empty:
            return
        retval=retval.reset_index()
        retval.loc[:,"trait"]=retval.loc[:,"MAPPED_TRAIT_URI"].apply(lambda x: parse_efo(x))
        rename={"CHR_ID":"chrom","CHR_POS":"pos","P-VALUE":"pval","PVALUE_MLOG":"pval_mlog","MAPPED_TRAIT":"trait_name","STUDY":"study","LINK":"study_link"}
        retval=retval.rename(columns=rename)
        retval=retval.astype(dtype={"chrom":str,"pos":int,"ref":str,"alt":str,"pval":float,"trait":str})
        retcols=["chrom","pos","ref","alt","pval","pval_mlog","trait","trait_name","study","study_link"]
        return retval.loc[:,retcols].to_dict("records")
    
    def get_trait(self, trait_code):
        return get_trait_name(trait_code)

class GwasApi(ExtDB):
    """ 
    Gwas Catalog + ensembl api, returning values identical to the gwas catalog website.
    """

    def get_associations(self,chromosome,start,end,pval=5e-8,size=1000):
        pval=float(pval)
        start=int(start)
        end=int(end)
        size=int(size)
        url="https://www.ebi.ac.uk/gwas/api/search/downloads?q=chromosomeName: {} AND chromosomePosition:[ {} TO {}"\
            "]&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true".format(
            chromosome,start,end)
        #gwcat_response=requests.get(url)
        try:
            gwcat_response=try_request("GET",url=url)
        except ResourceNotFound:
            return
        except ResponseFailure as e:
            print(e,e.parameters)
            return
        s_io=StringIO(gwcat_response.text)
        df=pd.read_csv(s_io,sep="\t")
        if df.empty:
            return
        rsids=list(df["SNPS"])
        rsids=[a for a in rsids if ' x ' not in a] #filter out variant*variant-interaction associations
        rsids=[a.split(";")[0].strip() for a in rsids]#some have multiple rsid codes.
        out = get_rsid_alleles_ensembl(rsids)
        rsid_df=pd.DataFrame(out,columns=["rsid","ref","alt"])
        df_out=df.merge(rsid_df,how="inner",left_on="SNPS",right_on="rsid")
        cols=["SNPS","CHR_ID","CHR_POS","ref","alt","P-VALUE","PVALUE_MLOG","MAPPED_TRAIT","MAPPED_TRAIT_URI","LINK","STUDY"]
        tmpdf=df_out.loc[:,cols].copy()
        #deal with multiple efo codes in retval trait uri column
        retval = split_traits(tmpdf)
        if retval.empty:
            return
        retval=retval.reset_index()
        retval.loc[:,"trait"]=retval.loc[:,"MAPPED_TRAIT_URI"].apply(lambda x: parse_efo(x))
        rename={"CHR_ID":"chrom","CHR_POS":"pos","P-VALUE":"pval","PVALUE_MLOG":"pval_mlog","MAPPED_TRAIT":"trait_name","STUDY":"study","LINK":"study_link"}
        retval=retval.rename(columns=rename)
        retval=retval.astype(dtype={"chrom":str,"pos":int,"ref":str,"alt":str,"pval":float,"trait":str})
        retcols=["chrom","pos","ref","alt","pval","pval_mlog","trait","trait_name","study","study_link"]
        return retval.loc[:,retcols].to_dict("records")

    def get_trait(self, trait_code):
        return get_trait_name(trait_code)


