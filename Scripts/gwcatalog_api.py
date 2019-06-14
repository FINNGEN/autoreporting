#! /usr/bin/python3
import json, requests
import time

def try_request(url, params=None,timeout=5):
    r=requests.get(url, params=params)
    tries=1
    while r.status_code not in (200,400,404):
        time.sleep(5)
        print("Request returned code {} with attempt {}".format(r.status_code,tries))
        r=requests.get(url, params=params)
        if tries > timeout:
            break
        tries+=1
    return r

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