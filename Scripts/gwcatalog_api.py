#! /usr/bin/python3
import json, requests

def get_all_associations(chromosome=1,bp_lower=1000000,bp_upper=2000000,p_upper=5e-8,p_lower=2.4704e-324):
    base_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/"
    association="/associations"
    payload={"p_upper":p_upper,"p_lower":p_lower,
        "reveal":"all","bp_lower":bp_lower,"bp_upper":bp_upper,"start":0}
    url="{}{}{}".format(base_url,chromosome,association)
    r=requests.get(url,params=payload)
    print(r.url)
    if r.status_code != 200:
        raise Exception("Request {} with params {} returned status code {}".format(url,payload,r.status_code))
    dump=r.json()
    dumplst=[]
    if "_links" not in dump.keys():
        return None
    dumplst.append(dump["_embedded"]["associations"])
    i=1
    while "next" in dump["_links"].keys():
        payload["start"]+=20
        print(i)
        #r=requests.get(dump["_links"]["next"]["href"],params={"reveal":"all"})
        r=requests.get(url,params=payload)
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
        raise Exception("Trait {} not found".format(trait))
    elif r.status_code != 200:
        raise Exception("Request for trait {} returned status code {}".format(trait,r.status_code))
    else:
        return r.json()["trait"]