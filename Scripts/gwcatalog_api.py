#! /usr/bin/python3
import json, requests

def get_all_associations(chromosome=1,bp_lower=1000000,bp_upper=2000000,p_upper=5e-8,p_lower=1e-323):
    base_url="https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/"
    association="/associations"
    payload={"p_upper":p_upper,"p_lower":p_lower,
        "reveal":"all","bp_lower":bp_lower,"bp_upper":bp_upper}
    url="{}{}{}".format(base_url,chromosome,association)
    r=requests.get(url,params=payload)
    print(r.url)
    dump=r.json()
    dumplst=[]
    if "_links" not in dump.keys():
        return None
    dumplst.append(dump["_embedded"]["associations"])
    i=1
    while "next" in dump["_links"].keys():
        print(i)
        r=requests.get(dump["_links"]["next"]["href"],params={"reveal":"all"})
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
