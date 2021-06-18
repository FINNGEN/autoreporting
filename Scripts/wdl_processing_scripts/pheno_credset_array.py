import argparse
from typing import List, Dict, Optional
import itertools

def construct_path_dict(phenolist: List[str], phenoprefix: str = "") -> Dict[str, str]:
    """Construct path dictionary from phenotype list and phenotype prefix
    Args:
        phenolist (List[str]): list of phenotype files
        phenoprefix (str): prefix in each filename
    Returns:
        (Dict[str, str]): Dictionary with 'phenotype': 'path' entries
    """
    phenoname_list = construct_phenoname_list(phenolist, phenoprefix)
    path_list = ["/".join(a.split("/")[: (len(a.split("/"))-1)]) for a in phenolist]
    return {a: b for (a, b) in itertools.zip_longest(phenoname_list, path_list)}

def construct_suffix_dict(phenolist: List[str], phenoprefix: str = "") -> Dict[str, str]:
    """Construct suffix dictionary from phenotype list and phenotype prefix
    Args:
        phenolist (List[str]): list of phenotype files
        phenoprefix (str): prefix in each filename
    Returns:
        (Dict[str, str]): Dictionary with 'phenotype': 'path' entries
    """
    phenoname_list = construct_phenoname_list(phenolist, phenoprefix)
    suffix_list = [".".join(a.split("/")[-1].split(".")[1:]) for a in phenolist]
    return {a: b for (a, b) in itertools.zip_longest(phenoname_list, suffix_list)}


def construct_phenoname_list(phenolist: List[str], phenoprefix: str = "") -> List[str]:
    """Construct phenotype name list from phenotype paths
    Args:
        phenolist (List[str]): List of phenotype files
        phenoprefix (str): prefix in each filename
    Returns:
        (List[str]): List of phenotype names, in same order as input
    """
    phenos = [a.split("/")[-1].split(".")[0] for a in phenolist]
    if phenoprefix != "":
        phenos = [a.replace(phenoprefix, "") for a in phenos]
    return phenos

def main(args):
    #phenotype list
    with open(args.phenotype_list,"r") as f:
        phenolist = f.readlines()
    #credible set list
    with open(args.credset_list,"r") as f:
        credlist = f.readlines()
    #previous release list
    with open(args.prev_release_list,"r") as f:
        prevlist = f.readlines()
    #clean up
    phenolist = [a.strip() for a in phenolist]
    credlist = [a.strip() for a in credlist]
    prevlist= [a.strip() for a in prevlist]
    #get paths
    phenopathdict = construct_path_dict(phenolist, args.phenotype_prefix)
    phenosuffixdict = construct_suffix_dict(phenolist,args.phenotype_prefix)
    phenos = construct_phenoname_list(phenolist, args.phenotype_prefix)
    #if args.phenotype_prefix != "":

    credpathdict = construct_path_dict(credlist, args.credset_prefix)
    credsuffixdict = construct_suffix_dict(credlist, args.credset_prefix)
    creds = construct_phenoname_list(credlist, args.credset_prefix)

    prevpathdict = construct_path_dict(prevlist, args.prev_release_prefix)
    prevsuffixdict = construct_suffix_dict(prevlist, args.prev_release_prefix)
    prevs = construct_phenoname_list(prevlist, args.prev_release_prefix)
    #match the rows
    output_list=[]
    if args.only_cred:
        phenos=[a for a in phenos if a in creds]
    for p in phenos:
        line="{}\t{}/{}{}.{}\t".format(p,phenopathdict[p], args.phenotype_prefix,p,phenosuffixdict[p])
        if p in creds:
            line=line+"{}/{}{}.{}".format(credpathdict[p], args.credset_prefix,p,credsuffixdict[p])
        else:
            line = line+"{}".format(args.empty_file_path)
        line=line+"\t"
        if p in prevs:
            line=line+"{}/{}{}.{}".format(prevpathdict[p], args.prev_release_prefix,p,prevsuffixdict[p])
        else:
            line = line+"{}".format(args.empty_file_path)
        line=line+"\n"
        output_list.append(line)
    #write to file
    with open(args.out, "w") as f:
        f.writelines(output_list)

if __name__ == "__main__":
    parser=argparse.ArgumentParser("Create a phenotype, summary statistic file, credible set file, previous release file-array from two file list files")
    parser.add_argument("--phenotype-list",type=str,required=True,help="phenotype list (e.g. output from 'gsutil ls gs://pheno_folder'")
    parser.add_argument("--credset-list",type=str,required=True,help="SuSiE credible set list (e.g. output from 'gsutil ls gs://credset_folder'")
    parser.add_argument("--prev-release-list",type=str,required=True,help="List of previous release summary statistics")
    parser.add_argument("--phenotype-prefix",type=str,default="",help="If the phenotypes have a prefix that is not part of the phenotype name, e.g. version number use this flag to include it.")
    parser.add_argument("--credset-prefix",type=str,default="",help="If the credible sets have a prefix that is not part of the phenotype name, e.g. version number use this flag to include it.")
    parser.add_argument("--prev-release-prefix",type=str,default="",help="If the previous release summary statistics have a prefix that is not part of the phenotype name, e.g. version number use this flag to include it.")
    parser.add_argument("--empty-file-path",type=str,required=True,help="File path for an empty dummy file in google cloud. Used to fill missing entries.")
    parser.add_argument("--only-cred",action="store_true",help="Include only phenotypes with credible set information")
    parser.add_argument("--out", type=str, required=True, help="Output file")
    args=parser.parse_args()
    main(args)