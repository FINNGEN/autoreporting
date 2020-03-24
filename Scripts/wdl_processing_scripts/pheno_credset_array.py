import argparse

def main(args):
    #phenotype list
    with open(args.phenotype_list,"r") as f:
        phenolist = f.readlines()
    with open(args.credset_list,"r") as f:
        credlist = f.readlines()
    #clean up
    phenolist = [a.strip() for a in phenolist]
    credlist = [a.strip() for a in credlist]
    #get paths
    a=phenolist[0]
    b=credlist[0]
    phenopath = "/".join(a.split("/")[ : ( len(a.split("/"))-1 ) ])
    credpath = "/".join(b.split("/")[ : ( len(b.split("/"))-1 ) ])
    #get suffixes
    phenosuffix = ".".join(a.split("/")[-1].split(".")[1:])
    credsuffix = ".".join(b.split("/")[-1].split(".")[1:])
    #get phenotype names
    phenos = [a.split("/")[-1].split(".")[0] for a in phenolist ]
    creds = [a.split("/")[-1].split(".")[0] for a in credlist ]
    #if prefixes, trim them
    if args.phenotype_prefix != "":
        phenos=[a.replace(args.phenotype_prefix, "") for a in phenos ]
    if args.credset_prefix != "":
        phenos=[a.replace(args.credset_prefix, "") for a in creds ]
        

    #match the rows
    output_list=[]
    for p in phenos:
        line="{}\t{}/{}{}.{}\t".format(p,phenopath, args.phenotype_prefix,p,phenosuffix)
        if p in creds:
            line=line+"{}/{}{}.{}".format(credpath, args.credset_prefix,p,credsuffix)
        line=line+"\n"
        output_list.append(line)
    #write to file
    with open(args.out, "w") as f:
        f.writelines(output_list)

if __name__ == "__main__":
    parser=argparse.ArgumentParser("Create a phenotype, summary statistic file, credible set file-array from two file list files")
    parser.add_argument("--phenotype-list",type=str,required=True,help="phenotype list (e.g. output from 'gsutil ls gs://pheno_folder'")
    parser.add_argument("--credset-list",type=str,required=True,help="SuSiE credible set list (e.g. output from 'gsutil ls gs://credset_folder'")
    parser.add_argument("--phenotype-prefix",type=str,default="",help="If the phenotypes have a prefix that is not part of the phenotype name, e.g. version number use this flag to include it.")
    parser.add_argument("--credset-prefix",type=str,default="",help="If the credible sets have a prefix that is not part of the phenotype name, e.g. version number use this flag to include it.")
    parser.add_argument("--out", type=str, required=True, help="Output file")
    args=parser.parse_args()
    main(args)