task credset_filter{
    #divide phenotypelist and credsetlist into
    # output_phenolist: phenotypes for which there is a credible set file
    # output_credlist: credible sets for the aforementioned phenotypes
    # output_other_phenolist: phenotypes for which there is no credible sets.
    String additional_prefix
    String docker
    File phenotypelist
    File credsetlist
    command <<<
        python3 <<CODE
        with open("${phenotypelist}","r") as f:
            with open("${credsetlist}","r") as f2:
                ph_list = f.readlines()
                cr_list = f2.readlines()
                ph_list= [x.strip("\n") for x in ph_list]
                cr_list= [x.strip("\n") for x in cr_list]
                ph_path = "/".join( ph_list[0].split("/")[:-1] )
                ph_suffix = ".".join( ph_list[0].split(".")[1:] )
                cr_path = "/".join( cr_list[0].split("/")[:-1] )
                cr_suffix = ".".join( cr_list[0].split(".")[1:] )
                ph_list = [x.split("/")[-1] for x in ph_list]
                cr_list = [x.split("/")[-1] for x in cr_list]
                ph_list = [x.split(".")[0] for x in ph_list]
                cr_list = [x.split(".")[0] for x in cr_list]
                if "${additional_prefix}" != '':
                    ph_list = [x.split("${additional_prefix}")[1] for x in ph_list ]#remove finngen_R4_ from ph_list
                common_phenos=[x for x in ph_list if x in cr_list]
                other_phenos = [x for x in ph_list if x not in cr_list]
                pheno_list = ["{}/{}{}.{}\n".format(ph_path,"${additional_prefix}",x,ph_suffix) for x in common_phenos ]
                other_pheno_list = ["{}/{}{}.{}\n".format(ph_path,"${additional_prefix}",x,ph_suffix) for x in other_phenos ]
                credible_set_list=["{}/{}.{}\n".format(cr_path,x,cr_suffix) for x in common_phenos]
                with open("output_phenolist","w") as f3:
                    f3.writelines(pheno_list)
                with open("output_credlist","w") as f4:
                    f4.writelines(credible_set_list)
                with open("output_other_phenolist","w") as f5:
                    f5.writelines(other_pheno_list) 
        CODE
    >>>

    output {
        File filtered_pheno_list = "output_phenolist"
        File filtered_cred_list = "output_credlist"
        File other_pheno_list = "output_other_phenolist"
    }
    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }

}

task simple_preprocess_chunks{
    File phenotypelist
    Int num
    String docker
    command <<<
        python3 <<CODE
        with open("${phenotypelist}","r") as f:
            lines = f.readlines()
            lines=[l.strip("\n") for l in lines]
            n=int(${num})
            f2=open("phenofiles.tbi.tsv","w")
            with open("phenofiles.tsv","w") as w:
                for i in range(len(lines)//n+1):
                    tmparr=lines[(i*n):(i+1)*n]
                    tbiarr=["".join([a.strip("\n"),".tbi"]) for a in tmparr]
                    if tmparr==[]:
                        break
                    w.write("\t".join( tmparr ))
                    f2.write("\t".join(tbiarr ) )
                    if len(tmparr)<(len(lines)//n+1):
                        w.write("".join(["\t"]*(len(lines)//n+1-len(tmparr)  ) ) )
                        f2.write("".join(["\t"]*(len(lines)//n+1-len(tmparr)  ) ) )
                    w.write("\n")
                    f2.write("\n")
            f2.close()
        CODE
    >>>
    output{
        File out_tsv="phenofiles.tsv"
        File out_tbi_tsv="phenofiles.tbi.tsv"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }
}

task fill_efo_map{
    String additional_prefix
    String docker
    File phenotypelist
    File efomap
    command <<<
        python3 <<CODE
        import pandas as pd
        pheno_fname="${phenotypelist}"
        map_fname="${efomap}"
        pheno_data=pd.read_csv(pheno_fname,header=None, names=["phenotype"])
        pheno_data["phenotype"]=pheno_data["phenotype"].apply(lambda x: x.split(".")[0].split("/")[-1].replace("${additional_prefix}","") )
        map_data=pd.read_csv(map_fname, sep="\t",header=None, names=["phenotype_efo","efos"],na_values="")
        #combine data, discard phenotype column with missing values. Since it's a left join, key order is kept.
        combined_data=pheno_data.merge(map_data,how="left",left_on="phenotype",right_on="phenotype_efo").loc[:,["phenotype","efos"]]
        combined_data.to_csv("output_map",sep="\t",header=False,index=False,na_rep="\"\"")
        combined_data["efos"].to_csv("efo_arr",sep="\t",header=False, index=False, na_rep="\"\"")
        CODE
    >>>

    output{
        File complete_efo_map = "output_map"
        File efo_arr = "efo_arr"
    }
    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }
}