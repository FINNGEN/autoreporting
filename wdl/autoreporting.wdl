task report {
    #variables
    String docker
    File summ_stat
    File summ_stat_tb=summ_stat+".tbi"
    File gnomad_exome
    File gnomad_exome_tb=gnomad_exome+".tbi"
    File gnomad_genome
    File gnomad_genome_tb=gnomad_genome+".tbi"
    String ld_panel
    File ld_panel_bed=ld_panel+".bed"
    File ld_panel_bim=ld_panel+".bim"
    File ld_panel_fam=ld_panel+".fam"
    File finngen_annotation
    File finngen_annotation_tb=finngen_annotation+".tbi"
    File functional_annotation
    File functional_annotation_tb=functional_annotation+".tbi"
    File? credible_set 
    File? summary_stat_listing
    File? endpoint_listing
    #trick wdl to write the external summary stats paths as a file
    Array[File]? ext_summary_stats = if defined(summary_stat_listing) then read_lines(summary_stat_listing  ) else []
    File? local_gwcatalog

    Float sign_treshold
    Float alt_sign_treshold
    Int grouping_locus_width
    Float ld_r2
    Int plink_memory
    Float gwascatalog_pval
    Int gwascatalog_width_kb
    Float ld_treshold
    Int cpus
    Int gwascatalog_threads
    Int docker_memory
    Map[String,String?] efo_map

    Boolean group
    Boolean overlap
    Boolean include_batch_freq
    Boolean check_for_ld

    String grouping_method
    String ignore_region
    String ignore_cmd = if ignore_region != "" then "--ignore-region" else ""
    String map_pheno=sub(basename(basename(summ_stat),".gz" ),"finngen_R4_","")#might work
    String efo_codes = efo_map[map_pheno]
    String compare_style
    String db_choice
    String summary_cmd=if defined(ext_summary_stats) then "--summary-fpath" else ""
    String efo_cmd = if efo_codes != "" then "--efo-codes" else ""
    String dollar = "$"  
    String annotation_version
    command {
        mod_ld=$( echo ${ld_panel_bed} | sed 's/.bed//g' )
        
        output_filename=$(basename ${summ_stat})

        main.py ${summ_stat} --sign-treshold ${sign_treshold} --alt-sign-treshold ${alt_sign_treshold}  \
        ${true='--group' false='' group} --grouping-method ${grouping_method} --locus-width-kb ${grouping_locus_width} \
        --ld-panel-path ${dollar}mod_ld --ld-r2 ${ld_r2} --plink-memory ${plink_memory} ${true='--overlap' false='' overlap} \
        ${ignore_cmd} ${ignore_region} --ld-api plink \
        --gnomad-genome-path ${gnomad_genome} --gnomad-exome-path ${gnomad_exome} ${true='--include-batch-freq' false='' include_batch_freq} --finngen-path ${finngen_annotation} \
        --functional-path ${ functional_annotation} ${"--credible-set-file " + credible_set} --finngen-annotation-version ${annotation_version} \
        --compare-style ${compare_style} ${true='--check-for-ld' false='' check_for_ld} --ld-treshold ${ld_treshold}  \
        --ldstore-threads ${cpus} --gwascatalog-threads ${gwascatalog_threads} \
        ${summary_cmd} ${write_lines(ext_summary_stats)} ${"--endpoint-fpath " + endpoint_listing} \
        --gwascatalog-pval ${gwascatalog_pval} --gwascatalog-width-kb ${gwascatalog_width_kb} ${efo_cmd} ${efo_codes} --db ${db_choice} ${"--local-gwascatalog " + local_gwcatalog} \
        --fetch-out $output_filename.fetch.out --annotate-out $output_filename.annotate.out --report-out $output_filename.report.out --top-report-out $output_filename.top.out --ld-report-out $output_filename.ld.out
    }
    #output
    output {
        Array[File] out=glob("*.out")
    }
    runtime {
        docker: "${docker}"
        cpu: "${cpus}"
        memory: "${docker_memory} GB"
        disks: "local-disk 150 HDD"
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
        CODE
    >>>

    output{
        File complete_efo_map = "output_map"
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

workflow autoreporting{
    File phenotypelist
    File credsetlist
    String docker
    Int memory
    Int cpus

    String gnomad_exome
    String gnomad_genome
    String ld_panel
    String finngen_annotation
    String functional_annotation
    String local_gwcatalog
    String annotation_version
    String compare_style
    String db_choice
    File efo_code_file
    String ignore_region
    String primary_grouping_method
    String secondary_grouping_method
    Boolean include_batch_freq
    Float ld_treshold
    Float sign_treshold
    Float alt_sign_treshold
    Int grouping_locus_width
    Float ld_r2
    Int plink_memory
    Float gwascatalog_pval
    Int gwascatalog_width_kb
    Int gwascatalog_threads
    Boolean group
    Boolean overlap
    Boolean check_for_ld

    call fill_efo_map{
        input:additional_prefix="finngen_R4_",docker=docker,phenotypelist=phenotypelist,efomap=efo_code_file
    }

    call credset_filter {
        input:additional_prefix="finngen_R4_",phenotypelist=phenotypelist,credsetlist=credsetlist,docker=docker
    }
    Map[String,String] efo_code_map  = read_map(fill_efo_map.complete_efo_map)
    Array[File] phenotypes = read_lines(credset_filter.filtered_pheno_list)
    Array[File] crediblesets = read_lines(credset_filter.filtered_cred_list)
    Array[File] nocred_phenotypes = read_lines(credset_filter.other_pheno_list)
    scatter (i in range(length( phenotypes)) ){
        call report {
            input: summ_stat=phenotypes[i],credible_set=crediblesets[i],
            grouping_method=primary_grouping_method,docker=docker,
            docker_memory=memory, cpus=cpus, gnomad_exome=gnomad_exome,
            gnomad_genome=gnomad_genome, ld_panel=ld_panel,
            finngen_annotation=finngen_annotation, functional_annotation=functional_annotation,
            local_gwcatalog=local_gwcatalog, annotation_version=annotation_version,
            compare_style=compare_style, db_choice=db_choice, efo_map=efo_code_map,
            include_batch_freq=include_batch_freq, ignore_region=ignore_region,
            ld_treshold=ld_treshold, sign_treshold=sign_treshold, 
            alt_sign_treshold=alt_sign_treshold, grouping_locus_width=grouping_locus_width, 
            ld_r2=ld_r2, plink_memory=plink_memory, gwascatalog_pval=gwascatalog_pval, 
            gwascatalog_width_kb=gwascatalog_width_kb, gwascatalog_threads=gwascatalog_threads, 
            group=group, overlap=overlap, check_for_ld=check_for_ld
        }
    }
    scatter (i in range(length( nocred_phenotypes)) ){
        call report as nocred_report{
            input: summ_stat=nocred_phenotypes[i], docker=docker, docker_memory=memory, cpus=cpus, 
            gnomad_exome=gnomad_exome,gnomad_genome=gnomad_genome, ld_panel=ld_panel, 
            finngen_annotation=finngen_annotation, functional_annotation=functional_annotation, 
            local_gwcatalog=local_gwcatalog, annotation_version=annotation_version, 
            compare_style=compare_style, db_choice=db_choice, efo_map=efo_code_map, 
            include_batch_freq=include_batch_freq, ignore_region=ignore_region, 
            ld_treshold=ld_treshold, sign_treshold=sign_treshold, alt_sign_treshold=alt_sign_treshold, 
            grouping_locus_width=grouping_locus_width, ld_r2=ld_r2, plink_memory=plink_memory, 
            gwascatalog_pval=gwascatalog_pval, gwascatalog_width_kb=gwascatalog_width_kb, 
            gwascatalog_threads=gwascatalog_threads, group=group, overlap=overlap, 
            check_for_ld=check_for_ld, grouping_method=secondary_grouping_method
        }
    }

    

}
