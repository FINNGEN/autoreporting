task preprocess_serial{
    File input_array
    Int phenos_per_worker
    String docker
    command{
        process_serial.py ${input_array} --n ${phenos_per_worker}
    }
    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 5 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }
    output {
        File pheno_array = "pheno_array"
        File summ_array = "summ_array"
        File summ_tb_array = "summ_tb_array"
        File credset_array = "credset_array"
        File credset_map = "boolean_map"
    }
}


task report{
    #variables
    String docker
    Int memory
    Array[String] pheno_ids
    Array[File] summ_stat
    Array[File] summ_stat_tb
    Map[String, Boolean] credset_status
    Array[File?] credsets
    Int len = length(summ_stat)
    Array[File] selected_credsets=select_all(credsets)

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
    #functional annotations
    File functional_annotation
    File functional_annotation_tb=functional_annotation+".tbi"
    #credible set annotation, no tabix
    File efo_map

    #File? summary_stat_listing
    #File? endpoint_listing
    #trick wdl to write the external summary stats paths as a file
    #NOTE: does not work with partially serialized widdle. 
    #Array[File]? ext_summary_stats = if defined(summary_stat_listing) then read_lines(summary_stat_listing  ) else []
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
    Float strict_group_r2

    Boolean group
    Boolean overlap
    Boolean include_batch_freq
    Boolean check_for_ld

    String primary_grouping_method
    String secondary_grouping_method
    String ignore_region
    String compare_style
    String db_choice
    String annotation_version
    #String summary_cmd=if defined(ext_summary_stats) then "--summary-fpath" else ""
    String additional_prefix
    #String efo_cmd = if efo_codes != "" then "--efo-codes" else ""
    String dollar = "$"  

    command <<<
        python3 <<CODE
        import subprocess, shlex, json
        from subprocess import PIPE
        #do everything a wrapper should do
        #unchanged variables
        num_phenos = ${len}
        plink_path = "${ld_panel_bed}".replace(".bed","")
        gnomad_exome="${gnomad_exome}"
        gnomad_genome="${gnomad_genome}"
        finngen_annotation="${finngen_annotation}"
        functional_annotation="${functional_annotation}"
        local_gwascatalog="${"--local-gwascatalog "+ local_gwcatalog}"
        prefix="${additional_prefix}"

        sign_treshold=${sign_treshold}
        alt_sign_treshold=${alt_sign_treshold}
        grouping_locus_width=${grouping_locus_width}
        ld_r2=${ld_r2}
        plink_memory=${plink_memory}
        gwascatalog_pval=${gwascatalog_pval}
        gwascatalog_width_kb=${gwascatalog_width_kb}
        ld_treshold=${ld_treshold}
        cpus=${cpus}
        gwascatalog_threads=${gwascatalog_threads}
        strict_group_r2=${strict_group_r2}

        group="--group" if "${group}"=="true" else ""
        overlap="--overlap" if "${overlap}"=="true" else ""
        include_batch_freq="--include-batch-freq" if "${include_batch_freq}"=="true" else ""
        check_for_ld="--check-for-ld" if "${check_for_ld}"=="true" else ""
        ignore_cmd = "--ignore-region ${ignore_region}" if "${ignore_region}" != "" else ""
        compare_style = "${compare_style}"
        db_choice = "${db_choice}"
        annotation_version = "${annotation_version}"

        #process the summstats and credsets etc
        summstats="${sep=";" summ_stat}".split(";")
        phenotypes="${sep=";" pheno_ids}".split(";")

        credset_calls = []
        boolean_map = json.load(open("${write_json(credset_status)}","r"))
        credset_list="${sep=";" selected_credsets}".split(";")
        credset_count=0
        group_method=[]
        for i,pheno in enumerate(phenotypes):
            if boolean_map[pheno] != True:
                credset_calls.append("")
                group_method.append("${secondary_grouping_method}")
            else:
                credset_calls.append(credst_list[credset_count])
                group_method.append("${primary_grouping_method}")
                credset_count += 1
        #efo codes
        with open("${efo_map}","r") as f:
            efos = {a.strip().split("\t")[0] : a.strip().split("\t")[1]  for a in f.readlines()}
        efo_array=["--efo-codes {}".format(efos[a.split("/")[-1].split(".")[0].replace(prefix,"")]) if a.split("/")[-1].split(".")[0].replace(prefix,"") in efos.keys() else "" for a in summstats ]
        

        for i in range(num_phenos):
            phenotype_name=summstats[i].split("/")[-1].split(".")[0].replace(prefix,"")
            call_command=("main.py {} "
                        " --sign-treshold {} " 
                        "--alt-sign-treshold {} "
                        "{} "
                        "--grouping-method {} "
                        "--locus-width-kb {} "
                        "--ld-panel-path {} "
                        "--ld-r2 {} "
                        "--plink-memory {} "
                        "{} "
                        "--finngen-path {} "
                        "--functional-path {} "
                        "--gnomad-genome-path {} "
                        "--gnomad-exome-path {} "
                        "{} "
                        "--finngen-annotation-version {} "
                        "--compare-style {} "
                        "{} "
                        "--ld-treshold {} "
                        "--ldstore-threads {} "
                        "--gwascatalog-threads {} "
                        "--strict-group-r2 {} "
                        #"{} "
                        #"{} "
                        "--gwascatalog-pval {} "
                        "--gwascatalog-width-kb {} "
                        "--db {} "
                        "{} "
                        "{} "
                        "{} "
                        "--fetch-out {}.fetch.out "
                        "--annotate-out {}.annotate.out "
                        "--report-out {}.report.out "
                        "--top-report-out {}.top.out "
                        "--ld-report-out {}.ld.out "
                        ).format(summstats[i],
                            sign_treshold,
                            alt_sign_treshold,
                            group,
                            grouping_method,
                            grouping_locus_width,
                            plink_path,
                            ld_r2,
                            plink_memory,
                            include_batch_freq,
                            finngen_annotation,
                            functional_annotation,
                            gnomad_genome,
                            gnomad_exome,
                            credsets[i],
                            annotation_version,
                            compare_style,
                            check_for_ld,
                            ld_treshold,
                            cpus,
                            gwascatalog_threads,
                            strict_group_r2,
                            #summary_cmd,
                            #endpoint_listing,
                            gwascatalog_pval,
                            gwascatalog_width_kb,
                            db_choice,
                            local_gwascatalog,
                            efo_array[i],
                            ignore_cmd,
                            phenotype_name,
                            phenotype_name,
                            phenotype_name,
                            phenotype_name,
                            phenotype_name)
            print("--- phenotype {} COMMAND ---".format(phenotype_name))
            print(call_command)
            pr=subprocess.run(shlex.split(call_command),stdout=PIPE,stderr=subprocess.STDOUT,encoding="utf8")
            print("--- phenotype {} RETURN CODE: {} ---".format(phenotype_name,pr.returncode))
            print("--- phenotype {} STDOUT ---".format(phenotype_name))
            print(pr.stdout)
        CODE
    >>>
    #output
    output {
        Array[File] out=glob("*.out")
    }
    runtime {
        docker: "${docker}"
        cpu: "${cpus}"
        memory: "${memory} GB"
        disks: "local-disk 300 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }
}


workflow autoreporting{
    File input_array_file
    Int phenos_per_worker
    Int memory
    String docker 
    Int cpus
    Array[String] empty_credset=[]
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
    Float strict_group_r2
    Boolean group
    Boolean overlap
    Boolean check_for_ld
    String additional_prefix
    String primary_grouping_method
    String secondary_grouping_method

    call preprocess_serial{
        input: input_array = input_array_file, phenos_per_worker = phenos_per_worker, docker=docker
    }

    Array[Array[File]] pheno_arr = read_tsv(preprocess_serial.summ_array)
    Array[Array[File]] pheno_tbi_arr = read_tsv(preprocess_serial.summ_tb_array)
    Array[Array[File]] credset_arr = read_tsv(preprocess_serial.credset_array)
    Array[Array[String]] pheno_ids = read_tsv(preprocess_serial.pheno_array)
    Map[String,Boolean] credset_available = read_map(preprocess_serial.credset_map)

    #reports
    scatter (i in range(length(pheno_arr))) {
        call report{
            input: 
            pheno_ids=pheno_ids[i],
            summ_stat=pheno_arr[i], 
            summ_stat_tb=pheno_tbi_arr[i], 
            credsets=credset_arr[i], 
            credset_status = credset_available,
            docker=docker, 
            memory=memory, 
            cpus=cpus, 
            gnomad_exome=gnomad_exome,
            gnomad_genome=gnomad_genome, 
            ld_panel=ld_panel, 
            finngen_annotation=finngen_annotation, 
            functional_annotation=functional_annotation, 
            local_gwcatalog=local_gwcatalog, 
            annotation_version=annotation_version, 
            compare_style=compare_style, 
            db_choice=db_choice, 
            include_batch_freq=include_batch_freq, 
            ignore_region=ignore_region, 
            ld_treshold=ld_treshold, 
            sign_treshold=sign_treshold, 
            alt_sign_treshold=alt_sign_treshold, 
            grouping_locus_width=grouping_locus_width, 
            ld_r2=ld_r2, 
            strict_group_r2=strict_group_r2,
            plink_memory=plink_memory, 
            gwascatalog_pval=gwascatalog_pval, 
            gwascatalog_width_kb=gwascatalog_width_kb, 
            gwascatalog_threads=gwascatalog_threads, 
            group=group, 
            overlap=overlap, 
            check_for_ld=check_for_ld,
            efo_map=efo_code_file,
            primary_grouping_method=primary_grouping_method,
            secondary_grouping_method=secondary_grouping_method,
            additional_prefix=additional_prefix
        }
    }
}