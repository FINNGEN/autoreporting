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
    }
}


task report{
    #variables
    String docker
    Int memory
    Array[String] pheno_ids
    Array[File] summ_stat
    Array[File] summ_stat_tb
    Array[String] credsets
    Boolean use_credsets = if length(credsets) > 0 then true else false
    Array[File] selected_credsets = if use_credsets then credsets else []
    Int len = length(pheno_ids)

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
    File custom_dataresource

    File? local_gwcatalog

    Float sign_treshold
    Float alt_sign_treshold
    Int grouping_locus_width
    Float ld_r2
    Int plink_memory
    Float gwascatalog_pval
    Int gwascatalog_width_kb
    Int cpus
    Int gwascatalog_threads
    Float strict_group_r2

    Boolean group
    Boolean overlap
    Boolean include_batch_freq

    String primary_grouping_method
    String secondary_grouping_method
    String ignore_region
    String db_choice
    Array[String] column_names
    String extra_columns

    command <<<
        python3 <<CODE
        import subprocess, shlex, json, sys, itertools
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

        sign_treshold=${sign_treshold}
        alt_sign_treshold=${alt_sign_treshold}
        grouping_locus_width=${grouping_locus_width}
        ld_r2=${ld_r2}
        plink_memory=${plink_memory}
        gwascatalog_pval=${gwascatalog_pval}
        gwascatalog_width_kb=${gwascatalog_width_kb}
        gwascatalog_threads=${gwascatalog_threads}
        strict_group_r2=${strict_group_r2}

        group="--group" if "${group}"=="true" else ""
        overlap="--overlap" if "${overlap}"=="true" else ""
        include_batch_freq="--include-batch-freq" if "${include_batch_freq}"=="true" else ""
        ignore_cmd = "--ignore-region ${ignore_region}" if "${ignore_region}" != "" else ""
        db_choice = "${db_choice}"
        grouping_method = "${primary_grouping_method}" if "true" == "${use_credsets}" else "${secondary_grouping_method}"
        #process the summstats and credsets etc
        summstats="${sep=";" summ_stat}".split(";")
        phenotypes="${sep=";" pheno_ids}".split(";")
        custom_dataresource="${custom_dataresource}"
        column_names = "${sep=" " column_names}"
        extra_columns = "${extra_columns}"
        credset_calls = []
        credset_cmds=["--credible-set-file {}".format(a) for a in "${sep=";" selected_credsets}".split(";")] if "true" == "${use_credsets}" else [""] * ${len}
        #efo codes
        with open("${efo_map}","r") as f:
            efos = {a.strip().split("\t")[0] : a.strip().split("\t")[1]  for a in f.readlines()}
        efo_array=["--efo-codes {}".format(efos[a]) if a in efos.keys() else "" for a in phenotypes ]
        
        exit_codes = []
        for i in range(num_phenos):
            phenotype_name=phenotypes[i]
            call_command=("main.py {} "
                        "--pheno-name {} "
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
                        "--use-gwascatalog "
                        "--gwascatalog-threads {} "
                        "--strict-group-r2 {} "
                        "--gwascatalog-pval {} "
                        "--gwascatalog-width-kb {} "
                        "--db {} "
                        "--column-labels {} "
                        "--extra-cols {} "
                        "{} "
                        "{} "
                        "{} "
                        "--custom-dataresource {} "
                        "--fetch-out {}.fetch.out "
                        "--annotate-out {}.annotate.out "
                        "--report-out {}.report.out "
                        "--top-report-out {}.top.out "
                        ).format(summstats[i],
                            phenotype_name,
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
                            credset_cmds[i],
                            gwascatalog_threads,
                            strict_group_r2,
                            gwascatalog_pval,
                            gwascatalog_width_kb,
                            db_choice,
                            column_names,
                            extra_columns,
                            local_gwascatalog,
                            efo_array[i],
                            ignore_cmd,
                            custom_dataresource,
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
            exit_codes.append(pr.returncode)
        if any([a for a in exit_codes if a != 0]):
            print("One or more of the reports did not run successfully. Check the logs.")
            print("Exit codes:")
            print(list(itertools.zip_longest(phenotypes, exit_codes) ) )
            sys.exit(1)
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
        disks: "local-disk 200 HDD"
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
    String db_choice
    File efo_code_file
    String ignore_region
    Boolean include_batch_freq
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
    String primary_grouping_method
    String secondary_grouping_method
    File custom_dataresource
    Array[String] column_names
    String extra_columns

    call preprocess_serial{
        input: input_array = input_array_file, phenos_per_worker = phenos_per_worker, docker=docker
    }

    Array[Array[File]] pheno_arr = read_tsv(preprocess_serial.summ_array)
    Array[Array[File]] pheno_tbi_arr = read_tsv(preprocess_serial.summ_tb_array)
    Array[Array[File]] credset_arr = read_tsv(preprocess_serial.credset_array)
    Array[Array[String]] pheno_ids = read_tsv(preprocess_serial.pheno_array)

    #reports
    scatter (i in range(length(pheno_arr))) {
        Array[String] credset_input = if length(credset_arr[i]) == length(pheno_arr[i]) then credset_arr[i]  else []
        call report{
            input: 
            pheno_ids=pheno_ids[i],
            summ_stat=pheno_arr[i], 
            summ_stat_tb=pheno_tbi_arr[i], 
            credsets=credset_input, 
            docker=docker, 
            memory=memory, 
            cpus=cpus, 
            gnomad_exome=gnomad_exome,
            gnomad_genome=gnomad_genome, 
            ld_panel=ld_panel, 
            finngen_annotation=finngen_annotation, 
            functional_annotation=functional_annotation, 
            local_gwcatalog=local_gwcatalog, 
            db_choice=db_choice, 
            include_batch_freq=include_batch_freq, 
            ignore_region=ignore_region, 
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
            efo_map=efo_code_file,
            primary_grouping_method=primary_grouping_method,
            secondary_grouping_method=secondary_grouping_method,
            custom_dataresource=custom_dataresource,
            column_names=column_names,
            extra_columns=extra_columns
        }
    }
}