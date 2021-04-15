task report {
    #variables
    String docker
    Array[String] input_file_list
    Int arr_len = length(input_file_list)
    String phenotype_name = input_file_list[0]
    File summ_stat = input_file_list[1]
    File summ_stat_tb=summ_stat+".tbi"
    File credible_set = input_file_list[2]
    File credible_set_cred = sub(sub(input_file_list[2],".snp.filter",".cred.summary"),".snp",".cred") #if there is snp.filter, it gets subbed to .cred.summary, and the latter sub doesn't do anything, otherwise the latter sub does the sub. clunky but should work.
    File previous_release = input_file_list[3]
    File previous_release_tbi =previous_release+".tbi" 

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

    File? local_gwcatalog

    Float sign_treshold
    Float alt_sign_treshold
    Int grouping_locus_width
    String ld_opts
    Int plink_memory
    Float gwascatalog_pval
    Int gwascatalog_width_kb
    Int cpus
    Int gwascatalog_threads
    Int docker_memory
    File efo_map
    File custom_dataresource

    Boolean group
    Boolean overlap
    Boolean include_batch_freq

    String primary_grouping_method
    String secondary_grouping_method
    String ignore_region

    String db_choice
    Array[String] column_names
    String extra_columns
    Float strict_group_r2
    String phenoname = basename(phenotype_name,".gz")
    File dummy_file
    File phenotype_info_file

    File allele_vcf_file
    File allele_vcf_tbi = allele_vcf_file+".tbi"

    command <<<
        python3 <<CODE
        import subprocess, shlex, sys
        from subprocess import PIPE
        #do everything a wrapper should do
        #unchanged variables
        empty_file = "${dummy_file}"
        pheno_id="${phenotype_name}"
        plink_path = "${ld_panel_bed}".replace(".bed","")
        gnomad_exome="${gnomad_exome}"
        gnomad_genome="${gnomad_genome}"
        finngen_annotation="${finngen_annotation}"
        functional_annotation="${functional_annotation}"
        previous_release="--previous-release-path ${previous_release}" if "${previous_release}" != empty_file else "" 
        local_gwascatalog="${"--local-gwascatalog "+ local_gwcatalog}"

        sign_treshold=${sign_treshold}
        alt_sign_treshold=${alt_sign_treshold}
        grouping_locus_width=${grouping_locus_width}
        ld_opts="${ld_opts}"
        plink_memory=${plink_memory}
        gwascatalog_pval=${gwascatalog_pval}
        gwascatalog_width_kb=${gwascatalog_width_kb}
        gwascatalog_threads=${gwascatalog_threads}
        strict_group_r2=${strict_group_r2}

        group="--group" if "${group}"=="true" else ""
        overlap="--overlap" if "${overlap}"=="true" else ""
        include_batch_freq="--include-batch-freq" if "${include_batch_freq}"=="true" else ""
        grouping_method = "${primary_grouping_method}" if "${credible_set}" != empty_file else "${secondary_grouping_method}" 
        ignore_cmd = "--ignore-region ${ignore_region}" if "${ignore_region}" != "" else ""
        db_choice = "${db_choice}"
        custom_dataresource="${custom_dataresource}"
        column_names = "${sep=" " column_names}"
        extra_columns = "${extra_columns}"
        phenotype_info = "${phenotype_info_file}"

        alleledb_file = "${allele_vcf_file}"

        #changing variables
        #summ stat
        summstat="${summ_stat}"
        #credible set
        credset=""
        if "${credible_set}" != empty_file:
            credset="--credible-set-file ${credible_set}"

        #efo codes
        with open("${efo_map}","r") as f:
            efos = {a.strip().split("\t")[0] : a.strip().split("\t")[1]  for a in f.readlines()}
        efo_cmd=""
        if pheno_id in efos.keys():
            if efos[pheno_id] != "NA" and efos[pheno_id] != "":
                efo_cmd="--efo-codes {}".format(efos[pheno_id])

        call_command=("main.py {} "
                    "--pheno-name {} "
                    " --sign-treshold {} " 
                    "--alt-sign-treshold {} "
                    "{} "
                    "--grouping-method {} "
                    "--locus-width-kb {} "
                    "--ld-panel-path {} "
                    "{} "#ld opts
                    "--plink-memory {} "
                    "{} " #include batch freq
                    "--finngen-path {} "
                    "--functional-path {} "
                    "--gnomad-genome-path {} "
                    "--gnomad-exome-path {} "
                    "{} "
                    "{} "
                    "--use-gwascatalog "
                    "--gwascatalog-threads {} "
                    "--strict-group-r2 {} "
                    "--gwascatalog-pval {} "
                    "--gwascatalog-width-kb {} "
                    "--db {} "
                    "--column-labels {} "
                    "--extra-cols {} "
                    "--pheno-info-file {} "
                    "--gwascatalog-allele-file {}"
                    "{} " #local gwascatalog
                    "{} " #efo
                    "{} " #ignore
                    "--custom-dataresource {} "
                    "--fetch-out {}.fetch.out "
                    "--annotate-out {}.annotate.out "
                    "--report-out {}.report.out "
                    "--top-report-out {}.top.out "
                    ).format(summstat,
                        pheno_id,
                        sign_treshold,
                        alt_sign_treshold,
                        group,
                        grouping_method,
                        grouping_locus_width,
                        plink_path,
                        ld_opts,
                        plink_memory,
                        include_batch_freq,
                        finngen_annotation,
                        functional_annotation,
                        gnomad_genome,
                        gnomad_exome,
                        previous_release,
                        credset,
                        gwascatalog_threads,
                        strict_group_r2,
                        gwascatalog_pval,
                        gwascatalog_width_kb,
                        db_choice,
                        column_names,
                        extra_columns,
                        phenotype_info,
                        alleledb_file,
                        local_gwascatalog,
                        efo_cmd,
                        ignore_cmd,
                        custom_dataresource,
                        pheno_id,
                        pheno_id,
                        pheno_id,
                        pheno_id)
        print("--- phenotype {} COMMAND ---".format(pheno_id))
        print(call_command)
        pr=subprocess.run(shlex.split(call_command),stdout=PIPE,stderr=subprocess.STDOUT,encoding="utf8")
        print("--- phenotype {} RETURN CODE: {} ---".format(pheno_id,pr.returncode))
        print("--- phenotype {} STDOUT ---".format(pheno_id))
        print(pr.stdout)
        if pr.returncode != 0:
            print("The report did not run successfully. Check the logs.")
            print("${phenoname}","exit code:",pr.returncode)
            sys.exit(1)
        CODE
    >>>

    output {
        Array[File] out=glob("*.out")
    }
    runtime {
        docker: "${docker}"
        cpu: "${cpus}"
        memory: "${docker_memory} GB"
        disks: "local-disk 100 HDD"
        zones: "europe-west1-b"
        preemptible: 2
    }
}

workflow autoreporting{
    File input_array_file
    Array[Array[String]] input_array = read_tsv(input_array_file)
    String docker
    Int memory
    Int cpus

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
    String ld_opts
    Int plink_memory
    Float gwascatalog_pval
    Int gwascatalog_width_kb
    Int gwascatalog_threads
    Boolean group
    Boolean overlap
    Float strict_group_r2
    String primary_grouping_method
    String secondary_grouping_method
    File custom_dataresource
    Array[String] column_names
    String extra_columns
    File dummy_file
    File phenotype_info_file

    scatter (arr in  input_array ){
        call report {
            input: input_file_list = arr,
            docker=docker,
            primary_grouping_method=primary_grouping_method,
            secondary_grouping_method=secondary_grouping_method,
            docker_memory=memory,
            cpus=cpus,
            gnomad_exome=gnomad_exome,
            gnomad_genome=gnomad_genome,
            ld_panel=ld_panel,
            strict_group_r2=strict_group_r2,
            finngen_annotation=finngen_annotation,
            functional_annotation=functional_annotation,
            local_gwcatalog=local_gwcatalog,
            db_choice=db_choice,
            efo_map=efo_code_file,
            include_batch_freq=include_batch_freq,
            ignore_region=ignore_region,
            sign_treshold=sign_treshold,
            alt_sign_treshold=alt_sign_treshold,
            grouping_locus_width=grouping_locus_width,
            ld_opts=ld_opts,
            plink_memory=plink_memory,
            gwascatalog_pval=gwascatalog_pval,
            gwascatalog_width_kb=gwascatalog_width_kb,
            gwascatalog_threads=gwascatalog_threads,
            group=group,
            overlap=overlap,
            custom_dataresource=custom_dataresource,
            column_names=column_names,
            extra_columns=extra_columns,
            dummy_file=dummy_file,
            phenotype_info_file=phenotype_info_file
        }
    }

}
