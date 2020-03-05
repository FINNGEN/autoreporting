import "report_utils.wdl" as utils

task report_credible_set {
    #variables
    String docker
    Int memory
    Array[File] summ_stat
    Int len = length(summ_stat)
    Array[File] summ_stat_tb
    Array[File?] credsets
    Array[File] selected_credsets=select_all(credsets)
    Boolean empty_credset = length(selected_credsets)==0
    #File credset_file = if defined(credsets) then write_lines(credsets) else write_lines(summ_stat_tb)
    String real_credset = if !empty_credset then "True" else "False"
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

    File? summary_stat_listing
    File? endpoint_listing
    #trick wdl to write the external summary stats paths as a file
    #NOTE: does not work with partially serialized widdle. 
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
    Float strict_group_r2

    Boolean group
    Boolean overlap
    Boolean include_batch_freq
    Boolean check_for_ld

    String grouping_method
    String ignore_region
    String compare_style
    String db_choice
    String annotation_version
    String summary_cmd=if defined(ext_summary_stats) then "--summary-fpath" else ""
    String additional_prefix
    #String efo_cmd = if efo_codes != "" then "--efo-codes" else ""
    String dollar = "$"  

    command <<<
        python3 <<CODE
        import subprocess, shlex
        from subprocess import PIPE
        #do everything a wrapper should do
        #unchanged variables
        num_phenos = ${len}
        plink_path = "${ld_panel_bed}".replace(".bed","")
        gnomad_exome="${gnomad_exome}"
        gnomad_genome="${gnomad_genome}"
        finngen_annotation="${finngen_annotation}"
        functional_annotation="${functional_annotation}"
        ext_summary_stats="${write_lines(ext_summary_stats)}"
        endpoint_listing="${"--endpoint-fpath " + endpoint_listing}"
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
        grouping_method = "${grouping_method}"
        ignore_cmd = "--ignore-region ${ignore_region}" if "${ignore_region}" != "" else ""
        compare_style = "${compare_style}"
        db_choice = "${db_choice}"
        annotation_version = "${annotation_version}"
        summary_cmd="${summary_cmd} ${write_lines(ext_summary_stats)}"

        #changing variables
        #summ stat
        with open("${write_lines( summ_stat )}","r") as f:
            summstats=[l.strip() for l in f]
        #credible set
        credsets=[""]*num_phenos
        if ${real_credset} == True:
            with open("${write_lines(selected_credsets)}","r") as f2:
                credsets=["--credible-set-file {}".format(l.strip()) for l in f2]
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
                        "{} "
                        "{} "
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
                            summary_cmd,
                            endpoint_listing,
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
    File phenotypelist
    File credsetlist
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

    call utils.credset_filter as credset_filter {
        input:additional_prefix=additional_prefix,phenotypelist=phenotypelist,credsetlist=credsetlist,docker=docker
    }

    #common pheno chunks
    call utils.simple_preprocess_chunks as common_phenos {
        input: phenotypelist=credset_filter.filtered_pheno_list, num=phenos_per_worker,docker=docker
    }
    #common credsets, tbi files are not valid but no matter
    call utils.simple_preprocess_chunks as common_creds {
        input: phenotypelist=credset_filter.filtered_cred_list, num=phenos_per_worker,docker=docker
    }

    call utils.simple_preprocess_chunks as other_phenos {
        input: phenotypelist=credset_filter.other_pheno_list, num=phenos_per_worker,docker=docker
    }
    Array[Array[File]] c_pheno_arr = read_tsv(common_phenos.out_tsv)
    Array[Array[File]] c_pheno_tbi_arr = read_tsv(common_phenos.out_tbi_tsv)

    Array[Array[File]] credset_arr = read_tsv(common_creds.out_tsv)


    Array[Array[File]] o_pheno_arr = read_tsv(other_phenos.out_tsv)
    Array[Array[File]] o_pheno_tbi_arr = read_tsv(other_phenos.out_tbi_tsv)
    #credset reports
    scatter (i in range(length(c_pheno_arr))) {
        call report_credible_set{
            input: summ_stat=c_pheno_arr[i], 
            summ_stat_tb=c_pheno_tbi_arr[i], 
            credsets=credset_arr[i], 
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
            additional_prefix=additional_prefix
        }
    }
    #reports with no credible sets
    scatter (i in range(length(o_pheno_arr))) {
        call report_credible_set as report_nocred{
            input: summ_stat=o_pheno_arr[i], 
            summ_stat_tb=o_pheno_tbi_arr[i], 
            credsets=empty_credset,
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
            additional_prefix=additional_prefix
        }
    }
}