import "report_utils.wdl" as utils

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
    String real_credset = if length(select_all([credible_set]))==0 then "False" else "True"
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
    String efo_codes = select_first([efo_map[map_pheno],""])
    String compare_style
    String db_choice
    String summary_cmd=if defined(ext_summary_stats) then "--summary-fpath" else ""
    String efo_cmd = if efo_codes != "" then "--efo-codes" else ""
    String dollar = "$"  
    String annotation_version
    Float strict_group_r2


    command <<<
        python3 <<CODE
        import subprocess, shlex
        from subprocess import PIPE
        #do everything a wrapper should do
        #unchanged variables
        plink_path = "${ld_panel_bed}".replace(".bed","")
        gnomad_exome="${gnomad_exome}"
        gnomad_genome="${gnomad_genome}"
        finngen_annotation="${finngen_annotation}"
        functional_annotation="${functional_annotation}"
        ext_summary_stats="${write_lines(ext_summary_stats)}"
        endpoint_listing="${"--endpoint-fpath " + endpoint_listing}"
        local_gwascatalog="${"--local-gwascatalog "+ local_gwcatalog}"

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
        summstat="${summ_stat}"
        #credible set
        credset=""
        if ${real_credset} == True:
            credset="--credible-set-file ${select_first([credible_set,""])}"
        #efo codes
        efo_cmd = ""
        if "${efo_codes}".replace("\"","") != "":
            efo_cmd= "--efo-codes ${efo_codes}"


        phenotype_name=summstat.split("/")[-1].split(".")[0].replace("finngen_R4_","")
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
                    ).format(summstat,
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
                        credset,
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
                        efo_cmd,
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
    Float strict_group_r2

    call utils.fill_efo_map as fill_efo_map{
        input:additional_prefix="finngen_R4_",docker=docker,phenotypelist=phenotypelist,efomap=efo_code_file
    }

    call utils.credset_filter as credset_filter {
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
            gnomad_genome=gnomad_genome, ld_panel=ld_panel, strict_group_r2=strict_group_r2,
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
            gnomad_exome=gnomad_exome,gnomad_genome=gnomad_genome, ld_panel=ld_panel, strict_group_r2=strict_group_r2, 
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
