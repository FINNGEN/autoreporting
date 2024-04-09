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

    File gnomad
    File gnomad_tbi = gnomad+".tbi"
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
    Int locus_width_kb
    String ld_opts
    Int plink_memory
    Float gwascatalog_pval
    Int gwascatalog_width_kb
    Int cpus
    Int gwascatalog_threads
    Int memory
    File efo_map
    File custom_dataresource

    Boolean group
    Boolean overlap

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

    Float disk_size_f = size(summ_stat,"G")+
        size(previous_release,"G")+
        size(gnomad,"G")+
        size(finngen_annotation,"G")+
        size(functional_annotation,"G")+
        size(ld_panel_bed,"G")+
        size(allele_vcf_file,"G")
    Int disk_size = ceil(disk_size_f*1.25)+25

    command <<<
        set -euxo pipefail
        python3 <<CODE
        import subprocess, shlex, sys
        from subprocess import PIPE
        #do everything a wrapper should do
        #unchanged variables
        empty_file = "${dummy_file}"
        pheno_id="${phenotype_name}"
        plink_path = "${ld_panel_bed}".replace(".bed","")
        gnomad="--gnomad-path ${gnomad}" if "${gnomad}" != empty_file else ""
        finngen_annotation="--finngen-path  ${finngen_annotation}" if "${finngen_annotation}" != empty_file else "" 
        functional_annotation="--functional-path ${functional_annotation}" if "${functional_annotation}" != empty_file else "" 
        previous_release="--previous-release-path ${previous_release}" if "${previous_release}" != empty_file else "" 
        local_gwascatalog="${"--local-gwascatalog "+ local_gwcatalog}"

        sign_treshold=${sign_treshold}
        alt_sign_treshold=${alt_sign_treshold}
        locus_width_kb=${locus_width_kb}
        ld_opts="${ld_opts}"
        plink_memory=${plink_memory}
        gwascatalog_pval=${gwascatalog_pval}
        gwascatalog_width_kb=${gwascatalog_width_kb}
        gwascatalog_threads=${gwascatalog_threads}
        strict_group_r2=${strict_group_r2}

        group="--group" if "${group}"=="true" else ""
        overlap="--overlap" if "${overlap}"=="true" else ""
        grouping_method = "${primary_grouping_method}" if "${credible_set}" != empty_file else "${secondary_grouping_method}" 
        ignore_cmd = "--ignore-region ${ignore_region}" if "${ignore_region}" != "" else ""
        db_choice = "${db_choice}"
        custom_dataresource="${custom_dataresource}"
        column_names = "${sep=" " column_names}"
        extra_columns = "${extra_columns}"
        phenotype_info = "--pheno-info-file ${phenotype_info_file}" if "${phenotype_info_file}" != empty_file else "" 

        alleledb_file = "${allele_vcf_file}"

        #changing variables
        #summ stat
        summstat="${summ_stat}"
        #credible set
        credset=""
        if "${credible_set}" != empty_file:
            credset="--credible-set-file ${credible_set}"

        
        efofile = "${efo_map}"
        #efo codes
        efo_cmd=""
        if efofile != empty_file:
            with open(efofile,"r") as f:
                efos = {a.strip().split("\t")[0] : a.strip().split("\t")[1]  for a in f.readlines()}
            if pheno_id in efos.keys():
                if efos[pheno_id] != "NA" and efos[pheno_id] != "":
                    efo_cmd="--efo-codes {}".format(efos[pheno_id])

        call_command=(f"main.py {summstat} "
                    f"--pheno-name {pheno_id} "
                    f" --sign-treshold {sign_treshold} " 
                    f"--alt-sign-treshold {alt_sign_treshold} "
                    f"{group} "
                    f"--grouping-method {grouping_method} "
                    f"--locus-width-kb {locus_width_kb} "
                    f"--ld-panel-path {plink_path} "
                    f"{ld_opts} "#ld opts
                    f"--plink-memory {plink_memory} "
                    f"{finngen_annotation} "
                    f"{functional_annotation} "
                    f"{gnomad} "
                    f"{previous_release} "
                    f"{credset} "
                    f"--use-gwascatalog "
                    f"--gwascatalog-threads {gwascatalog_threads} "
                    f"--strict-group-r2 {strict_group_r2} "
                    f"--gwascatalog-pval {gwascatalog_pval} "
                    f"--gwascatalog-width-kb {gwascatalog_width_kb} "
                    f"--db {db_choice} "
                    f"--column-labels {column_names} "
                    f"--extra-cols {extra_columns} "
                    f"{phenotype_info} "
                    f"--gwascatalog-allele-file {alleledb_file} "
                    f"{local_gwascatalog} " #local gwascatalog
                    f"{efo_cmd} " #efo
                    f"{ignore_cmd} " #ignore
                    f"--custom-dataresource {custom_dataresource} "
                    f"--report-out {pheno_id}.report.out "
                    f"--top-report-out {pheno_id}.top.out "
                    )
        print(f"--- phenotype {pheno_id} COMMAND ---")
        print(call_command)
        pr=subprocess.run(shlex.split(call_command),stdout=PIPE,stderr=subprocess.STDOUT,encoding="utf8")
        print(f"--- phenotype {pheno_id} RETURN CODE: {pr.returncode} ---")
        print(f"--- phenotype {pheno_id} STDOUT ---")
        print(pr.stdout)
        if pr.returncode != 0:
            print("The report did not run successfully. Check the logs.")
            print("${phenoname}","exit code:",pr.returncode)
            sys.exit(1)
        CODE
        if test -f "${phenotype_name}.report.out"; then
            echo "True" > had_results
        else
            echo "False" > had_results
            touch ${phenotype_name}.report.out
            touch ${phenotype_name}.top.out
        fi
    >>>

    output {
        Boolean had_results = read_boolean("had_results")
        File variant_report = phenotype_name+".report.out" 
        File group_report = phenotype_name+".top.out"
    }
    runtime {
        docker: "${docker}"
        cpu: "${cpus}"
        memory: "${memory} GB"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
    }
}

task post_process_top_reports {
    #from workflow
    File top_report
    String docker #
    String ld_panel
    Int locus_width_kb

    #directly to this task
    Float r2_threshold
    String af_col
    Float af_threshold
    String in_fg_col
    
    #derived
    File plink_bed = ld_panel +".bed"
    File plink_bim = ld_panel +".bim"
    File plink_fam = ld_panel +".fam"
    String top = basename(top_report,".top.out")
    String dollar="$"
    command <<<
        #filter the top report
        meta_filter_top.py ${top_report} --width-kb ${locus_width_kb} --r2-threshold ${r2_threshold} --af-column ${af_col} \
            --af-threshold ${af_threshold} --fg-specific-column ${in_fg_col} --output ${top}.top
        #check r2 between hits
        bed_path="${plink_bed}"
        plink_path="${dollar}{bed_path%.*}"
        post_process_hits.py  ${top}.top --ld-panel-path $plink_path  --region-width-kb ${locus_width_kb} --output ${top}.top.post.tsv
    >>>
    output {
        File top_filtered = top+".top"
        File not_in_fg = top+".top.not_in_finngen"
        File removed_hits = top+".top.filtered"
        File post_processed = top+".top.post.tsv"
    }
    runtime {
        docker: "${docker}"
        cpu: "4"
        memory: "30 GB"
        disks: "local-disk 100 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
    }
}

workflow autoreporting{
    File input_array_file
    String docker
    String ld_panel
    Int locus_width_kb
    Array[Array[String]] input_array = read_tsv(input_array_file)

    scatter (arr in  input_array ){
        call report {
            input: input_file_list = arr, docker=docker, ld_panel=ld_panel, locus_width_kb=locus_width_kb
        }
        if( report.had_results ) {
            File variants = report.variant_report
            File groups = report.group_report

            call post_process_top_reports {
                input: top_report = report.group_report, docker=docker, ld_panel=ld_panel,locus_width_kb=locus_width_kb
            }
        }

    }

    output {
        Array[File] variant_report = select_all(variants)
        Array[File] group_report = select_all(groups)
        Array[File] filtered_group = select_all(post_process_top_reports.top_filtered)
        Array[File] not_in_finngen = select_all(post_process_top_reports.not_in_fg)
        Array[File] removed_hits = select_all(post_process_top_reports.removed_hits)
        Array[File] post_processed_group = select_all(post_process_top_reports.post_processed)
    }

}
