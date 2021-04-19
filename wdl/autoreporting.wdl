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
    Int memory
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
        phenotype_info = "--pheno-info-file ${phenotype_info_file}" if "${phenotype_info_file}" != empty_file else "" 

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

        call_command=(f"main.py {summstat} "
                    f"--pheno-name {pheno_id} "
                    f" --sign-treshold {sign_treshold} " 
                    f"--alt-sign-treshold {alt_sign_treshold} "
                    f"{group} "
                    f"--grouping-method {grouping_method} "
                    f"--locus-width-kb {grouping_locus_width} "
                    f"--ld-panel-path {plink_path} "
                    f"{ld_opts} "#ld opts
                    f"--plink-memory {plink_memory} "
                    f"{include_batch_freq} " #include batch freq
                    f"--finngen-path {finngen_annotation} "
                    f"--functional-path {functional_annotation} "
                    f"--gnomad-genome-path {gnomad_genome} "
                    f"--gnomad-exome-path {gnomad_exome} "
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
                    f"--fetch-out {pheno_id}.fetch.out "
                    f"--annotate-out {pheno_id}.annotate.out "
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
    >>>

    output {
        Array[File] out=glob("*.out")
    }
    runtime {
        docker: "${docker}"
        cpu: "${cpus}"
        memory: "${memory} GB"
        disks: "local-disk 100 HDD"
        zones: "europe-west1-b"
        preemptible: 2
    }
}

workflow autoreporting{
    File input_array_file
    Array[Array[String]] input_array = read_tsv(input_array_file)

    scatter (arr in  input_array ){
        call report {
            input: input_file_list = arr
        }
    }

}
