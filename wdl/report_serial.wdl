task report_credible_set {
    #variables
    String docker
    Int memory
    Array[File] summ_stat
    Array[File] summ_stat_tb
    Array[File] credsets
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

    Boolean group
    Boolean overlap
    Boolean include_batch_freq
    Boolean check_for_ld

    String grouping_method
    String ignore_region
    String ignore_cmd = if ignore_region != "" then "--ignore-region" else ""
    String efo_codes
    String compare_style
    String db_choice
    String annotation_version
    String summary_cmd=if defined(ext_summary_stats) then "--summary-fpath" else ""
    String efo_cmd = if efo_codes != "" then "--efo-codes" else ""
    String dollar = "$"  
    #command
    command <<<
        mod_ld=$( echo ${ld_panel_bed} | sed 's/.bed//g' )
        arr_s=(${sep=' ' summ_stat} )
        arr_c=(${sep=' ' credsets} )
        for ((i=0;i<${dollar}{#arr_s[@]};++i)) ; do
            output_filename=$(basename ${dollar}{arr_s[i]})
            main.py ${dollar}{arr_s[i]} --sign-treshold ${sign_treshold} --alt-sign-treshold ${alt_sign_treshold}  \
            ${true='--group' false='' group} --grouping-method ${grouping_method} --locus-width-kb ${grouping_locus_width} \
            --ld-panel-path ${dollar}mod_ld --ld-r2 ${ld_r2} --plink-memory ${plink_memory} ${true='--overlap' false='' overlap} \
            ${ignore_cmd} ${ignore_region}\
            --gnomad-genome-path ${gnomad_genome} --gnomad-exome-path ${gnomad_exome} ${true='--include-batch-freq' false='' include_batch_freq} --finngen-path ${finngen_annotation} \
            --functional-path ${ functional_annotation} --credible-set-file ${dollar}{arr_c[i]} --finngen-annotation-version ${annotation_version} \
            --finngen-annotation-version ${annotation_version}\
            --compare-style ${compare_style} ${true='--check-for-ld' false='' check_for_ld} --ld-treshold ${ld_treshold}  \
            --ldstore-threads ${cpus} --gwascatalog-threads ${gwascatalog_threads} \
            ${summary_cmd} ${write_lines(ext_summary_stats)} ${"--endpoint-fpath " + endpoint_listing} \
            --gwascatalog-pval ${gwascatalog_pval} --gwascatalog-width-kb ${gwascatalog_width_kb} ${efo_cmd} ${efo_codes} --db ${db_choice} ${"--local-gwascatalog " + local_gwcatalog} \
            --fetch-out $output_filename.fetch.out --annotate-out $output_filename.annotate.out --report-out $output_filename.report.out --top-report-out $output_filename.top.out --ld-report-out $output_filename.ld.out
        done  
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

task report_nocred {
    #variables
    String docker
    Int memory
    Array[File] summ_stat
    Array[File] summ_stat_tb
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

    Boolean group
    Boolean overlap
    Boolean include_batch_freq
    Boolean check_for_ld

    String grouping_method
    String grouping_method_ = if grouping_method == "cred" then "ld" else grouping_method #replace grouping method cred with ld if necessary.
    String ignore_region
    String ignore_cmd = if ignore_region != "" then "--ignore-region" else ""
    String efo_codes
    String compare_style
    String db_choice
    String annotation_version
    String summary_cmd=if defined(ext_summary_stats) then "--summary-fpath" else ""
    String efo_cmd = if efo_codes != "" then "--efo-codes" else ""
    String dollar = "$"  
    #command
    command <<<
        mod_ld=$( echo ${ld_panel_bed} | sed 's/.bed//g' )
        arr_s=(${sep=' ' summ_stat} )
        for ((i=0;i<${dollar}{#arr_s[@]};++i)) ; do
            output_filename=$(basename ${dollar}{arr_s[i]})
            main.py ${dollar}{arr_s[i]} --sign-treshold ${sign_treshold} --alt-sign-treshold ${alt_sign_treshold}  \
            ${true='--group' false='' group} --grouping-method ${grouping_method_} --locus-width-kb ${grouping_locus_width} \
            --ld-panel-path ${dollar}mod_ld --ld-r2 ${ld_r2} --plink-memory ${plink_memory} ${true='--overlap' false='' overlap} \
            ${ignore_cmd} ${ignore_region}\
            --gnomad-genome-path ${gnomad_genome} --gnomad-exome-path ${gnomad_exome} ${true='--include-batch-freq' false='' include_batch_freq} --finngen-path ${finngen_annotation} \
            --functional-path ${ functional_annotation} --finngen-annotation-version ${annotation_version} \
            --finngen-annotation-version ${annotation_version}\
            --compare-style ${compare_style} ${true='--check-for-ld' false='' check_for_ld} --ld-treshold ${ld_treshold}  \
            --ldstore-threads ${cpus} --gwascatalog-threads ${gwascatalog_threads} \
            ${summary_cmd} ${write_lines(ext_summary_stats)} ${"--endpoint-fpath " + endpoint_listing} \
            --gwascatalog-pval ${gwascatalog_pval} --gwascatalog-width-kb ${gwascatalog_width_kb} ${efo_cmd} ${efo_codes} --db ${db_choice} ${"--local-gwascatalog " + local_gwcatalog} \
            --fetch-out $output_filename.fetch.out --annotate-out $output_filename.annotate.out --report-out $output_filename.report.out --top-report-out $output_filename.top.out --ld-report-out $output_filename.ld.out
        done  
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
    Int phenos_per_worker
    Int memory
    String docker 
    call credset_filter {
        input:additional_prefix="finngen_R4_",phenotypelist=phenotypelist,credsetlist=credsetlist,docker=docker
    }
    #common pheno chunks
    call simple_preprocess_chunks as common_phenos {
        input: phenotypelist=credset_filter.filtered_pheno_list, num=phenos_per_worker,docker=docker
    }
    #common credsets, tbi files are not valid but no matter
    call simple_preprocess_chunks as common_creds {
        input: phenotypelist=credset_filter.filtered_cred_list, num=phenos_per_worker,docker=docker
    }
    #phenos without credsets
    call simple_preprocess_chunks as other_phenos {
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
            input: summ_stat=c_pheno_arr[i], summ_stat_tb=c_pheno_tbi_arr[i], credsets=credset_arr[i], docker=docker,memory=memory
        }
    }
    #reports with no credible sets
    scatter (i in range(length(o_pheno_arr))) {
        call report_nocred{
            input: summ_stat=o_pheno_arr[i], summ_stat_tb=o_pheno_tbi_arr[i], docker=docker,memory=memory
        }
    }
}