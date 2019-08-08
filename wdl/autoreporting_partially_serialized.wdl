task report {
    #variables
    String docker
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
    #File ld_chrom_input_file
    #Array[File] ld_chrom_files = read_lines(ld_chrom_input_file)
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
    String summary_cmd=if defined(ext_summary_stats) then "--summary-fpath" else ""
    String efo_cmd = if efo_codes != "" then "--efo-codes" else ""
    String dollar = "$"  

    #command
    command {
        mod_ld=$( echo ${ld_panel_bed} | sed 's/.bed//g' )
        for file in ${sep=' ' summ_stat} ; do
            output_filename=$(basename $file)
            main.py $file --sign-treshold ${sign_treshold} --alt-sign-treshold ${alt_sign_treshold}  \
            ${true='--group' false='' group} --grouping-method ${grouping_method} --locus-width-kb ${grouping_locus_width} \
            --ld-panel-path ${dollar}mod_ld --ld-r2 ${ld_r2} --plink-memory ${plink_memory} ${true='--overlap' false='' overlap} \
            ${ignore_cmd} ${ignore_region}\
            --gnomad-genome-path ${gnomad_genome} --gnomad-exome-path ${gnomad_exome} ${true='--include-batch-freq' false='' include_batch_freq} --finngen-path ${finngen_annotation} \
            --compare-style ${compare_style} ${true='--check-for-ld' false='' check_for_ld} --ld-treshold ${ld_treshold}  \
            --ldstore-threads ${cpus} --gwascatalog-threads ${gwascatalog_threads} \
            ${summary_cmd} ${write_lines(ext_summary_stats)} ${"--endpoint-fpath " + endpoint_listing} \
            --gwascatalog-pval ${gwascatalog_pval} --gwascatalog-width-kb ${gwascatalog_width_kb} ${efo_cmd} ${efo_codes} --db ${db_choice} ${"--local-gwascatalog " + local_gwcatalog} \
            --fetch-out $output_filename.fetch.out --annotate-out $output_filename.annotate.out --raport-out $output_filename.raport.out --top-report-out $output_filename.top.out --ld-raport-out $output_filename.ld.out
        done  
    }
    #output
    output {
        Array[File] out=glob("*.out")
    }
    runtime {
        docker: "${docker}"
        cpu: "${cpus}"
        memory: "16 GB"
        disks: "local-disk 300 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }
}

task preprocess_chunks{
    File phenotypelist
    Int num
    String docker
    command <<<
        python3 <<CODE
        with open("${phenotypelist}") as f:
            lines = f.readlines()
            lines=[l.strip("\n") for l in lines]
            n=int(${num})
            f2=open("phenofiles.tbi.tsv","w")
            with open("phenofiles.tsv","w") as w:
                for i in range(n):
                    tmparr=lines[i*(len(lines)//n+1):(i+1)*(len(lines)//n+1)]
                    if tmparr==[]:
                        break
                    w.write("\t".join( tmparr ))
                    tbiarr=["".join([a.strip("\n"),".tbi"]) for a in tmparr]
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

workflow autoreporting{
    File phenotypelist
    Int phenos_per_worker
    String docker 
    call preprocess_chunks{
        input: phenotypelist=phenotypelist, num=phenos_per_worker,docker=docker
    }
    Array[Array[File]] pheno_arr = read_tsv(preprocess_chunks.out_tsv)
    Array[Array[File]] pheno_tbi_arr = read_tsv(preprocess_chunks.out_tbi_tsv)
    scatter (i in range(length(pheno_arr))) {
        call report{
            input: summ_stat=pheno_arr[i], summ_stat_tb=pheno_tbi_arr[i], docker=docker
        }
    }
}
