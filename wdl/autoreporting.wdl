task report {
    #variables
    String docker
    File summ_stat
    File summ_stat_tb=summ_stat+".tbi"
    File gnomad_genome
    File gnomad_genome_tb=gnomad_genome+".tbi"
    File gnomad_exome
    File gnomad_exome_tb=gnomad_exome+".tbi"
    File finngen_annotation
    File finngen_annotation_tb=finngen_annotation+".tbi"
    File ld_chrom_input_file
    Array[File] ld_chrom_files = read_lines(ld_chrom_input_file)
    
    String ld_panel
    File ld_panel_bed=ld_panel+".bed"
    File ld_panel_bim=ld_panel+".bim"
    File ld_panel_fam=ld_panel+".fam"
    Float s_tresh
    Float s_tresh2
    Boolean group
    String gr_method
    Int grouping_locus_width
    Float ld_r2
    Int plink_mem
    Boolean overlap
    Boolean batch_freq
    String compare_style
    Boolean check_for_ld
    Float gw_pval
    Int gw_width
    Float ld_treshold
    String dollar = "$"
    #trick wdl to write the external summary stats paths as a file
    File? summary_stat_listing
    File? endpoint_listing
    Array[File]? ext_summary_stats = if defined(summary_stat_listing) then read_lines(summary_stat_listing  ) else []
    String summary_cmd=if defined(ext_summary_stats) then "--summary-fpath" else ""
    String efo_codes
    String efo_cmd = if efo_codes != "" then "--efo-codes" else ""
    String db_choice
    File? local_gwcatalog

    #output file names
    String fetch_out
    String annotate_out
    String raport_out
    String top_raport_out
    String ld_raport_out

    #command
    command {
        mod_ld=$( echo ${ld_panel_bed} | sed 's/.bed//g' )
        ld_chrom=$( echo ${ld_chrom_files[0]} | sed 's/_[0-9].[fb][aei][dm]//g' )

        main.py ${summ_stat} --sign-treshold ${s_tresh} --alt-sign-treshold ${s_tresh2}  \
        ${true='--group' false='' group} --grouping-method ${gr_method} --locus-width-kb ${grouping_locus_width} \
        --ld-panel-path ${dollar}mod_ld --ld-r2 ${ld_r2} --plink-memory ${plink_mem} ${true='--overlap' false='' overlap} \
        --gnomad-genome-path ${gnomad_genome} --gnomad-exome-path ${gnomad_exome} ${true='--include-batch-freq' false='' batch_freq} --finngen-path ${finngen_annotation} \
        --compare-style ${compare_style} ${true='--check-for-ld' false='' check_for_ld} --ld-treshold ${ld_treshold}  \
        ${summary_cmd} ${write_lines(ext_summary_stats)} ${"--endpoint-fpath " + endpoint_listing} \
        --ld-chromosome-panel-path ${dollar}ld_chrom \
        --gwascatalog-pval ${gw_pval} --gwascatalog-width-kb ${gw_width} ${efo_cmd} ${efo_codes} --db ${db_choice} ${"--local-gwascatalog " + local_gwcatalog} \
        --fetch-out ${fetch_out} --annotate-out ${annotate_out} --raport-out ${raport_out} --top-report-out ${top_raport_out} 
    }
    #output
    output {
        File gws_out_ = summ_stat+fetch_out
        File annotate_out_ = summ_stat+annotate_out
        File raport_out_ = summ_stat+raport_out
        File top_raport_out_ = summ_stat+top_raport_out
        File ld_out = summ_stat+ld_raport_out
    }
    runtime {
        docker: "${docker}"
        cpu: 3
        memory: "16 GB"
        disks: "local-disk 300 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }
}

workflow autoreporting{
    File phenotypelist
    Array[File] phenotypes = read_lines(phenotypelist)
    scatter (pheno in phenotypes){
        call report {
            input: summ_stat=pheno
        }
    }
    

}
