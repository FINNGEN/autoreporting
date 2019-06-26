task report {
    #variables
    String docker
    File summ_stat
    File gnomad_genome
    File gnomad_exome
    File ld_panel
    File finngen_annotation
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
    Array[File] compare_summary_stats
    Array[String] compare_endpoints
    String check_for_ld
    Float gw_pval
    Int gw_width
    Float ld_treshold

    #command
    command {
        python3 /Scripts/main.py ${summ_stat} --sign_treshold ${s_tresh} --alt-sign-treshold ${s_tresh2} ${true='--group' false='' group} \
        --grouping-method ${gr_method} --locus-width-kb ${grouping_locus_width} --ld-panel-path ${ld_panel} \
        --ld-r2 ${ld_r2} --plink-memory ${plink_mem} ${true='--overlap' false='' overlap} --gnomad-genome-path ${gnomad_genome} \
        --gnomad-exome-path ${gnomad_exome} ${true='--include-batch-freq' false='' batch_freq} --finngen-path ${finngen_annotation} \
        --compare-style ${compare_style} ${check_for_ld} --gwascatalog-pval ${gw_pval} \
        --gwascatalog-width-kb ${gw_width} --ld-treshold ${ld_treshold}
    }
    #output
    output {
        File gws_out = "fetch_out.csv"
        File annotate_out = "annotate_out.csv"
        File raport_out = "raport_out.csv"
        File ld_out = "ld_raport_out.csv"
    }
    runtime {
        docker: "${docker}"
        cpu: 3
        memory: "16 GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }
}

workflow autoreporting{
    call report {}

}