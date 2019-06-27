task report {
    #variables
    String docker
    File summ_stat
    File summ_stat_tb=summ_stat+".tbi"
    File gnomad_genome
    File gnomad_genome_tb=gnomad_genome+".tbi"
    File gnomad_exome
    File gnomad_exome_tb=gnomad_exome+".tbi"
    String ld_panel
    File ld_panel_bed=ld_panel+".bed"
    File ld_panel_bim=ld_panel+".bim"
    File ld_panel_fam=ld_panel+".fam"
    File finngen_annotation
    File finngen_annotation_tb=finngen_annotation+".tbi"
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
    #Array[File] compare_summary_stats
    #Array[String] compare_endpoints
    String check_for_ld
    Float gw_pval
    Int gw_width
    Float ld_treshold
    String dollar = "$"

    #command
    command {
        mod_ld=$( echo ${ld_panel_bed} | sed 's/.bed//g' )

        main.py ${summ_stat} --sign-treshold ${s_tresh} --alt-sign-treshold ${s_tresh2} ${true='--group' false='' group} \
        --grouping-method ${gr_method} --locus-width-kb ${grouping_locus_width} --ld-panel-path ${dollar}mod_ld \
        --ld-r2 ${ld_r2} --plink-memory ${plink_mem} ${true='--overlap' false='' overlap} --gnomad-genome-path ${gnomad_genome} \
        --gnomad-exome-path ${gnomad_exome} ${true='--include-batch-freq' false='' batch_freq} --finngen-path ${finngen_annotation} \
        --compare-style ${compare_style} ${check_for_ld} --gwascatalog-pval ${gw_pval} \
        --gwascatalog-width-kb ${gw_width} --ld-treshold ${ld_treshold}
    }
    #output
    output {
        File gws_out = "fetch_out.csv"
        File annotate_out = "annotate_out.csv"
        File raport_out = "raport_output.csv"
        File ld_out = "ld_raport_out.csv"
    }
    runtime {
        docker: "${docker}"
        cpu: 3
        memory: "16 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2 
    }
}

workflow autoreporting{
    call report {}

}
