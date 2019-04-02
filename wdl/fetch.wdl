workflow fetch_gws{
    call fetch {}

    output {
        Array[File] output_files=fetch.output_files 
    }
}

task fetch{
    #variables
    File stat_file
    #command
    command {
        python3 /Scripts/gws_fetch.py stat_file
    }
    #docker

    #output

}