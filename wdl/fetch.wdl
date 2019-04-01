workflow fetch_gws{
    call fetch {}

    output {
        Array[File] output_files=fetch.output_files 
    }
}

task fetch{
    #variables

    #command

    #docker

    #output

}