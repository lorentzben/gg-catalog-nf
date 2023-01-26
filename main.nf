#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "${projectDir}/reads"
params.metadata = "${projectDir}/metadata.tsv"
params.outdir = "results"


log.info """\
         V I S U A L I Z E   P I P E L I N E    
         ===================================
         input    : ${params.input}
         metadata : ${params.metadata}
         outdir   : ${params.outdir}
         profile : ${workflow.profile}
         """
         .stripIndent()

input_ch = Channel.fromPath(params.input, checkIfExists: true)


workflow {
    TEST(input_ch)
    
}

process TEST{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/tidyverse:4.2.0' : 'lorentzb/tidyverse:4.2.0' }"

    input: 
    path input

    output:
    
    script:

    '''
    #!/usr/bin/env bash

    echo "Test"
    ls $input

    '''

}

process FILTLONG{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/filtlong:1.0' : 'lorentzb/filtlong:1.0' }"

    input: 
    path reads

    output:
    path ("*.fastq.gz"), emit: filtered

    script:
    '''
    #!/usr/bin/env bash

    READS="*.fastq"
    

    for read in $READS; do
        READNAME=$read[::-6]
        filtlong --min_length 2000 --keep_percent 99 $read | gzip > $READNAME.fastq.gz
    done

    
    '''
}

process MINIMAP{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/minimap2:1.0' : 'lorentzb/minimap2:1.0' }"

}

process READTABLE{

}

