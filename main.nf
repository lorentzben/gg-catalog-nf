#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    INPUT AND VARIABLES
========================================================================================
*/

// Input

params.input = null
pacbio = false
iontorrent = false
single_end = params.single_end
if (params.pacbio || params.iontorrent) {
    single_end = true
}

params.multiple_sequencing_runs = false
params.extension = "/*_R{1,2}_001.fastq.gz"

// Set non-params Variables

String[] fasta_extensions = [".fasta", ".fna", ".fa"] // this is the alternative ASV fasta input
is_fasta_input = WorkflowGGCat.checkIfFileHasExtension( params.input.toString().toLowerCase(), fasta_extensions )

params.metadata = "${projectDir}/metadata.tsv"
params.outdir = "results"

log.info """\
         V I S U A L I Z E   P I P E L I N E    
         ===================================
         input    : ${params.input}
         single_end : ${params.single_end}
         pacbio : ${params.pacbio}
         iontorrent : ${params.iontorrent}
         multiple seq runs: ${params.multiple_sequencing_runs}
         extension : ${params.extension}
         metadata : ${params.metadata}
         outdir   : ${params.outdir}
         profile : ${workflow.profile}
         """
         .stripIndent()

/*
  Import processes from external files
  It is common to name processes with UPPERCASE strings, to make
  the program more readable (this is of course not mandatory)
*/
include { FILTLONG } from '../modules/nf-core/filtlong/main'
include { MINIMAP2_ALIGN } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../modules/nf-core/minimap2/index/main'


input_ch = Channel.fromPath(params.input, checkIfExists: true)


workflow {

    //
    // Create a channel for input read files
    //
    PARSE_INPUT ( params.input, is_fasta_input, single_end, params.multiple_sequencing_runs, params.extension )

    ch_reads = PARSE_INPUT.out.reads
    ch_fasta = PARSE_INPUT.out.fasta
    print(ch_reads)
    print(ch_fasta)
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
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/filtlong:2.0' : 'lorentzb/filtlong:2.0' }"

    input: 
    path reads

    //TODO if we get a locale issue, can we call the command through python or R?

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

    filtlong --min_length 2000 --keep_percent 99 ../zhang/reads/duodenum/SRR19683891.fastq | gzip > SRR19683891.fastq.gz
    '''
}

process MINIMAP{
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/minimap2:1.0' : 'lorentzb/minimap2:1.0' }"

    input:
    path filtered

    output:
    path uncontam

    script:
    '''
    #!/usr/bin/env bash

    
    '''
}

process READTABLE{

}

