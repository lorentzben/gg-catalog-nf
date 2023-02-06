#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    INPUT AND VARIABLES
========================================================================================
*/

// Input

params.input = null
params.contam = null
params.pacbio = false
params.iontorrent = false
params.single_end = false
single_end = params.single_end
if (params.pacbio || params.iontorrent) {
    single_end = true
}

params.extension = "/*_R{1,2}_001.fastq.gz"
params.metadata = null

// Set non-params Variables

String[] fasta_extensions = [".fasta", ".fna", ".fa"] // this is the alternative ASV fasta input
is_fasta_input = WorkflowGGCat.checkIfFileHasExtension( params.input.toString().toLowerCase(), fasta_extensions )


params.outdir = "results"

log.info """\
         V I S U A L I Z E   P I P E L I N E    
         ===================================
         input    : ${params.input}
         contam    : ${params.contam}
         single_end : ${params.single_end}
         pacbio : ${params.pacbio}
         iontorrent : ${params.iontorrent}
         extension : ${params.extension}
         metadata: ${params.metadata}
         outdir   : ${params.outdir}
         profile : ${workflow.profile}
         """
         .stripIndent()

/*
  Import processes from external files
  It is common to name processes with UPPERCASE strings, to make
  the program more readable (this is of course not mandatory)
*/
//if don't work use projectDir

include { MINIMAP2_ALIGN } from "${moduleDir}/modules/nf-core/minimap2/align/main"
include { MINIMAP2_INDEX } from "${moduleDir}/modules/nf-core/minimap2/index/main"
include { PARSE_INPUT } from "${projectDir}/subworkflows/local/parse_input"
include { CONTAM_INPUT } from "${projectDir}/subworkflows/local/contam_input"


input_ch = Channel.fromPath(params.input, checkIfExists: true)
contam_ch = Channel.fromPath(params.contam, checkIfExists: true)


workflow{

    //
    // Create a channel for input read files
    //
    PARSE_INPUT ( params.input, is_fasta_input, single_end, params.extension )

    ch_reads = PARSE_INPUT.out.reads
    ch_fasta = PARSE_INPUT.out.fasta

    id_ch = ch_reads.map{it.first()}
    path_ch = ch_reads.map{it.last()}

    FILTLONG(id_ch,path_ch)

    CONTAM_INPUT(params.contam, false, true, "*.fna.gz")

    ch_contam_reads = CONTAM_INPUT.out.reads
    ch_contam_fasta = CONTAM_INPUT.out.fasta

    id_contam_ch = ch_contam_reads.map{it.first()}
    path_contam_ch = ch_contam_reads.map{it.last()}

    MINIMAP2_INDEX(ch_contam_reads)

    ref_1 = MINIMAP2_ALIGN(tuple(id_ch,FILTLONG.out.filtered),MINIMAP2_INDEX.out.index, true, false, true)

    view(ref_1.bam)
    
    
}

process FILTLONG{
    label 'process_high'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/filtlong:2.0' : 'lorentzb/filtlong:2.0' }"

    input:
    val meta
    path reads

    output:
    path "*.fastq.gz", emit: filtered

    script:
    
    """
    #!/usr/bin/env bash

    ${reads} 
    ${meta} 

    filtlong --min_length 2000 --keep_percent 99 ${reads} | gzip > ${meta.id}.fastq.gz

    """

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

