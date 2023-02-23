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
include { FILTLONG } from "${projectDir}/modules/nf-core/filtlong/main"
include { SAMTOOLS_FASTQ } from "${projectDir}/modules/nf-core/samtools/fastq/main"
include { SAMTOOLS_FASTA } from "${projectDir}/modules/nf-core/samtools/fasta/main"
include { SEQKIT_STATS; SEQKIT_STATS as SEQKIT_STATS_FILT; SEQKIT_STATS as SEQKIT_STATS_UNMAP } from "${projectDir}/modules/nf-core/seqkit/stats/main"
include { CSVTK_CONCAT; CSVTK_CONCAT as CSVTK_CONCAT_FILT; CSVTK_CONCAT as CSVTK_CONCAT_UNMAP } from "${projectDir}/modules/nf-core/csvtk/concat/main"


input_ch = Channel.fromPath(params.input, checkIfExists: true)
contam_ch = Channel.fromPath(params.contam, checkIfExists: true)


workflow{

    //
    // Create a channel for input read files
    //
    PARSE_INPUT ( params.input, is_fasta_input, single_end, params.extension )

    ch_reads = PARSE_INPUT.out.reads
    ch_fasta = PARSE_INPUT.out.fasta

    ch_reads_mod = ch_reads.map{
        it ->  [ it[0], [], it[1].flatten() ]
    }

    id_ch = ch_reads.map{it.first()}
    path_ch = ch_reads.map{it.last()}
   
    FILTLONG(ch_reads_mod)
    
    CONTAM_INPUT(params.contam, false, true, "*.fna.gz")

    ch_contam_reads = CONTAM_INPUT.out.reads
    ch_contam_fasta = CONTAM_INPUT.out.fasta

    id_contam_ch = ch_contam_reads.map{it.first()}
    path_contam_ch = ch_contam_reads.map{it.last()}

    MINIMAP2_INDEX(ch_contam_reads)

    
    meta_ch = MINIMAP2_INDEX.out.index.map{it.first()}
    index_path_ch = MINIMAP2_INDEX.out.index.map{it.last()}

    ch_contam_reads
        .map{it.last()}
        .set{ contam_path_ch }

    FILTLONG.out.reads
        .map{ meta, reads -> [ meta, reads ] }
        .set{ ch_filtered }


    MINIMAP2_ALIGN(ch_filtered ,contam_path_ch.first(), true, false, true)

    SAMTOOLS_FASTQ(MINIMAP2_ALIGN.out.bam , false)

    SAMTOOLS_FASTA(MINIMAP2_ALIGN.out.bam , false)

    // raw reads process


    raw_reads = SEQKIT_STATS(ch_reads)

    raw_reads.stats
        .map{ file(it.last()) }
        .collect()
        .set{ch_raw_table_loc}

    Channel
        .of([id:"raw", single_end:true])
        .set{ch_meta_raw}

    
    
    ch_raw_table_loc
        .map{
            it -> [[id:"raw", single_end:true], it]
        }
        .set{ch_raw_table}

    
    CSVTK_CONCAT(ch_raw_table,'tsv','tsv')

    // filtlong filtered process

    filtlong_reads = SEQKIT_STATS_FILT(ch_filtered)

    filtlong_reads.stats
        .map{ file(it.last()) }
        .collect()
        .set{ch_filtlong_table_loc}

    ch_filtlong_table_loc
        .map{
            it -> [[id: "filtlong"], it]
        }
        .set{ch_filtlong_table}
    
    CSVTK_CONCAT_FILT(ch_filtlong_table,'tsv','tsv')

    // minimap2 reads

    unmapped_reads = SAMTOOLS_FASTQ.out.fastq.concat(SAMTOOLS_FASTQ.out.interleaved).concat(SAMTOOLS_FASTQ.out.singleton).concat(SAMTOOLS_FASTQ.out.other)

    SAMTOOLS_FASTQ.out.fastq
        .map{ 
            it -> [[id : "fastq-"it.first().id], it.last()]
            }
        .set{fastq_reads_ch}

    fastq_reads_ch.view()
    //unmapped_reads = unmapped_reads.concat(SAMTOOLS_FASTQ.out.singleton)

    minimap_reads = SEQKIT_STATS_UNMAP(unmapped_reads)

    minimap_reads.stats
        .map{ file(it.last()) }
        .collect()
        .set{ch_minimap_table_loc}

    ch_minimap_table_loc
        .map{
            it -> [[id: "filtlong"], it]
        }
        .set{ch_minimap_table}
    
    //CSVTK_CONCAT_UNMAP(ch_minimap_table,'tsv','tsv')


    
    
    
}

/*
process FILTLONG{
    label 'process_high'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/filtlong:2.0' : 'lorentzb/filtlong:2.0' }"

    input:
    val meta
    path reads

    output:
    path "*.fastq.gz", emit: filtered
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: dual
    

    script:
    
    """
    #!/usr/bin/env bash

    ${reads} 
    ${meta} 

    filtlong --min_length 2000 --keep_percent 99 ${reads} | gzip > ${meta.id}.fastq.gz

    """

}
*/

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

