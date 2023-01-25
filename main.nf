#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "${projectDir}/reads"
params.metadata = "${projectDir}/metadata.tsv"
params.outdir = "results"


log.info """\
         V I S U A L I Z E   P I P E L I N E    
         ===================================
         input    : ${params.input }
         metadata : ${params.metadata}
         outdir   : ${params.outdir}
         profile : ${workflow.profile}
         """
         .stripIndent()

input_ch = Channel.fromPath(params.input, checkIfExists: true)
metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)


workflow {
    ord_ioi = ORDERIOI(ioi_ch, metadata_ch, ord_ioi_ch)
    
}

process RAREFACTIONPLOT{
    publishDir "${params.outdir}/rarefaction", pattern: "*.png", mode: "copy"
    publishDir "${params.outdir}/rarefaction", pattern: "*.csv", mode: "copy"
    publishDir "${params.outdir}", pattern: "*.html", mode: "copy"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'docker://lorentzb/tidyverse:4.2.0' : 'lorentzb/tidyverse:4.2.0' }"

    input:
    path 'results'
    path report
 

    output:
    path("*.png"), emit: rare_images
    path("*.csv"), emit: rare_tabs
    path("*.html"), emit: rare_report

    script:

    '''
    #!/usr/bin/env bash

    Rscript -e "rmarkdown::render('rarefaction_report.Rmd', output_file='$PWD/rarefaction_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    '''

}

