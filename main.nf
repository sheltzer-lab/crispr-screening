#!/usr/bin/env nextflow

// params.SAMPLES
//   Format: CSV file containing Case/Control, Time Point, Adapter
//   Example (remove leading indent):
//     sample,time,sequence
//     Cont1,initial,ATCACG
//     Cont1,final,TTAGGC
//     Ts13,initial,CGATGT
//     Ts13,final,TGACCA
//
// Note: {sample}.{time} should match the contents of matrix{sample}.txt
// so matrixCont1.txt contains a line for Cont1.initial and Cont1.final
params.SAMPLES="screening-example-samples.csv"

// params.MATRICES
//   Format: command separated list of tab-separated matrices
//   Example (remove leading indent):
//     Samples <tab>  baseline <tab> Cont1
//     Cont1.initial  1              0
//     Cont1.final    1              1
//
// Note: {sample}.{time} should match the contents of params.SAMPLES
params.MATRICES="matrixCont1.txt,matrixTs13.txt"

// params.LIBRARY
//   Format: line-based, each line is the guide name and sequence
//   Example (remove leading indent):
//     gene name <tab> unique identifier <tab> sequence
//     SLC25A24        108157552               AAAAAAAATCCGGACAATGG
//     FASTKD3         7867646                 AAAAAAAGGATGGTGATCAA
//     BCAS2           114575706               AAAAAAATGACATTACTGCA
//     GPR18           99255304                AAAAAAATGTCAGTCGAGTG
//     ZNF470          56577092                AAAAAACACAAGCAAGACCG
//     ...
params.LIBRARY="WGL.txt"

// params.SEQ_FQS
//   Format: raw FASTQ reads from sequencer
//   IMPORTANT: The params.PRIMER must consistently be found in this file. If not, then the reverse file is needed instead, but only one direction is read
params.SEQ_FQS="Kinase_Library_DLD1_1-11320_S1_R1_001.fastq"

// params.SEQ_FQS
//   Format: the DNA primer sequence being used
//   IMPORTANT: This must consistently be found in params.SEQ_FQS.
params.PRIMER="GTGGAAAGGACGAAACACC"

workflow {
  def ch_seqfqs = Channel.fromPath(params.SEQ_FQS, checkIfExists: true)
  def ch_samples = Channel.fromPath(params.SAMPLES, checkIfExists: true)
  def ch_library = Channel.fromPath(params.LIBRARY, checkIfExists: true)
  def ch_matrices = Channel.fromPath(params.MATRICES.tokenize(','), checkIfExists: true)
  
  extract_reads(params.PRIMER, ch_seqfqs, ch_samples, ch_library)
  run_mle(extract_reads.out, ch_matrices)
}

process run_mle {
  debug true
  publishDir "results/", mode: 'copy'
  conda 'bioconda::mageck==0.5.9'

  input:
    path count_files
    each path(matrix)

  output:
    path 'mle_*'

  script:
    def name = matrix.simpleName.replaceAll('matrix', '')
    def cfile = count_files.findAll( { it.simpleName.contains(name) } ).first()
    """
    mageck mle -k ${cfile} -d ${matrix} -n mle_${name}
    """
}

process extract_reads {
  conda 'pandas==1.4.2 python==3.10.4'

  input:
    val primer
    path reads, stageAs: 'reads.fastq.gz'
    path samples, stageAs: 'samples.csv'
    path library, stageAs: 'library.csv'

  output:
    path 'count-*-i.f.csv'

  script:
    """
    extract-reads.py ${primer} ${reads} ${samples} ${library}
    """
}
