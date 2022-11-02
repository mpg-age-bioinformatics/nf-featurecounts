#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process featurecounts {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    tuple val(pair_id), path(fastq)

  when:
    ( ! file("/workdir/featureCounts_output/${pair_id}_biotype_counts_mqc.txt").exists() ) 
  
  script:
  def single = fastq instanceof Path

  if ( 'mapping_output' in params.keySet() ) {
    mapping_output=${params.mapping_output}
  } else {
    mapping_output="kallisto_output"
  }

  if ( 'gtf' in params.keySet() ) {
    gtf=${params.gtf}
  } else {
    gtf="/workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf"
  }

  if ( 'bam' in params.keySet() ) {
    bam=${params.bam}
  } else {
    bam="pseudoalignments.bam"
  }

  if ( 'strand_file' in params.keySet() ) {
    strand_file=${params.strand_file}
  } else {
    strand_file="/workdir/featureCounts_output/.strandness.txt"
  }

  if ( single ) {
    paired=""
  } else {
    paired="-p"
  }

  """
    strand=\$(cat ${strand_file})
    mkdir -p /workdir/featureCounts_output
    cd /workdir/${mapping_output}/${pair_id}
    featureCounts -a ${gtf} -T ${task.cpus} -g gene_id -o /workdir/featureCounts_output/${pair_id}_gene.featureCounts.txt ${paired} -s \${strand} ${bam}
  """

}

workflow {
  read_files=Channel.fromFilePairs( "${params.kallisto_raw_data}*.READ_{1,2}.fastq.gz", size: -1 )
  featurecounts(read_files)
}




