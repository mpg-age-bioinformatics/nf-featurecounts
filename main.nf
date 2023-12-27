#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.run_type}" == "r2d2" ]] || [[ "${params.run_type}" == "raven" ]] ; 

      then

        cd ${params.image_folder}

        if [[ ! -f subread-2.0.3.sif ]] ;
          then
            singularity pull subread-2.0.3.sif docker://index.docker.io/mpgagebioinformatics/subread:2.0.3
        fi

        if [[ ! -f rnaseq.python-3.8-1.sif ]] ;
          then
            singularity pull rnaseq.python-3.8-1.sif docker://index.docker.io/mpgagebioinformatics/rnaseq.python:3.8-1
        fi

    fi


    if [[ "${params.run_type}" == "local" ]] ; 

      then

        docker pull mpgagebioinformatics/subread:2.0.3
        docker pull mpgagebioinformatics/rnaseq.python:3.8-1

    fi

    """

}

process exomeBED_to_exomeGTF {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    path(bed)
  
  when:
    ( ! file("${params.project_folder}/tmp/target_exons.gtf").exists() )

  script:
    """
    mkdir -p /workdir/tmp/
    awk '{print \$1 "\\tillumina\\texon\\t" \$2 "\\t" \$3 "\\t.\\t.\\t.\\tgene_id \\\"" \$1 "_" \$2 "_" \$3 "\\\";" }' ${bed} > /workdir/tmp/target_exons.gtf
    """
}

process featurecounts_general {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val mapping_output
    val gtf
    val bam
    val strand_file
    tuple val(pair_id), path(fastq)

  output:
    val pair_id

  when:
    ( ! file("${params.project_folder}/featureCounts_output/${pair_id}_gene.featureCounts.txt").exists() ) 
  
  script:
  def single = fastq instanceof Path

  if ( single ) {
    paired=""
  } else {
    paired="-p"
  }
  """
    mkdir -p /workdir/featureCounts_output

    if [[ -e "${strand_file}" ]] ; then strand=\$(cat ${strand_file}) ; else strand="0" ; fi
    echo \${strand}

    if [[ "${mapping_output}" == "kallisto_output" ]]
    then
      cd /workdir/${mapping_output}/${pair_id}
      if [[ "${bam}" == "None" ]] ; then bam=${pair_id}/\$(ls *.bam | head -n 1 ) ; else bam=${pair_id}/${bam} ; fi
    
    elif [[ "${mapping_output}" == "bwa_output" ]]
    then
      if [[ "${bam}" == "None" ]] ; then bam=${pair_id}.sorted.bam ; else bam=${bam} ; fi
    fi

    featureCounts -a ${gtf} -T ${task.cpus} -g gene_id -o /workdir/featureCounts_output/${pair_id}_gene.featureCounts.txt ${paired} -s \${strand} /workdir/${mapping_output}/\${bam}

  """
}

process featurecounts_rnaSeq_specific {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val mapping_output
    val gtf
    val bam
    val strand_file
    tuple val(pair_id), path(fastq)

  output:
    val pair_id

  when:
    ( ! file("${params.project_folder}/featureCounts_output/${pair_id}_biotype_counts_mqc.txt").exists() ) 
  
  script:
  def single = fastq instanceof Path

  if ( single ) {
    paired=""
  } else {
    paired="-p"
  }
  """
    strand=\$(cat ${strand_file})
    mkdir -p /workdir/featureCounts_output
    cd /workdir/${mapping_output}/${pair_id}
    if [[ "${bam}" == "None" ]] ; then bam=\$(ls *.bam | head -n 1 ) ; else bam=${bam} ; fi
    featureCounts -a ${gtf} -T ${task.cpus} -g gene_biotype -o /workdir/featureCounts_output/${pair_id}_biotype.featureCounts.txt ${paired} -s \${strand} \${bam}
    cp /workdir/featureCounts_output/biotypes_header.txt /workdir/featureCounts_output/${pair_id}_biotype_counts_mqc.txt
    cut -f 1,7 /workdir/featureCounts_output/${pair_id}_biotype.featureCounts.txt | tail -n +3 | grep -v '^\\s' >> /workdir/featureCounts_output/${pair_id}_biotype_counts_mqc.txt
  """
}

process featurecounts_headers {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val bam
    val pair_id

  script:
    """
    #!/usr/local/bin/python
    import os
    import csv
    import pandas as pd
    import numpy as np
    files_path="/workdir/featureCounts_output/"
    files=os.listdir(files_path)
    files=[s for s in files if "summary" not in s and "mqc" not in s]
    gene_files=[ s for s in files if s [-len("_gene.featureCounts.txt"):] == "_gene.featureCounts.txt" ]
    biotype_files=[ s for s in files if s[-len("_biotype.featureCounts.txt"):] == "_biotype.featureCounts.txt" ]
    files=gene_files+biotype_files
    for f in files:
        file=pd.read_csv(files_path+f, sep="\\t")
        header=file.columns[0]
        file_=pd.read_csv(files_path+f, sep="\\t", comment="#")
        
        if "_gene.featureCounts.txt" in f:
            sample_name=f.split("_gene",1)[0]
        elif "_biotype.featureCounts.txt" in f:
            sample_name=f.split("_biotype",1)[0]
        
        file_=file_.rename(columns={"${bam}":sample_name})
        
        text_file = open(files_path+f, 'w')
        text_file.write(header+"\\n")
        text_file.close()
        
        file_.to_csv(files_path+f, mode='a', header=True, sep="\\t", index=None)
        
        file_sum=pd.read_csv(files_path+f+".summary", sep="\\t")
        file_sum=file_sum.rename(columns={"${bam}":f})
        file_sum.to_csv(files_path+f+".summary", index=None, sep="\\t")
    """
}


workflow images {
  main:
    get_images()
}

workflow exomeGTF {
  main:
  if ( 'exomebed' in params.keySet() ) {
    bed=params["exomebed"]
    if ( bed != "none" ) {
      exomeBED_to_exomeGTF( bed )
    }
  }
}

workflow {
  if ( 'mapping_output' in params.keySet() ) {
    // mapping_output=${params.mapping_output}
    mapping_output=params["mapping_output"]
  } else {
    mapping_output="kallisto_output"
  }

  if ( 'gtf' in params.keySet() ) {
    // gtf=${params.gtf}
    gtf=params["gtf"]
  } else {
    gtf="/workdir/kallisto_index/${params.organism}.${params.release}.no.rRNA.gtf"
  }

  if ( 'bam' in params.keySet() ) {
    // bam=${params.bam}
    bam=params["bam"]
  } else {
    bam="pseudoalignments.bam"
  }

  if ( 'strand_file' in params.keySet() ) {
    // strand_file=${params.strand_file}
    strand_file=params["strand_file"]
  } else {
    strand_file="/workdir/featureCounts_output/.strandness.txt"
  }

  read_files=Channel.fromFilePairs( "${params.featurecounts_raw_data}/*.READ_{1,2}.fastq.gz", size: -1 )

  folder=file("${params.project_folder}/featureCounts_output/")
  if( ! folder.isDirectory() ) {
    folder.mkdirs()
  }

  featurecounts_general(mapping_output, gtf, bam, strand_file, read_files)

  if( ! ('model' in params.keySet()) ) {
  
    file=file("${params.project_folder}/featureCounts_output/biotypes_header.txt")
    headers = """# id: 'biotype-counts'\n\
  # section_name: 'Biotype Counts'\n\
  # description: "shows reads overlapping genomic features of different biotypes,\n\
  #     counted by <a href='http://bioinf.wehi.edu.au/featureCounts'>featureCounts</a>."\n\
  # plot_type: 'bargraph'\n\
  # anchor: 'featurecounts_biotype'\n\
  # pconfig:\n\
  #     id: "featureCounts_biotype_plot"\n\
  #     title: "featureCounts: Biotypes"\n\
  #     xlab: "# Reads"\n\
  #     cpswitch_counts_label: "Number of Reads\n"""
    file.text = headers
  
    featurecounts_rnaSeq_specific(mapping_output, gtf, bam, strand_file, read_files)

  }

  if ( bam == "pseudoalignments.bam" ) {
    featurecounts_headers(bam, featurecounts_rnaSeq_specific.out.collect() )
  }

}

