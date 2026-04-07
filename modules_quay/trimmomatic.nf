process trimmomatic_pe {

  container "quay.io/biocontainers/trimmomatic:0.39--1"
  publishDir "analysis/trimmomatic", pattern: "*.fastq.gz"
  publishDir "report/logs/", pattern: "*.log"
  cpus params.trimmomatic_cpus

  input:
    tuple path(read1), path(read2), val(name)
    val trimmers
    val meta

  output:
    tuple path("${name}_trimmed_1.fastq.gz"), path("${name}_trimmed_2.fastq.gz"), val(name), emit: data
    path "${name}_trimmomatic.log", emit: log
    val meta, emit: meta

  script:
    _trimmers = trimmers.join(" ")
    """
    run.sh PE -threads ${params.trimmomatic_cpus} \
      ${read1} ${read2} \
      ${name}_trimmed_1.fastq.gz ${name}_orphans_1.fastq.gz \
      ${name}_trimmed_2.fastq.gz ${name}_orphans_2.fastq.gz \
      ${_trimmers} &> ${name}_trimmomatic.log
    """
}

def map_input_pe(Ch) {
  Ch.map(row -> { 
    ["${baseDir}/" + row.RF, "${baseDir}/" + row.RS, row.Id] 
  })
}

def bibtex() {
  """
  @article{bolger2014trimmomatic,
    title={Trimmomatic: a flexible trimmer for Illumina sequence data},
    author={Bolger, Anthony M and Lohse, Marc and Usadel, Bjoern},
    journal={Bioinformatics},
    volume={30},
    number={15},
    pages={2114--2120},
    year={2014},
    publisher={Oxford University Press}
  }
  """
}

def citekey() {
  "@bolger2014trimmomatic"
}