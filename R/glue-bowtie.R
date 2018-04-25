
#' Generate genome index for bowtie
#' @importFrom glue glue
#' @param fasta_path fasta file containing reference sequences to be aligned for
#' @param index_name Name of the index
#' @export
glue_bowtie_genome_generate <-
  function(fasta_path, index_name){
    glue("
bowtie-build -f {fasta_path} {index_name}
         ")
  }

#' Map with bowtie and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param index_name .tar file path
#' @param core_num number of threads to use
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_bowtie_bamsort <-
  function(
    head_label,
    index_name,
    in_dir = "./fastq",
    out_dir = "./mapped_by_bowtie",
    core_num = 2
  ){
    lineend <- "\\"
    glue(
      "

mkdir {out_dir}
bowtie {lineend}
  -p {core_num} {lineend}
  -S {lineend}
  -a --best --strata {line_end}
  {index_name} {lineend}
  {in_dir}/{head_label}.fastq {lineend}
| {lineend}
samtools view -@ {core_num} -bS {lineend}
| {lineend}
samtools sort -@ {core_num} > {out_dir}/{head_label}.sort.bam
samtools index {out_dir}/{head_label}_bowtie.sort.bam

      "
    )
  }
