
#' Generate genome index for HISAT2
#' @importFrom glue glue
#' @param chr_fasta_path .tar file path
#' @param index_name .tar file path
#' @export
glue_hisat_genome_generate <-
  function(
    chr_fasta_path = "./allChr.fas",
    index_name = "TAIR10"
  ){
    lineend <- "\\"
    glue(
      "

hisat2-build {lineend}
  -p 2 {lineend}
  --seed 1 {lineend}
  {chr_fasta_path} {lineend}
  {index_name}

    "
    )
  }

#' Map with STAR and output a sorted bam file
#' @importFrom glue glue
#' @param read1 .tar file path
#' @param read2 .tar file path
#' @param core_num .tar file path
#' @param index_name .tar file path
#' @param out .tar file path
#' @export
glue_hisat_bamsort <-
  function(read1, read2,
           core_num = 2,
           index_name = "TAIR10",
           out = "hoge"
  ){
    lineend <- "\\"
    glue(
      "

hisat2 {lineend}
  -p {core_num} {lineend}
  -x {index_name} {lineend}
  -1 {read1} {lineend}
  -2 {read2} {lineend}
| {lineend}
samtools view -@ {core_num} -bS {lineend}
| {lineend}
samtools sort -@ {core_num} > {out}.sort.bam
samtools index {out}.sort.bam

    "
    )
  }
