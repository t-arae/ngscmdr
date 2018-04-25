
#' Generate genome index for bowtie2
#' @importFrom glue glue
#' @param fasta_path .tar file path
#' @param index_name .tar file path
#' @export
glue_bowtie2_genome_generate <-
  function(fasta_path, index_name){
    glue("
bowtie2-build -f {fasta_path} {index_name}
         ")
  }

#' Map with bowtie2 and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param index_name .tar file path
#' @param core_num .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_bowtie2_bamsort <-
  function(
    head_label,
    index_name,
    core_num = 2,
    in_dir = "./fastq",
    out_dir = "./mapped_by_bowtie2"
  ){
    lineend <- "\\"
    glue(
      "
mkdir {out_dir}
bowtie2 {lineend}
  -p {core_num} {lineend}
  -a --best --strata {lineend}
  --trim5 8 {lineend}
  --trim3 40 {lineend}
  -x {index_name} {lineend}
  -U {in_dir}/{head_label}.fastq {lineend}
| {lineend}
samtools view -@ {core_num} -bS {lineend}
| {lineend}
samtools sort -@ {core_num} > {out_dir}/{head_label}_bowtie2.sort.bam
samtools index {out_dir}/{head_label}_bowtie2.sort.bam && {line_end}

      "
    )
  }
