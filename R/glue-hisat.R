
#' Generate genome index for HISAT2
#' @importFrom glue glue
#' @param chr_fasta_path .tar file path
#' @param index_name .tar file path
#' @param index_dir .tar file path
#' @param core_num .tar file path
#' @export
glue_hisat_genome_generate <-
  function(
    chr_fasta_path = "./allChr.fas",
    index_name = "TAIR10",
    index_dir = "./genomeidx_hisat",
    core_num = 2
  ){
    lineend <- "\\"
    glue(
      "

mkdir {index_dir}
hisat2-build {lineend}
  -p {core_num} {lineend}
  --seed 1 {lineend}
  {chr_fasta_path} {lineend}
  {index_dir}/{index_name}

    "
    )
  }

#' Map single-end reads with hisat and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param core_num .tar file path
#' @param index_name .tar file path
#' @param index_dir .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_se_hisat_bamsort <-
  function(head_label,
           core_num = 2,
           index_name = "TAIR10",
           index_dir = "./genomeidx_hisat",
           in_dir = "./fastq",
           out_dir = "./mapped_by_hisat"
  ){
    lineend <- "\\"
    glue(
      "

mkdir {out_dir}
hisat2 {lineend}
  --dta-cufflinks {lineend}
  -p {core_num} {lineend}
  -x {index_dir}/{index_name} {lineend}
  -U {in_dir}/{head_label}.fastq.gz {lineend}
| {lineend}
samtools view -@ {core_num} -bS {lineend}
| {lineend}
samtools sort -@ {core_num} > {out_dir}/{head_label}.sort.bam
samtools index {out_dir}/{head_label}.sort.bam

    "
    )
  }


#' Map paired-end reads with hisat and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param core_num .tar file path
#' @param index_name .tar file path
#' @param index_dir .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_pe_hisat_bamsort <-
  function(head_label,
           core_num = 2,
           index_name = "TAIR10",
           index_dir = "./genomeidx_hisat",
           in_dir = "./fastq",
           out_dir = "./mapped_by_hisat"
  ){
    lineend <- "\\"
    glue(
      "

mkdir {out_dir}
hisat2 {lineend}
  --dta-cufflinks {lineend}
  -p {core_num} {lineend}
  -x {index_dir}/{index_name} {lineend}
  -1 {in_dir}/{head_label}_1.fastq.gz {lineend}
  -2 {in_dir}/{head_label}_2.fastq.gz {lineend}
| {lineend}
samtools view -@ {core_num} -bS {lineend}
| {lineend}
samtools sort -@ {core_num} > {out_dir}/{head_label}.sort.bam
samtools index {out_dir}/{head_label}.sort.bam

    "
    )
  }
