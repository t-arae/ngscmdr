
#' Clipping adaptor from single end reads
#' @importFrom glue glue
#' @param label .tar file path
#' @param adaptor .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_se_cutadapt <-
  function(
    label,
    adaptor,
    in_dir = "./fastq",
    out_dir = "./fastq"
  ){
    lineend <- "\\"
    glue("

cutadapt {lineend}
  -a {adaptor} {lineend}
  -o {out_dir}/{label}_cl.fastq.gz {lineend}
  {in_dir}/{label}.fastq.gz {lineend}

         ")
  }

#' Clipping adaptor from pair end reads
#' @importFrom glue glue
#' @param label .tar file path
#' @param f_adaptor .tar file path
#' @param r_adaptor .tar file path
#' @param overlap .tar file path
#' @param minimum_length .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_pe_cutadapt <-
  function(
    label,
    f_adaptor,
    r_adaptor,
    overlap = 3,
    minimum_length = 20,
    in_dir = "./fastq",
    out_dir = "./fastq"
  ){
    lineend <- "\\"
    glue("

cutadapt {lineend}
  --pair-filter=any {lineend}
  -O {overlap} {lineend}
  -m {minimum_length} {lineend}
  -a {f_adaptor} {lineend}
  -A {r_adaptor} {lineend}
  -o {out_dir}/{label}_cl_1.fastq.gz {lineend}
  -p {out_dir}/{label}_cl_2.fastq.gz {lineend}
  {in_dir}/{label}_1.fastq.gz {in_dir}/{label}_2.fastq.gz

         ")
  }
