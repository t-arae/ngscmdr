
#' Clipping adaptor from single end reads
#' @importFrom glue glue
#' @inheritParams param_general
#' @param adaptor .tar file path
#' @export
glue_se_cutadapt <-
  function(
    head_label,
    adaptor,
    in_dir = "./fastq",
    out_dir = "./fastq"
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

cutadapt {le}
  -a {adaptor} {le}
  -o {out_dir}/{head_label}_cl.{fq_ext}.gz {le}
  {in_dir}/{head_label}.{fq_ext}.gz {le}

         ")
  }

#' Clipping adaptor from pair end reads
#' @importFrom glue glue
#' @inheritParams param_general
#' @param f_adaptor .tar file path
#' @param r_adaptor .tar file path
#' @param overlap .tar file path
#' @param minimum_length .tar file path
#' @export
glue_pe_cutadapt <-
  function(
    head_label,
    f_adaptor,
    r_adaptor,
    overlap = 3,
    minimum_length = 20,
    in_dir = "./fastq",
    out_dir = "./fastq"
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

cutadapt {le}
  --pair-filter=any {le}
  -O {overlap} {le}
  -m {minimum_length} {le}
  -a {f_adaptor} {le}
  -A {r_adaptor} {le}
  -o {out_dir}/{head_label}_cl_1.{fq_ext}.gz {le}
  -p {out_dir}/{head_label}_cl_2.{fq_ext}.gz {le}
  {in_dir}/{head_label}_1.{fq_ext}.gz {in_dir}/{head_label}_2.{fq_ext}.gz

         ")
  }
