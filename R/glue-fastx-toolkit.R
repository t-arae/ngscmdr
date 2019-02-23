
#' Filter fastq file by quality and Trim adaptor sequence
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_fastq_qf_cl <-
  function(
    head_label,
    in_dir = "./fastq",
    out_dir = "./fastq_qf_cl"
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

mkdir {out_dir}
gunzip -c {in_dir}/{head_label}.{fq_ext}.gz {le}
| {le}
fastq_quality_filter -q 20 -p 85 {le}
| {le}
fastx_clipper -a AGATCGGAAGAGCACACGTCT -z -o {out_dir}/{head_label}.{fq_ext}.gz

    ")
  }

#' Filter fastq file by quality
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_fastq_qf <-
  function(
    head_label,
    in_dir = "./fastq",
    out_dir = "./fastq_qf"
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

mkdir {out_dir}
gunzip -c {in_dir}/{head_label}.{fq_ext}.gz {le}
| {le}
fastq_quality_filter -q 20 -p 85 -z -o {out_dir}/{head_label}.{fq_ext}.gz

    ")
  }
