
#' Filter fastq file by quality
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_fastp_qf <-
  function(
    head_label,
    in_dir = "./fastq",
    out_dir = "./fastq_qf"
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

mkdir {out_dir}
{in_dir}/{head_label}.{fq_ext}.gz {le}
fastp -q 20 -u 15 {le}
  -i {in_dir}/{head_label}.{fq_ext}.gz {le}
  -o {out_dir}/{head_label}.{fq_ext}.gz

    ")
  }

#' Filter fastq file by quality
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_fastp_qf_pe <-
  function(
    head_label,
    in_dir = "./fastq",
    out_dir = "./fastq_qf"
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

mkdir {out_dir}
./fastp -q 20 -u 15 {le}
  -i {in_dir}/{head_label}_1.{fq_ext}.gz {le}
  -I {in_dir}/{head_label}_2.{fq_ext}.gz {le}
  -o {out_dir}/{head_label}_1.{fq_ext}.gz {le}
  -O {out_dir}/{head_label}_2.{fq_ext}.gz

    ")
  }
