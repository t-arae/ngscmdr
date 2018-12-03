
#' Create fastqc reports
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_bamstat <-
  function(
    head_label,
    in_dir,
    out_dir = in_dir
  ){
    le <- "\\"
    glue("

mkdir {out_dir}/samtoools_stat
samtools stats {in_dir}/{head_label}.sort.bam > {le}
  {out_dir}/samtoools_stat/{head_label}_samtools_stats.txt
samtools flagstat {in_dir}/{head_label}.sort.bam > {le}
  {out_dir}/samtoools_stat/{head_label}_samtools_flagstat.txt
samtools idxstats {in_dir}/{head_label}.sort.bam > {le}
  {out_dir}/samtoools_stat/{head_label}_samtools_idxstats.txt

       ")
  }

