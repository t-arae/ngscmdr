
#' samtools stats/flagstat/idxstats
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


#' samtools sort
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_bam_sort <-
  function(
    head_label,
    in_dir,
    out_dir = in_dir
  ){
    le <- "\\"
    glue("

mkdir {out_dir}
samtools sort {in_dir}/{head_label}.bam > {out_dir}/{head_label}.sort.bam
samtools index {out_dir}/{head_label}.sort.bam

       ")
  }

#' samtools index
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_bam_index <-
  function(
    head_label,
    in_dir
  ){
    le <- "\\"
    glue("

samtools index {in_dir}/{head_label}.sort.bam

       ")
  }

#' bam2wig.py
#' @importFrom glue glue
#' @inheritParams param_general
#' @param chr_name_len_txt text file. chromatin name and length
#' @export
glue_bam2wig <-
  function(
    head_label,
    in_dir,
    chr_name_len_txt,
    out_dir = in_dir
  ){
    le <- "\\"
    glue("

bam2wig.py {le}
  -i {in_dir}/{head_label}.sort.bam {le}
  -s {chr_name_len_txt} {le}
  -o {out_dir}/{head_label}

       ")
  }

#' convert fastq.gz to rfq.xz
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_fqgz2rfqxz <-
  function(
    head_label,
    in_dir,
    out_dir = in_dir,
    core_num = 4
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

gzip -dc {in_dir}/{head_label}.{fq_ext}.gz {le}
| {le}
./repaq -c --stdin --stdout {le}
| {le}
xz -z  -9 -T{core_num} -c > {out_dir}/{head_label}.rfq.xz && rm {in_dir}/{head_label}.{fq_ext}.gz

       ")
  }

#' convert rqf.gz to  fastq.gz
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_rfqxz2fqgz <-
  function(
    head_label,
    in_dir,
    out_dir = in_dir,
    core_num = 4
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

xz -d -T{core_num} -c {in_dir}/{head_label}.rfq.xz {le}
| {le}
./repaq -d --stdin --stdout {le}
| {le}
gzip -c - > {out_dir}/{head_label}.{fq_ext}.gz

       ")
  }

#' Check the decompressed file size of gzipped file
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_check_gz_size <-
  function(in_dir){
    le <- "\\"
    glue("

gzet -dc {in_dir}/*.gz$ | wc -c

       ")
  }
