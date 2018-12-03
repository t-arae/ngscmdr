if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))


#' Extract files from arachive file (.tar)
#' @importFrom glue glue
#' @param archive_file .tar file path.
#' @export
glue_extract_file_from_tar <-
  function(archive_file){
    extract_to <- dirname(archive_file)
    glue("

tar xvf {archive_file} -C {extract_to}

      ")
  }


#' Merge some fastq files to a fastq file
#' @importFrom glue glue
#' @param fpath_in_fastq_group .tar file path
#' @param fpath_out_fastq .tar file path
#' @export
glue_merge_fastq <-
  function(fpath_in_fastq_group, fpath_out_fastq){
    glue("

cat {fpath_in_fastq_group} | gzip -c > {fpath_out_fastq}
rm {fpath_in_fastq_group}

    ")
  }

#' Merge some fastq files to a fastq file
#' @importFrom glue glue
#' @param gziped_files .tar file path
#' @param gzip_out .tar file path
#' @export
glue_merge_gz <-
  function(gziped_files, gzip_out){
    glue("

cat {gziped_files} > {gzip_out}
rm {gziped_files}

    ")
  }


#' Create fastqc reports
#' @importFrom glue glue
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_fastqc <-
  function(in_dir = "./fastq", out_dir = "./fastqc_out"){
    glue("

mkdir {out_dir}
fastqc -o {out_dir} {in_dir}/*.fastq.gz

        ")
  }


#' Calculate genome coverage
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_genomecov <-
  function(head_label, in_dir, out_dir = in_dir){
    le <- "\\"
    glue("

mkdir {out_dir}/genomecov
bedtools genomecov -bga -split {le}
  -ibam {in_dir}/{head_label}.sort.bam > {out_dir}/genomecov/{head_label}.bedgraph

      ")
  }

#' Convert gff3 file to bed file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_gff2bed <-
  function(head_label, in_dir = ".", out_dir = "."){
    glue("

gff2bed < {in_dir}/{head_label}.gff > {out_dir}/{head_label}.bed

      ")
  }

glue_samtools_view <-
  function(head_label, bedfile, out_label, in_dir = ".", out_dir = "."){
    lineend <- "\\"
    glue(
      "
samtools view -@ 2 -b -L {bedfile} {in_dir}/{head_label}.sort.bam {lineend}
| {lineend}
samtools sort -@ 2 > {out_dir}/{head_label}_{out_label}.sort.bam
samtools index {out_dir}/{head_label}_{out_label}.sort.bam
      "
    )
  }


glue_rseqc_clipping_profile <-
  function(
    head_label,
    in_dir,
    out_dir = in_dir,
    seq_layout = '"SE"'
  ){
    le <- "\\"
    seq_layout <- '"SE"'
    glue("

mkdir {out_dir}/clip_prof
clipping_profile.py {le}
  -i {in_dir}/{head_label}.sort.bam {le}
  -s {seq_layout} {le}
  -o {out_dir}/clip_prof/{head_label}

       ")
  }

glue_rseqc_deletion_profile <-
  function(
    head_label,
    in_dir,
    out_dir = in_dir
  ){
    le <- "\\"
    seq_layout <- '"SE"'
    glue("

mkdir {out_dir}/clip_prof
clipping_profile.py {le}
  -i {in_dir}/{head_label}.sort.bam {le}
  -s {seq_layout} {le}
  -o {out_dir}/clip_prof/{head_label}

       ")
  }

glue_multiqc <-
  function(
    out_dir = "."
  ){
    glue("
multiqc -f -o {out_dir} .
         ")
  }

glue_sra2fastq <-
  function(
    head_label,
    is_pe = F,
    in_dir = "./sra",
    out_dir = "./fastq"
  ){
    le <- "\\"
    pe_option <- ifelse(is_pe, "\n  --split-files \\", "")
glue("

fastq-dump {le}{pe_option}
  -O {out_dir} {le}
  --gzip {le}
  {in_dir}/{head_label}.sra

     ")
  }
