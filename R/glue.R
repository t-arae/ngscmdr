
#' Create fastqc reports
#' @importFrom glue glue
#' @inheritParams param_general
#' @export
glue_fastqc <-
  function(in_dir = "./fastq", out_dir = "./fastqc_out", core_num = 4){
    fq_ext <- get_file_ext()$fastq
    glue("

mkdir {out_dir}
fastqc --extract -t {core_num} -o {out_dir} {in_dir}/*.{fq_ext}.gz

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
