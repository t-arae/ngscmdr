if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))


#' Extract files from arachive file (.tar)
#' @importFrom glue glue
#' @param archive_file .tar file path.
#' @export
glue_extract_file_from_tar <-
  function(archive_file){
    extract_to <- dirname(archive_file)
    glue(
      "
tar xvf {archive_file} -C {extract_to}
      "
    )
  }


#' Merge some fastq files to a fastq file
#' @importFrom glue glue
#' @param fpath_in_fastq_group .tar file path
#' @param fpath_out_fastq .tar file path
#' @export
glue_merge_fastq <-
  function(fpath_in_fastq_group, fpath_out_fastq){
    glue("
cat {fpath_in_fastq_group} > {fpath_out_fastq}
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
  function(head_label, in_dir = ".", out_dir = "."){
    glue(
      "
bedtools genomecov -bga -split -ibam {in_dir}/{head_label}.sort.bam > {out_dir}/{head_label}.bedgraph
      "
    )
  }

glue_gff2bed <-
  function(head_label, in_dir = ".", out_dir = "."){
    glue(
      "
gff2bed < {in_dir}/{head_label}.gff > {out_dir}/{head_label}.bed
      "
    )
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

#' Calculate FPKM by Cufflinks
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @param gff_path .tar file path
#' @export
glue_cufflinks <-
  function(head_label, in_dir, out_dir,
           gff_path = "TAIR10_GFF3_genes_transposons.gff"
  ){
    line_end <- "\\"
    glue(
      "

mkdir {out_dir}
cufflinks {line_end}
  -o {out_dir} {line_end}
  -p 2 {line_end}
  -G {gff_path} {line_end}
  {in_dir}/{head_label}.sort.bam

mv {out_dir}/genes.fpkm_tracking {out_dir}/{head_label}_genes.fpkm_tracking
mv {out_dir}/isoforms.fpkm_tracking {out_dir}/{head_label}_isoforms.fpkm_tracking

    "
    )
  }


#' featureCounts
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @param gtf_path .tar file path
#' @export
glue_featurecounts <-
  function(
    head_label,
    in_dir,
    out_dir = "./readcount",
    gtf_path = "./TAIR10_GFF3_genes_transposons.gtf"
  ){
    line_end <- "\\"
    glue("
mkdir {out_dir}
featureCounts {line_end}
  -s 1 {line_end}
  -t exon {line_end}
  -g gene_id {line_end}
  -a {gtf_path} {line_end}
  -o {out_dir}/{head_label}_counts.txt {line_end}
  {in_dir}/{head_label}.sort.bam
         ")
  }
