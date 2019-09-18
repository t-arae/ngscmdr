
#' Generate genome index for bowtie2
#' @importFrom glue glue
#' @param fasta_path .tar file path
#' @param idx_dir .tar file path
#' @export
glue_bowtie2_genome_generate <-
  function(fasta_path, idx_dir = "./idx_bw2"){
    bn <- basename(fasta_path)
    glue("

mkdir {idx_dir}
bowtie2-build -f {fasta_path} {idx_dir}/{bn}

         ")
  }

#' Remove reads from contamination with bowtie2
#' @importFrom glue glue
#' @inheritParams param_general
#' @inheritParams param_genomeidx
#' @param fasta_path .tar file path
#' @export
glue_bowtie2_rm_contami <-
  function(
    head_label,
    in_dir,
    fasta_path,
    out_dir = "fastq_rm",
    idx_dir = "./idx_bw2",
    core_num = 4
  ){
    le <- "\\"
    fq_ext <- fq_ext
    bn <- basename(fasta_path)
    glue("

mkdir {out_dir}
bowtie2 {le}
  -L 20 {le}
  -p {core_num} {le}
  -t {le}
  --quiet {le}
  --no-unal {le}
  -x {idx_dir}/{bn} {le}
  -U {in_dir}/{head_label}.{fq_ext}.gz {le}
  --un-gz {out_dir}/{head_label}.{fq_ext}.gz {le}
  --al-gz {out_dir}/{head_label}_al.{fq_ext}.gz {le}
  > /dev/null

    ")
  }

#' Remove paired end reads from contamination with bowtie2
#' @importFrom glue glue
#' @inheritParams param_general
#' @inheritParams param_genomeidx
#' @param fasta_path .tar file path
#' @export
glue_bowtie2_rm_contami_pe <-
  function(
    head_label,
    in_dir,
    fasta_path,
    out_dir = "fastq_rm",
    idx_dir = "./idx_bw2",
    core_num = 4
  ){
    le <- "\\"
    fq_ext <- fq_ext
    bn <- basename(fasta_path)
    glue("

mkdir {out_dir}
bowtie2 {le}
  -L 20 {le}
  -p {core_num} {le}
  -t {le}
  --quiet {le}
  --no-unal {le}
  -x {idx_dir}/{bn} {le}
  -1 {in_dir}/{head_label}_1.{fq_ext}.gz {le}
  -2 {in_dir}/{head_label}_2.{fq_ext}.gz {le}
  --un-conc-gz {out_dir}/{head_label}.{fq_ext}.gz {le}
  --al-conc-gz {out_dir}/{head_label}_al.{fq_ext}.gz {le}
  > /dev/null

    ")

  }
# bowtie2 {le}
#   -L 20 {le} # length of seed substrings; must be >3, <32 (22)
#   -p {core_num} {le}
#   -t {le} # print wall-clock time taken by search phases
#   --quiet {le} # print nothing to stderr except serious errors
#   -x {idx_dir} {le}
#   -U {in_dir}/{head_label}.fastq.gz {le}
#   --un-gz {out_dir}/{head_label}.fastq.gz

#' Map with bowtie2 and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param idx_dir .tar file path
#' @param core_num .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_bowtie2_bamsort <-
  function(
    head_label,
    idx_dir = "./idx_bw2",
    core_num = 2,
    in_dir = "./fastq",
    out_dir = "./mapped_by_bowtie2"
  ){
    lineend <- "\\"
    fq_ext <- fq_ext
    glue(
      "
mkdir {out_dir}
bowtie2 {lineend}
  -p {core_num} {lineend}
  -a --best --strata {lineend}
  --trim5 8 {lineend}
  --trim3 40 {lineend}
  -x {index_name} {lineend}
  -U {in_dir}/{head_label}.{fq_ext} {lineend}
| {lineend}
samtools view -@ {core_num} -bS {lineend}
| {lineend}
samtools sort -@ {core_num} > {out_dir}/{head_label}_bowtie2.sort.bam
samtools index {out_dir}/{head_label}_bowtie2.sort.bam && {line_end}

      "
    )
  }
