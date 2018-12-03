
#' Generate genome index for HISAT2
#' @importFrom glue glue
#' @inheritParams param_general
#' @inheritParams param_genomeidx
#' @export
glue_hisat_genome_generate <-
  function(
    fasta_path = "./allChr.fas",
    idx_name = "TAIR10",
    idx_dir = "./idx_hisat",
    core_num = 2
  ){
    le <- "\\"
    glue(
      "

mkdir {idx_dir}
hisat2-build {le}
  -p {core_num} {le}
  --seed 1 {le}
  {fasta_path} {le}
  {idx_dir}/{idx_name}

    "
    )
  }

#' Map single-end reads with hisat and output a sorted bam file
#' @importFrom glue glue
#' @inheritParams param_general
#' @inheritParams param_genomeidx
#' @export
glue_se_hisat_bamsort <-
  function(
    head_label,
    core_num = 2,
    idx_name = "TAIR10",
    idx_dir = "./idx_hisat",
    in_dir = "./fastq",
    out_dir = "./mapped_by_hisat"
  ){
    le <- "\\"
    glue(
      "

mkdir {out_dir}
hisat2 {le}
  --dta-cufflinks {le}
  --new-summary {le}
  --summary-file {out_dir}/report_{head_label}.txt {le}
  -p {core_num} {le}
  -x {idx_dir}/{idx_name} {le}
  -U {in_dir}/{head_label}.fastq.gz {le}
| {le}
samtools view -@ {core_num} -bS {le}
| {le}
samtools sort -@ {core_num} > {out_dir}/{head_label}.sort.bam
samtools index {out_dir}/{head_label}.sort.bam

    "
    )
  }


#' Map paired-end reads with hisat and output a sorted bam file
#' @importFrom glue glue
#' @inheritParams param_general
#' @inheritParams param_genomeidx
#' @export
glue_pe_hisat_bamsort <-
  function(head_label,
           core_num = 2,
           idx_name = "TAIR10",
           idx_dir = "./idx_hisat",
           in_dir = "./fastq",
           out_dir = "./mapped_by_hisat"
  ){
    le <- "\\"
    glue(
      "

mkdir {out_dir}
hisat2 {le}
  --dta-cufflinks {le}
  --new-summary {le}
  --summary-file {out_dir}/report_{head_label}.txt {le}
  -p {core_num} {le}
  -x {idx_dir}/{idx_name} {le}
  -1 {in_dir}/{head_label}_1.fastq.gz {le}
  -2 {in_dir}/{head_label}_2.fastq.gz {le}
| {le}
samtools view -@ {core_num} -bS {le}
| {le}
samtools sort -@ {core_num} > {out_dir}/{head_label}.sort.bam
samtools index {out_dir}/{head_label}.sort.bam

    "
    )
  }
