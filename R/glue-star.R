

#' Generate genome index for STAR
#' @importFrom glue glue
#' @inheritParams param_general
#' @inheritParams param_genomeidx
#' @param gff_path .tar file path
#' @export
glue_star_genome_generate <-
  function(
    idx_dir = "./idx_star",
    fasta_path = "./allChr.fas",
    gff_path = "./TAIR10_GFF3_genes_transposons.gff",
    core_num = 4
  ){
    le <- "\\"
    glue("

mkdir {idx_dir}
STAR {le}
  --runMode genomeGenerate {le}
  --runThreadN {core_num} {le}
  --genomeDir {idx_dir} {le}
  --genomeFastaFiles {fasta_path} {le}
  --sjdbGTFfile {gff_path} {le}
  --sjdbOverhang 35 {le}
  --sjdbGTFtagExonParentTranscript Parent

    ")
  }

#' Map with STAR and output a sorted bam file
#' @importFrom glue glue
#' @inheritParams param_general
#' @inheritParams param_genomeidx
#' @param gff_path .tar file path
#' @export
glue_se_star_bamsort <-
  function(
    head_label,
    in_dir = "./fastq",
    out_dir = "./mapped_by_star",
    idx_dir = "./idx_star",
    gff_path = "TAIR10_GFF3_genes_transposons.gff",
    core_num = 4
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

STAR {le}
  --runThreadN {core_num} {le}
  --genomeDir {idx_dir} {le}
  --sjdbGTFfile {gff_path} {le}
  --sjdbGTFtagExonParentTranscript Parent {le}
  --alignIntronMin 20 {le}
  --alignIntronMax 3000 {le}
  --outSAMstrandField intronMotif {le}
  --outSAMmultNmax 1 {le}
  --outMultimapperOrder Random {le}
  --outFilterMultimapNmax 20 {le}
  --outFilterMismatchNmax 3 {le}
  --sjdbOverhang 35 {le}
  --seedSearchStartLmax 15 {le}
  --outSAMtype BAM SortedByCoordinate {le}
  --readFilesCommand zcat {le}
  --readFilesIn {in_dir}/{head_label}.{fq_ext}.gz {le}
  --outFileNamePrefix {head_label}

mkdir {out_dir}
mv {head_label}Aligned.sortedByCoord.out.bam {out_dir}/{head_label}.sort.bam
mv {head_label}Log.final.out {out_dir}/{head_label}Log.final.out
mv {head_label}Log.out {out_dir}/{head_label}Log.out
rm -R {head_label}_STARgenome
rm {head_label}Log.progress.out
rm {head_label}SJ.out.tab

    ")
  }

#' Map with STAR and output a sorted bam file
#' @importFrom glue glue
#' @inheritParams param_general
#' @inheritParams param_genomeidx
#' @param gff_path .tar file path
#' @export
glue_pe_star_bamsort <-
  function(
    head_label,
    in_dir = "./fastq",
    out_dir = "./mapped_by_star",
    idx_dir = "./idx_star",
    gff_path = "TAIR10_GFF3_genes_transposons.gff",
    core_num = 4
  ){
    le <- "\\"
    fq_ext <- fq_ext
    glue("

STAR {le}
  --runThreadN {core_num} {le}
  --genomeDir {idx_dir} {le}
  --sjdbGTFfile {gff_path} {le}
  --sjdbGTFtagExonParentTranscript Parent {le}
  --alignIntronMin 20 {le}
  --alignIntronMax 3000 {le}
  --outSAMstrandField intronMotif {le}
  --outSAMmultNmax 1 {le}
  --outMultimapperOrder Random {le}
  --outFilterMultimapNmax 20 {le}
  --outFilterMismatchNmax 3 {le}
  --sjdbOverhang 35 {le}
  --seedSearchStartLmax 15 {le}
  --outSAMtype BAM SortedByCoordinate {le}
  --readFilesCommand zcat {le}
  --readFilesIn {in_dir}/{head_label}_1.{fq_ext}.gz {in_dir}/{head_label}_2.{fq_ext}.gz{le}
  --outFileNamePrefix {head_label}

mkdir {out_dir}
mv {head_label}Aligned.sortedByCoord.out.bam {out_dir}/{head_label}.sort.bam
mv {head_label}Log.final.out {out_dir}/{head_label}Log.final.out
mv {head_label}Log.out {out_dir}/{head_label}Log.out
rm -R {head_label}_STARgenome
rm {head_label}Log.progress.out
rm {head_label}SJ.out.tab

    ")
  }
