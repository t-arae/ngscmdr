

#' Generate genome index for STAR
#' @importFrom glue glue
#' @param stargenome_dir .tar file path
#' @param chr_fasta_path .tar file path
#' @param gff_path .tar file path
#' @param core_num .tar file path
#' @export
glue_star_genome_generate <-
  function(
    stargenome_dir = "./stargenome",
    chr_fasta_path = "./allChr.fas",
    gff_path = "./TAIR10_GFF3_genes_transposons.gff",
    core_num = 2
  ){
    lineend <- "\\"
    glue(
      "

mkdir {stargenome_dir}
STAR {lineend}
  --runThreadN {core_num} {lineend}
  --runMode genomeGenerate {lineend}
  --genomeDir {stargenome_dir} {lineend}
  --genomeFastaFiles {chr_fasta_path} {lineend}
  --sjdbGTFfile {gff_path} {lineend}
  --sjdbOverhang 35 {lineend}
  --sjdbGTFtagExonParentTranscript Parent

    "
    )
  }

#' Map with STAR and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @param stargenome_dir .tar file path
#' @param gff_path .tar file path
#' @export
glue_star_bamsort <-
  function(head_label, in_dir = "./fastq", out_dir = "./mapped_by_star",
           stargenome_dir = "./stargenome",
           gff_path = "TAIR10_GFF3_genes_transposons.gff"){
    lineend <- "\\"
    glue(
      "

STAR {lineend}
  --runThreadN 2 {lineend}
  --genomeDir {stargenome_dir} {lineend}
  --sjdbGTFfile {gff_path} {lineend}
  --sjdbGTFtagExonParentTranscript Parent {lineend}
  --alignIntronMin 20 {lineend}
  --alignIntronMax 3000 {lineend}
  --outSAMstrandField intronMotif {lineend}
  --outSAMmultNmax 1 {lineend}
  --outMultimapperOrder Random {lineend}
  --outFilterMultimapNmax 20 {lineend}
  --outFilterMismatchNmax 3 {lineend}
  --sjdbOverhang 35 {lineend}
  --seedSearchStartLmax 15 {lineend}
  --outSAMtype BAM SortedByCoordinate {lineend}
  --readFilesIn {in_dir}/{head_label}.fastq {lineend}
  --outFileNamePrefix {head_label}

mkdir {out_dir}/{head_label}
mv {head_label}Aligned.sortedByCoord.out.bam {out_dir}/{head_label}/{head_label}.sort.bam
mv {head_label}Log.final.out {out_dir}/{head_label}/{head_label}Log.final.out
rm -R {head_label}_STARgenome
rm {head_label}Log.progress.out
rm {head_label}Log.out
rm {head_label}SJ.out.tab

    "
    )
  }


#' Pipeline: STAR -> Cufflinks
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param stargenome_dir .tar file path
#' @param gff_path .tar file path
#' @export
glue_star_pipeline <-
  function(
    head_label,
    in_dir = "./fastq",
    stargenome_dir = "./stargenome",
    gff_path = "TAIR10_GFF3_genes_transposons.gff"
  ){
    fastq_file_dir <- in_dir
    star_output_dir <- "./mapped_by_star"

    temp <-
      glue_star_bamsort(
        head_label = head_label,
        in_dir = fastq_file_dir,
        out_dir = star_output_dir,
        stargenome_dir = stargenome_dir,
        gff_path = gff_path
      )

    bam_file_dir <- paste(star_output_dir, head_label, sep = "/")
    cufflinks_out_dir <- "./fpkm"
    temp2 <-
      glue_cufflinks(
        head_label = head_label,
        in_dir = bam_file_dir,
        out_dir = cufflinks_out_dir,
        gff_path = gff_path
      )

    cat(temp, temp2)
  }
