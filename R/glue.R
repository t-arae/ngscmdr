
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
    glue(
      "
cat {fpath_in_fastq_group} > {fpath_out_fastq}
    "
    )
  }

#' Map with bowtie and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_bowtie_bamsort <-
  function(head_label, in_dir = "./fastq", out_dir = "./mapped_by_bowtie"){
    line_end <- "\\"
    glue(
      "
mkdir {out_dir}
bowtie -p 2 -S {line_end}
  TAIRIDX {in_dir}/{head_label}.fastq {in_dir}/{head_label}.sam && {line_end}
mv {in_dir}/{head_label}.sam {out_dir} && {line_end}
samtools sort -O bam -o {out_dir}/{head_label}_bowtie.sort.bam {line_end}
  {out_dir}/{head_label}.sam && {line_end}
samtools index {out_dir}/{head_label}_bowtie.sort.bam && {line_end}
rm {out_dir}/{head_label}.sam

      "
    )
  }


#' Generate genome index for STAR
#' @importFrom glue glue
#' @param stargenome_dir .tar file path
#' @param chr_fasta_path .tar file path
#' @param gff_path .tar file path
#' @export
glue_star_genome_generate <-
  function(stargenome_dir = "./stargenome",
           chr_fasta_path = "./allChr.fas",
           gff_path = "./TAIR10_GFF3_genes_transposons.gff"
  ){
    glue(
      "
mkdir {stargenome_dir}
STAR \
--runThreadN 2 \
--runMode genomeGenerate \
--genomeDir {stargenome_dir} \
--genomeFastaFiles {chr_fasta_path} \
--sjdbGTFfile {gff_path} \
--sjdbOverhang 35 \
--sjdbGTFtagExonParentTranscript Parent
    "
    )
  }


#' Map with bowtie and output a sorted bam file
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
    line_end <- "\\"
    glue(
      "

STAR {line_end}
  --runThreadN 2 {line_end}
  --genomeDir {stargenome_dir} {line_end}
  --sjdbGTFfile {gff_path} {line_end}
  --sjdbGTFtagExonParentTranscript Parent {line_end}
  --alignIntronMin 20 {line_end}
  --alignIntronMax 3000 {line_end}
  --outSAMstrandField intronMotif {line_end}
  --outSAMmultNmax 1 {line_end}
  --outMultimapperOrder Random {line_end}
  --outFilterMultimapNmax 20 {line_end}
  --outFilterMismatchNmax 3 {line_end}
  --sjdbOverhang 35 {line_end}
  --seedSearchStartLmax 15 {line_end}
  --outSAMtype BAM SortedByCoordinate {line_end}
  --readFilesIn {in_dir}/{head_label}.fastq {line_end}
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

#' Calculate genome coverage
#' @importFrom glue glue
#' @param fpath_sorted_bam .tar file path
#' @param fpath_out .tar file path
#' @export
glue_genomecov <-
  function(fpath_sorted_bam, fpath_out){
    glue(
      "
bedtools genomecov -bga -split -ibam {fpath_sorted_bam} > {fpath_out}
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


#' Calculate FPKM by Cufflinks
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
