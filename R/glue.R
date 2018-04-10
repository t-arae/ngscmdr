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
    glue(
      "
cat {fpath_in_fastq_group} > {fpath_out_fastq}
    "
    )
  }


#' Create fastqc reports
#' @importFrom stringr str_detect
#' @importFrom magrittr %>%
#' @importFrom glue glue
#' @param jissai_in_dir .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_fastqc <-
  function(jissai_in_dir, in_dir = "./fastq", out_dir = "./fastqc_out"){
    all_files <- list.files(jissai_in_dir)
    fastq_path <-
      str_detect(all_files, ".*(.fastq$)|(.fastq.gz$)") %>%
      all_files[.] %>%
      paste(in_dir, ., sep = "/") %>%
      paste(collapse = " ")
    glue("
fastqc -o {out_dir} {fastq_path}
           ")
  }

#' Generate genome index for bowtie
#' @importFrom glue glue
#' @param fasta_path .tar file path
#' @param index_name .tar file path
#' @export
glue_bowtie_genome_generate <-
  function(fasta_path, index_name){
    glue("
bowtie-build -f {fasta_path} {index_name}
         ")
  }

#' Generate genome index for bowtie2
#' @importFrom glue glue
#' @param fasta_path .tar file path
#' @param index_name .tar file path
#' @export
glue_bowtie2_genome_generate <-
  function(fasta_path, index_name){
    glue("
bowtie2-build -f {fasta_path} {index_name}
         ")
  }


#' Map with bowtie2 and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param index_name .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_bowtie2_bamsort <-
  function(
    head_label,
    index_name,
    in_dir = "./fastq",
    out_dir = "./mapped_by_bowtie2"
  ){
    line_end <- "\\"
    glue(
      "
mkdir {out_dir}
bowtie2 -p 2 -a --best --strata --trim5 8 --trim3 40 {line_end}
  -x {index_name} -U {in_dir}/{head_label}.fastq -S {in_dir}/{head_label}.sam && {line_end}
mv {in_dir}/{head_label}.sam {out_dir} && {line_end}
samtools sort -O bam -o {out_dir}/{head_label}_bowtie2.sort.bam {line_end}
  {out_dir}/{head_label}.sam && {line_end}
samtools index {out_dir}/{head_label}_bowtie2.sort.bam && {line_end}
rm {out_dir}/{head_label}.sam

      "
    )
  }


#' Map with bowtie and output a sorted bam file
#' @importFrom glue glue
#' @param head_label .tar file path
#' @param index_name .tar file path
#' @param in_dir .tar file path
#' @param out_dir .tar file path
#' @export
glue_bowtie_bamsort <-
  function(
    head_label,
    index_name,
    in_dir = "./fastq",
    out_dir = "./mapped_by_bowtie"
  ){
    line_end <- "\\"
    glue(
      "
mkdir {out_dir}
bowtie -p 2 -S -a --best --strata {line_end}
  {index_name} {in_dir}/{head_label}.fastq {in_dir}/{head_label}.sam && {line_end}
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
