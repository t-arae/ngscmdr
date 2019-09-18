
#' Create new cmdopt class object
#' @param ... ...
#' @param cmd_id command type specifier
new_cmdopt <-
  function(..., cmd_id = NULL){
    cmdopt_obj <- list(...)
    class(cmdopt_obj) <- c(cmd_id, "cmdopt", class(cmdopt_obj))
    cmdopt_obj
  }

#' Print cmdopt class object
#' @param x cmdopt class object
#' @param ... ...
#' @export
print.cmdopt <-
  function(x, ...){
    cat("cmdopt class object\n")
    cat("input file: ", paste(as.character(x$inf), collapse = " "), "\n")
    cat("output file directory: ", x$out_dir, "\n")
    cat("output file name: ", x$outf, "\n")
  }

#' cmd_get methods for cmdopt class object
#' @param x cmdopt class object
#' @export
cmd_get <- function(x) UseMethod("cmd_get")

cmd_set_star <-
  function(
    fastq_obj,
    out_dir = "./mapped_by_star",
    fasta_path = "./allChr.fas",
    idx_dir = "./idx_star",
    gff_path = "./TAIR10_GFF3_genes_transposons.gff",
    core_num = 4
  ){
    new_cmdopt(
      fastq_obj = fastq_obj,
      out_dir = out_dir,
      fasta_path = fasta_path,
      idx_dir = idx_dir,
      gff_path = gff_path,
      core_num = core_num,
      cmd_id = "star"
    )
  }

cmd_get.star <-
  function(x){
    x$le <- "\\"
    x$fq_ext <- get_file_ext()$fastq
    x$fasta_path <- paste(x$fasta_path, collapse = " ")
    x$inf <- paste(x$fastq_obj, collapse = " ")
    genome_indexing_string <- "
mkdir {idx_dir}
STAR {le}
  --runMode genomeGenerate {le}
  --runThreadN {core_num} {le}
  --genomeDir {idx_dir} {le}
  --genomeFastaFiles {fasta_path} {le}
  --sjdbGTFfile {gff_path} {le}
  --sjdbOverhang 35 {le}
  --sjdbGTFtagExonParentTranscript Parent
"
    mapping_string <- "
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
  --readFilesIn {inf} {le}
  --outFileNamePrefix {label}

mkdir {out_dir}
mv {label}Aligned.sortedByCoord.out.bam {out_dir}/{label}.sort.bam
mv {label}Log.final.out {out_dir}/{label}Log.final.out
mv {label}Log.out {out_dir}/{label}Log.out
rm -R {label}_STARgenome
rm {label}Log.progress.out
rm {label}SJ.out.tab
"

    glue(paste0(genome_indexing_string, mapping_string), .envir = as.environment(x))
  }

# cmd_set_star(
#   new_fastq(c("hoge"), "./fastq/hoge.gz"),
#   fasta_path = paste0("./fasta/Chr", c(1:5, "M", "C"), ".fas"),
#   gff_path = "./TAIR10_GFF3_genes.gff"
# ) %>%
#   cmd_get()
