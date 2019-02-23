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
    cat("input file directory: ", x$in_dir, "\n")
    cat("input file: ", x$inf, "\n")
    cat("output file directory: ", x$out_dir, "\n")
    cat("output file name: ", x$outf, "\n")
  }

#' cmd_get methods for cmdopt class object
#' @param x cmdopt class object
#' @export
cmd_get <- function(x) UseMethod("cmd_get")

cmd_set_star <-
  function(
    sample_label,
    in_dir = "./fastq",
    out_dir = "./mapped_by_star",
    fasta_path = "./allChr.fas",
    idx_dir = "./idx_star",
    gff_path = "./TAIR10_GFF3_genes_transposons.gff",
    core_num = 4
  ){
    new_cmdopt(
      label = sample_label,
      in_dir = in_dir,
      out_dir = out_dir,
      fasta_path = "./allChr.fas",
      idx_dir = "./idx_star",
      gff_path = "./TAIR10_GFF3_genes_transposons.gff",
      core_num = 4,
      cmd_id = "star"
    )
  }

cmd_get.star <-
  function(x){
    x$le <- "\\"
    x$fq_ext <- get_file_ext()$"fastq"
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
  --readFilesIn {in_dir}/{label}.{fq_ext}.gz {le}
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
